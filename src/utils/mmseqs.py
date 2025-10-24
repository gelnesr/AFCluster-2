import os
import time
import tqdm
import random
import shutil
import pathlib
import logging
import requests
import tarfile
import tempfile
import subprocess
from typing import Tuple, List
import hashlib


logger = logging.getLogger(__name__)
def _have_local_mmseqs():
    mmseqs_bin = os.environ.get("MMSEQS_BIN") or shutil.which("mmseqs")
    return mmseqs_bin is not None

def _mmseqs_bin():
    return os.environ.get("MMSEQS_BIN") or shutil.which("mmseqs") or "mmseqs"

def _db_prefixes(use_env: bool):
    dbs = []
    uniref = os.environ.get("MMSEQS_DB_UNIREF")
    envdb  = os.environ.get("MMSEQS_DB_ENV")
    if uniref:
        dbs.append(uniref)
    if use_env and envdb:
        dbs.append(envdb)
    return dbs

def _write_query_fasta(seqs: List[str], start_n: int, fasta_path: str):
    with open(fasta_path, "w") as f:
        n = start_n
        for s in seqs:
            f.write(f">{n}\n{s}\n")
            n += 1

def _run_one_local_msa(query_fa: str, db_prefix: str, out_a3m: str, threads: int = 4):
    """
    Minimal local pipeline:
      createdb -> search -> result2msa (A3M)
    Assumes 'db_prefix' points to an existing mmseqs2 DB (prefixed files).
    """
    mm = _mmseqs_bin()
    tmpd = tempfile.mkdtemp(prefix="mmseqs_tmp_")
    try:
        qdb = os.path.join(tmpd, "qdb")
        res = os.path.join(tmpd, "res")

        subprocess.run([mm, "createdb", query_fa, qdb], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([mm, "search", qdb, db_prefix, res, tmpd,
                        "--threads", str(threads),
                        "-s", "7.5"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([mm, "result2msa", db_prefix, qdb, res, out_a3m,
                        "--msa-format-mode", "6",
                        "--threads", str(threads)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)

def _gather_a3m_blocks(a3m_path: str) -> dict:
    """Return {header_int: [lines]} mapping like your remote parser."""
    blocks = {}
    if not os.path.isfile(a3m_path):
        return blocks
    update_M, M = True, None
    with open(a3m_path) as fh:
        for line in fh:
            if not line:
                continue
            if "\x00" in line:
                line = line.replace("\x00","")
                update_M = True
            if line.startswith(">") and update_M:
                try:
                    M = int(line[1:].strip())
                except ValueError:
                    M = None
                update_M = False
                if M is not None and M not in blocks:
                    blocks[M] = []
            if M is not None:
                blocks[M].append(line)
    return blocks

def _run_mmseqs2_local(x, prefix, use_env=True, **kwargs):
    """
    Local implementation:
      - de-duplicate queries (preserve your >N mapping)
      - search against available DBs (UNIREf and optionally ENV)
      - concatenate A3Ms per query (UNIREf first, then ENV)
    Returns list[str] like the remote path.
    """

    seqs = [x] if isinstance(x, str) else x

    seqs_unique = []
    for s in seqs:
        if s not in seqs_unique:
            seqs_unique.append(s)

    N = 101
    Ms = [N + seqs_unique.index(seq) for seq in seqs]

    dbs = _db_prefixes(use_env=use_env)
    if not dbs:
        raise RuntimeError("Local mmseqs2 found, but no databases configured. "
                           "Set MMSEQS_DB_UNIREF (and optionally MMSEQS_DB_ENV).")

    work = pathlib.Path(f"{prefix}_local")
    work.mkdir(parents=True, exist_ok=True)
    qfa = str(work / "queries.fasta")
    _write_query_fasta(seqs_unique, N, qfa)

    per_db_blocks = []
    for i, db in enumerate(dbs):
        out_a3m = str(work / f"out_{i}.a3m")
        _run_one_local_msa(qfa, db, out_a3m, threads=int(os.environ.get("MMSEQS_THREADS", "4")))
        per_db_blocks.append(_gather_a3m_blocks(out_a3m))

    a3m_lines = []
    for m in Ms:
        merged = []
        for blocks in per_db_blocks:
            if m in blocks:
                merged.extend(blocks[m])
        a3m_lines.append("".join(merged) if merged else "")

    return a3m_lines

def get_hash(x):
    return hashlib.sha1(x.encode()).hexdigest()

def run_mmseqs2(x, prefix, use_env=True, use_filter=True,
                filter=None, pairing_strategy="greedy",
                host_url="https://api.colabfold.com", user_agent="") :
    '''
    Generate MSA by running MMSeqs-2 via API
    '''
    
    submission_endpoint = "ticket/msa"
    headers = {}
    if user_agent != "":
        headers['User-Agent'] = user_agent
    #else:
    #    logger.warning("No user agent specified. Please set a user agent (e.g., 'toolname/version contact@email') to help us debug in case of problems. This warning will become an error in the future.")

    def submit(seqs, mode, N=101):
        nx = N
        query = ""
        for seq in seqs:
            query += f'>{nx}\n{seq}\n'
            nx += 1

        while True:
            error_count = 0
            try:
                res = requests.post(f'{host_url}/{submission_endpoint}', data={'q': query, 'mode': mode}, timeout=6.02, headers=headers)
            except requests.exceptions.Timeout:
                logger.warning("Timeout while submitting to MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break

        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        while True:
            error_count = 0
            try:
                res = requests.get(f'{host_url}/ticket/{ID}', timeout=6.02, headers=headers)
            except requests.exceptions.Timeout:
                logger.warning("Timeout while fetching status from MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        error_count = 0
        while True:
            try:
                res = requests.get(f'{host_url}/result/download/{ID}', timeout=6.02, headers=headers)
            except requests.exceptions.Timeout:
                logger.warning("Timeout while fetching result from MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        with open(path, "wb") as out: out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # compatibility to old option
    if filter is not None:
        use_filter = filter

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    path = f"{prefix}_{mode}"
    if not os.path.isdir(path): os.mkdir(path)

    tar_gz_file = f'{path}/out.tar.gz'
    N, REDO = 101, True

    seqs_unique = []

    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]

    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique)
        with tqdm.tqdm(total=TIME_ESTIMATE) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

                if out["status"] == "MAINTENANCE":
                    raise Exception(f'MMseqs2 API is undergoing maintenance. Please try again in a few minutes.')

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 2)
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)
                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE-TIME))
                    REDO = False
                if out["status"] == "ERROR":
                    REDO = False
                    raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

            # Download results
            download(ID, tar_gz_file)

    a3m_files = [f"{path}/uniref.a3m"]
    if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    
    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines: a3m_lines[M] = []
                a3m_lines[M].append(line)

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]

    return a3m_lines


def run_mmseqs(seq, temp_dir='./tmp', use_env=True, use_filter=True, filter=None,
               host_url="https://api.colabfold.com"):

    def get_hash(x):
        return hashlib.sha1(x.encode()).hexdigest()

    prefix = get_hash(seq)
    prefix = os.path.join(temp_dir, prefix)
    
    if _have_local_mmseqs():
        try:
            return _run_mmseqs2_local(seq, prefix, use_env=use_env)
        except Exception as e:
            logger.warning(f"Local mmseqs2 failed ({e}); falling back to ColabFold API.")

    return run_mmseqs2(seq, prefix, 
                        use_env=use_env, 
                        use_filter=use_filter, 
                        filter=filter,   
                        host_url=host_url)
