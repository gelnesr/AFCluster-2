import os
import sys
import glob
import copy
import yaml
import argparse
import subprocess
import numpy as np
import pandas as pd

from src.cluster import *
from src.utils.msa import *
from src.utils.seqs import *
from src.utils.mmseqs import *
from src.utils.helpers import *

def get_labels(args, df):
    if args.cluster_method == "dbscan":
        return cluster_DBSCAN(args, df)
    raise ValueError(f"Unknown clustering method")


def generate_command(args, return_seeds=False):
    run_command = ['colabfold_batch']
    afcluster = dict_to_namespace(args.afcluster)
    if return_seeds:
        return afcluster.num_seeds
    if afcluster.amber:
        run_commmand.extend(['--amber'])
        if torch.cuda.is_available():
            run_command.extend(['--use-gpu-relax'])
        run_command.extend(['--num-relax', f'{afcluster.num_relax}'])
    if afcluster.use_dropout:
        run_command.extend(['--use-dropout'])
    run_command.extend(['--num-recycle', f'{afcluster.num_recycle}'])
    if afcluster.templates:
        run_command.extend(['--templates'])
    return run_command


def run_cluster(args, subfolder, input):

    with open(f"{subfolder}/{args.keyword}.log", "w") as f:
        IDs, seqs_incl_del = load_fasta(input)
        seqs = clean_seqs(seqs_incl_del)

        df = pd.DataFrame({'SequenceName': IDs, 'sequence': seqs, 'seq_incl_del': seqs_incl_del})
        query_ = df.iloc[:1]
        df = df.iloc[1:]
        L = len(df.sequence.iloc[0])
        df['frac_gaps'] = [x.count('-') / L for x in df['sequence']]

        df = df.loc[df.frac_gaps < float(args.gap_cutoff)]
        f.write(f"Filtered sequences by gap_cutoff={args.gap_cutoff}\n")
        
        df, clusters = get_labels(args, df)

        f.write(f"Found {len(clusters)} clusters using {args.cluster_method}\n")

        cluster_dir = os.path.join(subfolder, "clusters",)
        os.makedirs(cluster_dir, exist_ok=True)
        for clust in clusters:
            tmp = df.loc[df.dbscan_label == clust]
            out = pd.concat([query_, tmp], axis=0)
            outpath = os.path.join(cluster_dir, f"{args.keyword}_{clust:03d}.a3m")
            write_fasta(out.SequenceName.tolist(), out.seq_incl_del.tolist(), outfile=outpath)
            f.write(f"Wrote {outpath} (n={len(out)})\n")
    os.remove(f"{subfolder}/{args.keyword}.log")

def main(args):
    if args.input is not None:
        ids, seqs = load_fasta(args.input);
    if args.seq is not None:
        seqs = [args.seq]
        ids = [args.jobid]
    if args.msa is not None and len(ids) > 0:
        print(f'Assuming ID associated with input MSA is {id[0]}.')
        subfolder = os.path.join(args.outdir, ids[0])
        
        print(f'Running clustering...')
        run_cluster(args, subfolder, args.msa)

    for id_, seq_ in zip(ids, seqs): 
        args.keyword = id_
        subfolder = os.path.join(args.outdir, id_)
        os.makedirs(subfolder, exist_ok=True)
        msa_file = os.path.join(subfolder, f'{id_}.a3m')
        
        print(f'Running generating MSA...')
        if not os.path.exists(msa_file): 
            msa_seqs = run_mmseqs(seq_, args.tmpdir)
            with open(msa_file, "w") as a3m:
                a3m.write(msa_seqs[0])
        
        print(f'Running clustering...')
        run_cluster(args, subfolder, msa_file)

        print(f'Running structure prediction...')
        pred_dir = os.path.join(subfolder, 'preds')
        os.makedirs(pred_dir, exist_ok=True)
        seeds = generate_command(args, return_seed=True)
        for i in np.arange(0, seeds):
            for fil in glob.glob(f"{subfolder}/clusters/*.a3m"):
                fil_name = fil.split('/')[-1].strip('.a3m')
                os.makedirs(f'{pred_dir}/{fil_name}/s{i}', exist_ok=True)

                if os.path.exists(f'{pred_dir}/{fil_name}/s{i}/{fil_name}_0.done.txt'):
                    continue

                run_command = generate_command(args)
                run_command.extend(['--random-seed', f'{i}'])
                run_command.extend(['--jobname-prefix', f'{fil_name}'])
                run_command.extend([f'{fil}', f'{pred_dir}/{fil_name}/s{i}'])
                subprocess.run(run_command, shell=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument("--input", type=str, default=None, help="Input fasta")
    p.add_argument('--msa', type=str, default=None, help='Initial input MSA')
    p.add_argument('--seq', type=str, default=None, help='Initial input sequence')
    p.add_argument('--jobid', type=str, default='default')
    args = p.parse_args()

    if args.input is None and args.msa is None and args.seq is None:
        exit('Must specify either input or MSA or sequence')
    
    with open('configs/afcluster.yml', "r") as f:
        cfg = yaml.safe_load(f)
    
    for k, v in cfg.items():
        setattr(args, k, v)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.tmpdir, exist_ok=True)

    np.random.default_rng(seed=args.random_seed)
    
    os.environ['PATH'] = args.path_vars['PATH']
    os.environ["XDG_CACHE_HOME"] = args.path_vars['XDG_CACHE_HOME']
    os.environ["MPLCONFIGDIR"] = args.path_vars['MPLCONFIGDIR']
    
    main(args)
