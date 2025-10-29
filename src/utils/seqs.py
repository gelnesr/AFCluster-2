import numpy as np

def encode_seqs(seqs, max_len=108, alphabet=None):
    if alphabet is None:
        alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    A = len(alphabet)
    lookup = {res: i for i, res in enumerate(alphabet)}

    arr = np.zeros((len(seqs), max_len, A), dtype=np.float32)
    for j, seq in enumerate(seqs):
        idxs = [lookup.get(c, A-1) for c in seq[:max_len]]  # unknowns -> '-'
        arr[j, np.arange(len(idxs)), idxs] = 1
    return arr.reshape(len(seqs), max_len * A)