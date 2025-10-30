import numpy as np
from src.utils.seqs import *
from src.utils.helpers import *
from sklearn.cluster import DBSCAN

def cluster_DBSCAN(args, df):
    dbscan = dict_to_namespace(args.dbscan)
    L = len(df.sequence.iloc[0])
    n_clusters=[]
    ohe_seqs = encode_seqs(df.sequence.tolist(), max_len=L)
    if dbscan.eps_val is None:
        eps_test_vals = np.arange(dbscan.min_eps, dbscan.max_eps + dbscan.eps_step, dbscan.eps_step)
        for eps in eps_test_vals:
            testset = encode_seqs(df.sample(frac=args.gap_cutoff).sequence.tolist(), max_len=L)
            clustering = DBSCAN(eps=eps, min_samples=dbscan.min_samples).fit(testset)
            n_clust = len(set(clustering.labels_))
            n_not_clustered = len(clustering.labels_[np.where(clustering.labels_==-1)])
            n_clusters.append(n_clust)

        eps_to_select = eps_test_vals[np.argmax(n_clusters)]
        if np.argmax(n_clusters) == 0:
            eps_to_try = dbscan.min_eps
            scanner, best_n, best_eps = 0, 0, 0
            while scanner < 50:
                eps_to_try = eps_to_try + dbscan.eps_step/2
                clustering = DBSCAN(eps=eps_to_try, min_samples=dbscan.min_samples).fit(ohe_seqs)
                n_clust = len(set(clustering.labels_))
                n_not_clustered = len(clustering.labels_[np.where(clustering.labels_==-1)])
                if n_clust > best_n:
                    best_n = n_clust
                    scanner = 1
                    best_eps = eps_to_try
                if n_clust == best_n:
                    best_eps = eps_to_try
                    scanner = scanner + 1
                if n_clust > 0:
                    scanner = scanner + 1
            eps_to_select = best_eps

    else:
        eps_to_select = dbscan.eps_val

    clustering = DBSCAN(eps=eps_to_select, min_samples=dbscan.min_samples).fit(ohe_seqs)

    df['dbscan_label'] = clustering.labels_
    clusters = [x for x in df.dbscan_label.unique() if x>=0]

    return df, clusters