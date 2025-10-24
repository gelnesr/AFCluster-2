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
            testset = encode_seqs(df.sample(frac=0.25).sequence.tolist(), max_len=L)
            clustering = DBSCAN(eps=eps, min_samples=dbscan.min_samples).fit(testset)
            n_clust = len(set(clustering.labels_))
            n_not_clustered = len(clustering.labels_[np.where(clustering.labels_==-1)])
            n_clusters.append(n_clust)

        eps_to_select = eps_test_vals[np.argmax(n_clusters)]
    else:
        eps_to_select = dbscan.eps_val

    clustering = DBSCAN(eps=eps_to_select, min_samples=dbscan.min_samples).fit(ohe_seqs)

    df['dbscan_label'] = clustering.labels_
    clusters = [x for x in df.dbscan_label.unique() if x>=0]

    return df, clusters