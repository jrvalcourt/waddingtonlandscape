import os
import pandas
import numpy as np
import Bio.motifs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn
import pickle
import pandas as pd
from scipy.cluster.hierarchy import linkage
import sys

np.random.seed(42)

tpm_cutoff = 3
deseq2_qval_cutoff = 0.000001

abundance_file = 'data/abundance_genes.tsv'
deseq2_file = 'data/biecto_results.tsv'
gene_to_motif_file = 'data/Human_TF_MotifList_v_1.01.txt'
motif_prepath = 'data/PWMs_lambert-et-al/'
ensid2genename_file = 'data/ensembl2genename.txt'
ensid2motif_patch_file = 'data/manual_tf_patch.txt'

mi_pickle = 'motif_mi_matrix.p'
clustering_method = 'average'
metric = 'braycurtis'

def complement(n):
    if n == 'A':
        return 'T'
    if n == 'T':
        return 'A'
    if n == 'G':
        return 'C'
    if n == 'C':
        return 'G'

def motif_dist(m1, m2):
    return m1.dist_product(m2)

def motif_corr(m1, m2):
    return m1.pssm.dist_pearson(m2.pssm)[0]

def motif_mi(m1, m2):
    len1 = len(m1.pwm.consensus)
    len2 = len(m2.pwm.consensus)

    max_mi = 0
    for offset in np.linspace(-(len1 - 1), (len2 - 1), 
                              num=len1+len2-1, dtype='int64'):
        curr_mi = 0
        curr_mi_rc = 0
        for ii in range(len1):
            if offset < 0:
                idx1 = ii
                idx2 = ii - offset
                idx1rc = min(len1, len2 + offset) - ii - 1
            else:
                idx1 = ii + offset
                idx2 = ii
                idx1rc = min(len1, len2) - ii - 1
            if idx1 >= len1 or idx2 >= len2:
                break
            for letter in list(m1.pwm):
                px = 0.25
                py = 0.25
                pxy    = float(m1.pwm[letter][idx1]) * float(m2.pwm[letter][idx2])
                pxy_rc = float(m1.pwm[complement(letter)][idx1rc]) * float(m2.pwm[letter][idx2])
                curr_mi    += pxy    * np.log2(pxy    / (px * py))
                curr_mi_rc += pxy_rc * np.log2(pxy_rc / (px * py))
        max_mi = max(max(curr_mi, max_mi), curr_mi_rc)
    return max_mi

def motif_kld(m1, m2):
    len1 = len(m1.pwm.consensus)
    len2 = len(m2.pwm.consensus)

    min_kld = float('inf')
    for offset in np.linspace(-(len1 - 1), (len2 - 1), 
                              num=len1+len2-1, dtype='int64'):
        curr_kld = 0
        curr_kld_rc = 0
        for ii in range(len1):
            if offset < 0:
                idx1 = ii
                idx2 = ii - offset
                idx1rc = min(len1, len2 + offset) - ii - 1
            else:
                idx1 = ii + offset
                idx2 = ii
                idx1rc = min(len1, len2) - ii - 1
            if idx1 >= len1 or idx2 >= len2:
                break
            for letter in list(m1.pwm):
                px    = float(m1.pwm[letter][idx1])
                px_rc = float(m1.pwm[complement(letter)][idx1rc])
                py    = float(m2.pwm[letter][idx2])
                curr_kld    += 0.5 * (px *    np.log2(px    / py) + py * np.log2(py / px))
                curr_kld_rc += 0.5 * (px_rc * np.log2(px_rc / py) + py * np.log2(py / px_rc))
        min_kld = min(min(curr_kld, min_kld), curr_kld_rc)
    return min_kld

def motif_from_ppm(ppm):
    ppm = ppm + 0.01
    minval = np.amin(np.amin(ppm))
    counts = ppm / minval
    m = Bio.motifs.Motif(counts=counts)
    return m

ensid2genename = {}
with open(ensid2genename_file) as fin:
    fin.readline()
    for line in fin:
        e = line.rstrip().split('\t')
        ensid2genename[e[0]] = e[1]

if os.path.exists(mi_pickle):
    df_mi = pickle.load(open(mi_pickle, 'rb'))
else:
    abundance_data = pandas.read_csv(abundance_file, sep='\t')
    
    ensid2motifid = {}
    ensid2motiffile = {}
    with open(gene_to_motif_file) as fin:
        fin.readline()
        for line in fin:
            e = line.rstrip().split('\t')
            if len(e) >= 8:
                ensid2motifid[e[0]] = e[6]
                ensid2motiffile[e[0]] = f"{motif_prepath}{e[6]}.txt"
            elif len(e) >= 7:
                pass
            else:
                pass
    
    with open(ensid2motif_patch_file) as fin:
        fin.readline()
        for line in fin:
            e = line.rstrip().split('\t')
            ensid2motifid[e[0]] = e[2]
            ensid2motiffile[e[0]] = '../grn_model/' + e[3]

    sig_ensids = []
    ensid2motif = {}
    with open(deseq2_file) as fin:
        fin.readline()
        for line in fin:
            e = line.rstrip().split('\t')
            try:
                qval = float(e[6])
            except ValueError:
                continue
            if qval >= 0.05:
                continue
            if e[0] in ensid2motifid:
                motif_path = ensid2motiffile[e[0]]
                if os.path.exists(motif_path):
                    sig_ensids.append(e[0])
                    tmp = pandas.read_csv(motif_path, index_col=0, header=0, sep='\t')
                    motif_ppm = tmp[['A', 'G', 'C', 'T']]
                    m = motif_from_ppm(motif_ppm)
                    ensid2motif[e[0]] = m
                else:
                    print(f"{motif_path} does not exist so skipping")
    
    ensids = []
    p4_cols = ['EXP0178_P4_S1', 'EXP0179_P4_S3', 'EXP0221_P4_S1', 'EXP0223_P4_S3']
    p7_cols = ['EXP0178_P7_S2', 'EXP0179_P7_S4', 'EXP0221_P7_S2', 'EXP0223_P7_S4']
    for query in sig_ensids:
        if ((abundance_data.loc[abundance_data['ens_gene'] == query, p4_cols] > tpm_cutoff).all(axis=None)):
            ensids.append(query)
        elif ((abundance_data.loc[abundance_data['ens_gene'] == query, p7_cols] > tpm_cutoff).all(axis=None)):
            ensids.append(query)
    for ensid in ensids:
        print(ensid2genename[ensid], ensid2motifid[ensid])
    
    print(len(ensids))
    mi_matrix = {}
    for ii, ensid1 in enumerate(ensids):
        mi_matrix[ensid1] = {}
        print(f"\r{ii / len(ensids) * 100:.2f}%       ", end="")
        for jj, ensid2 in enumerate(ensids):
            motif1 = ensid2motif[ensid1]
            motif2 = ensid2motif[ensid2]
            mi_matrix[ensid1][ensid2] = motif_mi(motif1, motif2)
    print("\rDone.                ")

    df_mi = pandas.DataFrame.from_dict(mi_matrix)
    pickle.dump(df_mi, open(mi_pickle, 'wb'))

mat = df_mi
mat = mat.rename(mapper=ensid2genename, axis=0)
mat = mat.rename(mapper=ensid2genename, axis=1)

link = linkage(mat, method=clustering_method, metric=metric)
seaborn.clustermap(mat, xticklabels=True, yticklabels=True, cmap="viridis",
                   robust=True, figsize=(20,20), row_linkage=link, col_linkage=link, vmax=16, vmin=2)
plt.savefig('motif_cluster.png', dpi=500)
