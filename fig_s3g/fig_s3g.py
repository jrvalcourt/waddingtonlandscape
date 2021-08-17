import os
import pandas as pd
import numpy as np
import pickle
import umap
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import seaborn as sb
import math

samplid2name = {
        'PREP0032_RHuanxxxxxA_v1_P4_sorted_S1' : "pre",
        'PREP0032_RHuanxxxxxA_v1_P7_sorted_S2' : "post",
        'ENCFF034KRQ' : "ectoderm",
        'ENCFF044YLS' : "mesendoderm",
        'ENCFF081JBX' : "neural",
        'ENCFF145AQN' : "neural progenitor",
        'ENCFF183XSM' : "neuronal stem",
        'ENCFF290ZZQ' : "neural crest",
        'ENCFF342LYI' : "trophoblast",
        'ENCFF419KMW' : "ectoderm",
        'ENCFF425FGL' : "excitatory neuron",
        'ENCFF466QUZ' : "mesendoderm",
        'ENCFF483MRL' : "excitatory neuron",
        'ENCFF538XVQ' : "ectoderm",
        'ENCFF684BKA' : "neural crest",
        'ENCFF691ZYQ' : "ectoderm",
        'ENCFF699LBP' : "ectoderm",
        'ENCFF760HDK' : "trophoblast",
        'ENCFF768SPT' : "ectoderm"
}

def add_delta(n):
    return n + 0.01

tfs = []
with open('ensembl_go0003700.txt') as fin:
    for line in fin:
        tfs.append(line.strip())

if not os.path.exists('data_pre_post.p'):
    data = pd.read_table('data/abundance_genes.tsv', sep='\t')
    ref_genes = set(data['ens_gene'])
    data = data.set_index('ens_gene')
    
    indir = 'encode_comparison_data'
    for f in os.listdir(indir):
        if not f.endswith('.tsv'):
            continue
        path = os.path.join(indir, f)
        enc_id = f.split('.')[0]
        temp = pd.read_table(path, sep="\t")
        temp = temp.rename(columns={"gene_id": "ens_gene", "TPM": enc_id})
        temp['ens_gene'] = [x.split('.')[0] for x in temp['ens_gene']]
        temp = temp.loc[[x.startswith('ENSG') for x in temp['ens_gene']]]
        temp = temp.loc[temp['ens_gene'].isin(ref_genes)]
        temp = temp.drop_duplicates(subset='ens_gene')
        temp = temp.set_index('ens_gene')
        data = data.merge(temp[[enc_id]], how='left', left_index=True, right_index=True)
    
    pickle.dump(data, open('data_pre_post.p', 'wb'))
data = pickle.load(open('data_pre_post.p', 'rb'))


# clean data
data = data.dropna().transpose()
samples_to_use = [x for x in samplid2name]

data = data.loc[data.index.isin(samples_to_use)]
common_tfs = set(tfs).intersection(set(data.columns))
data = data[common_tfs]
data = data.rename(index={x: samplid2name[x] for x in samplid2name})
data = data.applymap(add_delta)
scaled_data = StandardScaler().fit_transform(data)

# make the plots
sb.clustermap(data.transpose().corr())
plt.savefig('plots/corr_mat_pre_post.png')

sb.clustermap(data.transpose().corr())
plt.savefig('plots/corr_mat_pre_post.pdf')
