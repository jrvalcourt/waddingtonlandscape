import os
import pandas as pd
import numpy as np
import pickle
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import seaborn as sb
import math

samplid2name = {
        'PREP0032_RHuanxxxxxA_v1_P4_FOXB2_S3'  : "FOXB2-OE",
        'PREP0032_RHuanxxxxxA_v1_P4_sorted_S1' : "unperturbed_P4",
        'PREP0095_RHuan11781A_A02v1_9_RH_973_S1'   : "973-OE",
        'PREP0095_RHuan11781A_B02v1_10_RH_1018_S2' : "1018-OE",
        'PREP0095_RHuan11781A_C02v1_11_RH_1024_S3' : "1024-OE",
        'PREP0095_RHuan11781A_D02v1_12_RH_1028_S4' : "1028-OE",
        'PREP0095_RHuan11781A_E02v1_13_RH_1053_S5' : "1053-OE",
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
        'ENCFF567GQW' : "neural progenitor",
        'ENCFF663ARH' : "neural progenitor",
        'ENCFF672VVX' : "neural progenitor",
        'ENCFF684BKA' : "neural crest"
}

def add_delta(n):
    return n + 0.01

# grab a list of TFs by Ensembl gene id, e.g. ENSG00000204531
# this list is available from BioMart
tfs = []
with open('ensembl_go0003700.txt') as fin:
    for line in fin:
        tfs.append(line.strip())

if not os.path.exists('data.p'):
    data1 = pd.read_table('data/abundance_genes_oe1.tsv', sep='\t')
    data2 = pd.read_table('data/abundance_genes_oe2.tsv', sep='\t')                     
    data  = pd.merge(data1, data2, on="ens_gene")
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
    
    pickle.dump(data, open('data.p', 'wb'))
data = pickle.load(open('data.p', 'rb'))

# clean data
data = data.dropna().transpose()
samples_to_use = [x for x in samplid2name]
data = data.loc[data.index.isin(samples_to_use)]

# consider only TFs
common_tfs = set(tfs).intersection(set(data.columns))
data = data[common_tfs]

# rename with more readable gene names
data = data.rename(index={x: samplid2name[x] for x in samplid2name})

# normalize the data
data = data.applymap(add_delta)
data = data.applymap(math.log10)
scaled_data = StandardScaler().fit_transform(data)

# make the clustermap
sb.clustermap(data.transpose().corr())
plt.savefig('plots/corr_mat.png')
plt.savefig('plots/corr_mat.pdf')
