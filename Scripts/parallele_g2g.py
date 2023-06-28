# import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from pykeen.datasets.hetionet import Hetionet
from scipy import stats
import sys
import os
from pykeen.datasets.base import PathDataset
from pykeen.predict import predict_target
from scipy.stats import kstest
import pickle

sys.path.append(os.path.abspath('src/'))

from utils import get_predictions_tail

# start and end incides for the genes to run so this can be parallelized
start = int(sys.argv[1])
end = int(sys.argv[2])


TEST_PATH: str =  'ELs_for_Rotate/String_HPO_2019.all_hpo/test.txt'
TRAIN_PATH: str = 'ELs_for_Rotate/String_HPO_2019.all_hpo/train.txt'
VALID_PATH: str = 'ELs_for_Rotate/String_HPO_2019.all_hpo/valid.txt'

class STRINGHPO(PathDataset):
    def __init__(self, **kwargs):
        super().__init__(
            training_path=TRAIN_PATH,
            testing_path=TEST_PATH,
            validation_path=VALID_PATH,
            **kwargs,
        )

data = STRINGHPO()
data.summarize() 

# these are used to know if a triple is in the training, testing or validation set
train_triples = set([ str(x) for x in data.training.triples])
test_triples = set([ str(x) for x in data.testing.triples])
valid_triples = set([ str(x) for x in data.validation.triples])

# Load the STRING_HPO network edges
df: pd.DataFrame = pd.read_csv(
    "ELs_for_Rotate/String_HPO_2019.all_hpo/train.txt", sep="\t"
)
df.columns = ["source", "relation", "target"]

# Load these edges into a NX graph and compute the degree for each entity
G: nx.MultiGraph = nx.from_pandas_edgelist(
    df, "source", "target", create_using=nx.MultiGraph()
)
degs: dict = dict(G.degree())

# Load the pretrained model
model = torch.load(
    "PyKeenOut/stringhpo_rotate_trail_1.2/trained_model.pkl",
    map_location=torch.device("cpu"),
)

# this cell originally took 425 minutes to run - results have been pickled in the following cell
# load sex specific data
female_hpos: list = [ line.strip() for line in open('Resources/just_female_hpos.txt')]
print('# female terms', len(female_hpos))
male_hpos: list = [ line.strip() for line in open('Resources/just_male_hpos.txt')]
print('# male terms', len(male_hpos))

# get a set of genes
genes = set()
for line in open(TRAIN_PATH,'r'):
    row = line.strip().split('\t')
    if 'STRING:' in row[0]:
        genes.add(row[0])
    if 'STRING:' in row[2]:
        genes.add(row[2])
print('# genes', len(genes))

def update_hpo_percentiles(hp_dict: dict,df: pd.DataFrame, hpos: list) -> dict: 
    for h in hpos:
        if h not in hp_dict:
            hp_dict[h] = []
        row = df.loc[df['tail_label'] == h,'percentile']
        if row.size == 0:
            continue
        p = row.values[0]
        hp_dict[h].append(p)
    return hp_dict


xls = pd.ExcelFile('Resources/aba3066-table-s2.xlsx')

gene_tissue_counts_female = {}
gene_tissue_counts_male = {}
t2g2s = {'tissue':[],'gene':[],'sex':[]}
tissue_2_gene_by_sex = {}
for sheet_name in xls.sheet_names[2:]:
    xdf = xls.parse(sheet_name)
    for i,row in xdf.iterrows():
        g = row['HUGO_gene_id']
        is_female = row['MASH beta'] > 0
        if is_female:
            if g not in gene_tissue_counts_female:
                gene_tissue_counts_female[g] = 0
            gene_tissue_counts_female[g] += 1
            t2g2s['tissue'].append(sheet_name)
            t2g2s['gene'].append(g)
            t2g2s['sex'].append('female')
        else:
            if g not in gene_tissue_counts_male:
                gene_tissue_counts_male[g] = 0
            gene_tissue_counts_male[g] += 1
            t2g2s['tissue'].append(sheet_name)
            t2g2s['gene'].append(g)
            t2g2s['sex'].append('male')

# build a DF for these that are differentially expressed and with greater expression in females
genes_female = gene_tissue_counts_female.keys()
data_female = {'gene':genes_female,'tissue_count':[gene_tissue_counts_female[g] for g in genes_female]}
df_female = pd.DataFrame(data_female)
df_female['tissue_count'].max()
print('females',df_female[df_female['tissue_count'] == df_female['tissue_count'].max()])
df_female = df_female[df_female['tissue_count'] < 45]
print('females',df_female.shape)

# build a DF for these that are differentially expressed and with greater expression in males
genes_male = gene_tissue_counts_male.keys()
data_male = {'gene':genes_male,'tissue_count':[gene_tissue_counts_male[g] for g in genes_male]}
df_male = pd.DataFrame(data_male)
df_male['tissue_count'].max()
print('males',df_male[df_male['tissue_count'] == df_male['tissue_count'].max()])
df_male = df_male[df_male['tissue_count'] < 45]
print('males',df_male.shape)

female_genes = ['STRING:' + x for x in df_female['gene']]
male_genes =  ['STRING:' + x for x in df_male['gene']]

print(female_genes)

query_relation: str = "STRING2STRING"

all_female_g2g_percentiles = {}
all_male_g2g_percentiles = {}

predictions_df: pd.DataFrame = get_predictions_tail(
    list(genes)[100],
    query_relation,
    data,
    model,
    degs,
    train_triples=train_triples,
    test_triples=test_triples,
    valid_triples=valid_triples
)

for i,gene in enumerate(genes):
    print('gene: ', gene, i)
    if i % 1000 == 0:
        print(round(i/len(genes),2) * 100, '% done')
    if i < start:
        continue
    if i > end:
        break
    try:
        predictions_df: pd.DataFrame = get_predictions_tail(
            gene,
            query_relation,
            data,
            model,
            degs,
            train_triples=train_triples,
            test_triples=test_triples,
            valid_triples=valid_triples
        )
    except KeyError:
        print('key error', gene)
        continue
    predictions_df['percentile'] = predictions_df['score'].rank(pct=True)
    all_female_g2g_percentiles = update_hpo_percentiles(all_female_g2g_percentiles, predictions_df, female_genes)
    all_male_g2g_percentiles = update_hpo_percentiles(all_male_g2g_percentiles, predictions_df, male_genes)
    print(len(all_female_g2g_percentiles))
    print(len(all_male_g2g_percentiles))
    print(len(all_female_g2g_percentiles[female_genes[0]]))
    print(len(all_male_g2g_percentiles[male_genes[0]]))
    print()

pickle.dump(all_female_g2g_percentiles, open('G2G_Shards/all_female_g2g_percentiles.part_{}.pkl'.format(str(start) + '_' + str(end)),'wb'))
pickle.dump(all_male_g2g_percentiles, open('G2G_Shards/all_male_g2g_percentiles.part_{}.pkl'.format(str(start) + '_' + str(end)),'wb'))