import matplotlib.pyplot as plt
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
from scipy.stats import mannwhitneyu

sys.path.append(os.path.abspath('../src/'))

from utils import get_predictions_tail

def split_sort_return_str(line):
    # split a line on tabs, sort the elements and return a string
    return '\t'.join(sorted(line.split('\t')))

TEST_PATH: str =  '../ELs_for_Rotate/String_HPO_2019.all_hpo/test.txt'
TRAIN_PATH: str = '../ELs_for_Rotate/String_HPO_2019.all_hpo/train.txt'
VALID_PATH: str = '../ELs_for_Rotate/String_HPO_2019.all_hpo/valid.txt'

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

# get a set of genes
genes = set()
for line in open(TRAIN_PATH,'r'):
    row = line.strip().split('\t')
    if 'STRING:' in row[0]:
        genes.add(row[0])
    if 'STRING:' in row[2]:
        genes.add(row[2])
genes = list(genes)
# sort the genes
genes.sort()
print('# genes', len(genes))

# Load the STRING_HPO network edges
df: pd.DataFrame = pd.read_csv(
    "../ELs_for_Rotate/String_HPO_2019.all_hpo/train.txt", sep="\t"
)
df.columns = ["source", "relation", "target"]

# Load these edges into a NX graph and compute the degree for each entity
G: nx.MultiGraph = nx.from_pandas_edgelist(
    df, "source", "target", create_using=nx.MultiGraph()
)
degs: dict = dict(G.degree())

# Load the pretrained model
model = torch.load(
    "../PyKeenOut/stringhpo_rotate_trail_1.2/trained_model.pkl",
    map_location=torch.device("cpu"),
)

def get_percentile_dict(hpos,relation,data,model,degs,train_triples,test_triples,valid_triples):
    hpo2genes_percentiles = {}
    for i,term in enumerate(hpos):
        if i % 20 == 0:
            print(round(i/len(hpos),2) * 100, '% done')
        try:
            predictions_df: pd.DataFrame = get_predictions_tail(
                term,
                'STRING2HPO',
                data,
                model,
                degs,
                train_triples=train_triples,
                test_triples=test_triples,
                valid_triples=valid_triples
            )
        except KeyError:
            print('key error', term)
            continue
        # remove 'HP:' terms from the predictions
        predictions_df = predictions_df[~predictions_df['tail_label'].str.contains('HP:')]
        predictions_df['percentile'] = predictions_df['score'].rank(pct=True)
        # sort the predictions by the tail_label
        predictions_df = predictions_df.sort_values(by=['tail_label'])
        hpo2genes_percentiles[term] = predictions_df['percentile'].tolist()
    return hpo2genes_percentiles

def get_group_speficic_percentiles(all_genes,x_percentiles,new_edges_tp1):
    specific_phenotypes2genes = {}
    for term in x_percentiles.keys():
        specific_phenotypes2genes[term] = []
        tmp_seen_genes = set()
        for edge_line in new_edges_tp1:
            # if the term is in the edge and the edge is g2p
            if term in edge_line and 'STRING:' in edge_line and 'HP:' in edge_line:
                edge = edge_line.strip().split('\t')
                # get the gene
                gene = [x for x in edge if 'STRING:' in x][0]
                if gene in tmp_seen_genes:
                    continue
                tmp_seen_genes.add(gene)
                if gene not in all_genes:
                    print('gene not in genes', gene)
                    continue
                gi = all_genes.index(gene)
                specific_phenotypes2genes[term].append(x_percentiles[term][gi])
    return specific_phenotypes2genes

def two_group_g2p_ranking_test(prefix,hpos_a,hpos_b,label_a,label_b,data,model,degs,train_triples,test_triples,valid_triples,all_genes):
    # get the t edge list, this was the validation set for 2019
    t_el = '../ELs_for_Rotate/String_HPO_2019.all_hpo/test.txt'
    # load t_el as a set of string
    t_el_set = set([split_sort_return_str(line) for line in open(t_el,'r')])
    # ge the t+1 edge list, this was the validation set for 2020
    tp1_el = '../ELs_for_Rotate/String_HPO_2020.all_hpo/test.txt'
    # load tp1_el as a set of string
    tp1_el_set = set([split_sort_return_str(line) for line in open(tp1_el,'r')])
    new_edges_tp1 = tp1_el_set.difference(t_el_set)

    a_percentiles = get_percentile_dict(hpos_a,'STRING2HPO',data,model,degs,train_triples,test_triples,valid_triples)
    b_percentiles = get_percentile_dict(hpos_b,'STRING2HPO',data,model,degs,train_triples,test_triples,valid_triples)
    a_specific_percentiles = get_group_speficic_percentiles(all_genes,a_percentiles,new_edges_tp1)
    b_specific_percentiles = get_group_speficic_percentiles(all_genes,b_percentiles,new_edges_tp1)

    flattened_percentiles = {label_a:[],label_b:[]}

    for key in a_specific_percentiles.keys():
        flattened_percentiles[label_a] += list(a_specific_percentiles[key])

    for key in b_specific_percentiles.keys():
        flattened_percentiles[label_b] += list(b_specific_percentiles[key])

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].hist(flattened_percentiles[label_a],bins=10)
    axes[0].set_title(label_a)
    axes[0].set_xlabel('Percentile')
    axes[0].set_ylabel('Frequency')
    axes[1].hist(flattened_percentiles[label_b],bins=10)
    axes[1].set_title(label_b)
    axes[1].set_xlabel('Percentile')
    axes[1].set_ylabel('Frequency')

    # remove top and right borders
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    # match the y axis range
    max_y = max([max(axes[0].get_ylim()),max(axes[1].get_ylim())])
    axes[0].set_ylim([0,max_y])
    axes[1].set_ylim([0,max_y])


    plt.tight_layout()
    plt.savefig('../Figures/{}_{}_v_{}_g2p_rankings_hist.png'.format(prefix,label_a,label_b))
    # plt.show()

    # Mann Whitney is a non-parametric test intended for 2 groups
    print(len(flattened_percentiles[label_a]))
    print(len(flattened_percentiles[label_b]))
    flattened_percentiles[label_a].sort()
    flattened_percentiles[label_b].sort()
    print(flattened_percentiles[label_a])
    print(flattened_percentiles[label_b])
    U1, p = mannwhitneyu(flattened_percentiles[label_a], flattened_percentiles[label_b], method="exact")
    print('p-value', p)

    # what percent of the genes are in the top 10%?
    print('Percent of genes in top 10%')
    print(label_a, len([x for x in flattened_percentiles[label_a] if x > 0.8])/len(flattened_percentiles[label_a]))
    print(label_b, len([x for x in flattened_percentiles[label_b] if x > 0.8])/len(flattened_percentiles[label_b]))

    U1_8, p_8 = mannwhitneyu([x for x in flattened_percentiles[label_a] if x >= .8], [x for x in flattened_percentiles[label_b] if x >= .8], method="exact")
    print('p-value for scores > 0.8', p_8)

    return p
    

female_annotated_hpos: list = [ line.strip() for line in open('../Resources/just_female_hpos.txt')]
print('# female terms', len(female_annotated_hpos))
male_annotated_hpos: list = [ line.strip() for line in open('../Resources/just_male_hpos.txt')]
print('# female terms', len(male_annotated_hpos))

annotated_p = two_group_g2p_ranking_test(prefix='annotated',
                            hpos_a=female_annotated_hpos,
                            hpos_b=male_annotated_hpos,
                            label_a='female',
                            label_b='male',
                            data=data,
                            model=model,
                            degs=degs,
                            train_triples=train_triples,
                            test_triples=test_triples,
                            valid_triples=valid_triples,
                            all_genes=genes)