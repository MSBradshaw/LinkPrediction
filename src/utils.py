# -*- coding: utf-8 -*-

"""
This module originally came from https://github.com/AstraZeneca/biomedical-kg-topological-imbalance/
but has been modified to work with pykeen version 1.10.1 and my specific use cases.
"""

import pandas as pd
from pykeen.datasets.base import Dataset
from pykeen.models import Model
from pykeen import predict
import matplotlib.pyplot as plt

def create_type(row) -> str:    
    if row["in_training"] is True and row["in_testing"] is False:
        return "train"
    elif row["in_training"] is False and row["in_testing"] is True:
        return "test"
    elif row["in_training"] is False and row["in_testing"] is False:
        return "novel"
    else:
        return "unknown"

def annotate_predicted_df(
    df: pd.DataFrame,
    degs: dict,
    position: str,
) -> pd.DataFrame:
    """Annotate a pykeen predictions dataframe."""

    df["entity_type"] = df[position].str.split(":", expand=True)[0]
    df["triple_type"] = df.apply(lambda row: create_type(row), axis=1)
    df["deg"] = [degs[e] for e in list(df[position].values)]

    return df

def make_string_triple(q_entity, q_relation, tail, data) -> str:
     return str([q_entity, q_relation, data.training.entity_id_to_label[tail]])

def get_predictions_tail(
    q_entity: str,
    q_relation: str,
    data: Dataset,
    model: Model,
    degs: dict,
    train_triples: set,
    test_triples: set,
    valid_triples: set
) -> pd.DataFrame:
    """Make a prediction using a a partial triple (missing tail)."""

    pred_df = predict.predict_target(
        model=model,
        head=q_entity,
        relation=q_relation,
        triples_factory=data
    ).df

    pred_df['in_training'] = [ make_string_triple(q_entity, q_relation, x, data) in train_triples or make_string_triple(q_entity, q_relation, x, data)[::-1] in train_triples for x in pred_df['tail_id']]
    pred_df['in_testing'] = [ make_string_triple(q_entity, q_relation, x, data) in test_triples or make_string_triple(q_entity, q_relation, x, data)[::-1] in test_triples for x in pred_df['tail_id']]

    pred_df['tail_label'] = [data.training.entity_id_to_label[x] for x in pred_df['tail_id']]

    pred_df = annotate_predicted_df(pred_df, degs, "tail_label")

    return pred_df


def get_predictions_head(
    q_entity: str,
    q_relation: str,
    data: Dataset,
    model: Model,
    degs: dict,
) -> pd.DataFrame:
    """Make a prediction using a a partial triple (missing head)."""

    pred_df = model.get_head_prediction_df(
        q_relation,
        q_entity,
        triples_factory=data.training,
        testing=data.testing.mapped_triples,
    )
    pred_df = annotate_predicted_df(pred_df, degs, "head_label")

    return pred_df

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

def get_group_specific_percentiles(all_genes,x_percentiles,new_edges_tp1):
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

def plot_degs_vs_percentiles(flattened_percentiles,label_a_degs,label_b_degs,prefix, label_a, label_b):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].scatter(label_a_degs,flattened_percentiles[label_a])
    axes[0].set_title(label_a)
    axes[0].set_xlabel('Degree')
    axes[0].set_ylabel('Percentile')
    axes[1].scatter(label_b_degs,flattened_percentiles[label_b])
    axes[1].set_title(label_b)
    axes[1].set_xlabel('Degree')
    axes[1].set_ylabel('Percentile')

    # log scale the x axis
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    

    # remove top and right borders
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    # match the y axis range
    max_y = max([max(axes[0].get_ylim()),max(axes[1].get_ylim())])
    axes[0].set_ylim([0,max_y])
    axes[1].set_ylim([0,max_y])
    # match x axis range
    max_x = max([max(axes[0].get_xlim()),max(axes[1].get_xlim())])
    axes[0].set_xlim([0,max_x])
    axes[1].set_xlim([0,max_x])

    plt.tight_layout()
    plt.savefig('../Figures/{}_{}_v_{}_g2p_rankings_deg.png'.format(prefix,label_a,label_b))

def two_group_g2p_ranking_test(prefix,hpos_a,hpos_b,label_a,label_b,data,model,degs,train_triples,test_triples,valid_triples,all_genes,relation):
    # get the t edge list, this was the validation set for 2019
    t_el = '../ELs_for_Rotate/String_HPO_2019.all_hpo/test.txt'
    # load t_el as a set of string
    t_el_set = set([split_sort_return_str(line) for line in open(t_el,'r')])
    # ge the t+1 edge list, this was the validation set for 2020
    tp1_el = '../ELs_for_Rotate/String_HPO_2020.all_hpo/test.txt'
    # load tp1_el as a set of string
    tp1_el_set = set([split_sort_return_str(line) for line in open(tp1_el,'r')])
    new_edges_tp1 = tp1_el_set.difference(t_el_set)

    a_percentiles = get_percentile_dict(hpos_a,relation,data,model,degs,train_triples,test_triples,valid_triples)
    b_percentiles = get_percentile_dict(hpos_b,relation,data,model,degs,train_triples,test_triples,valid_triples)
    print('a_percentiles',len(a_percentiles))
    print('b_percentiles',len(b_percentiles))
    a_specific_percentiles = get_group_specific_percentiles(all_genes,a_percentiles,new_edges_tp1)
    b_specific_percentiles = get_group_specific_percentiles(all_genes,b_percentiles,new_edges_tp1)
    print('a_specific_percentiles',len(a_specific_percentiles))
    print('b_specific_percentiles',len(b_specific_percentiles))

    flattened_percentiles = {label_a:[],label_b:[]}

    label_a_degs = []
    label_b_degs = []

    for key in a_specific_percentiles.keys():
        flattened_percentiles[label_a] += list(a_specific_percentiles[key])
        if key in degs:
            d = [degs[key]] * len(a_specific_percentiles[key])
        else:
            d = [0] * len(a_specific_percentiles[key])
        label_a_degs += d


    for key in b_specific_percentiles.keys():
        flattened_percentiles[label_b] += list(b_specific_percentiles[key])
        if key in degs:
            d = [degs[key]] * len(b_specific_percentiles[key])
        else:
            d = [0] * len(b_specific_percentiles[key])
        label_b_degs += d
    
    if len(flattened_percentiles[label_a]) == 0 or len(flattened_percentiles[label_b]) == 0:
        print('one of the groups has new g2p edges added in year t+1')
        return 1.0

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    # make bins 0 -  1 in increments of .1
    bins = np.arange(0, 1.1, 0.1)
    axes[0].hist(flattened_percentiles[label_a],bins=bins)
    axes[0].set_title(label_a)
    axes[0].set_xlabel('Percentile')
    axes[0].set_ylabel('Frequency')
    axes[1].hist(flattened_percentiles[label_b],bins=bins)
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
    plt.show()

    plot_degs_vs_percentiles(flattened_percentiles,label_a_degs, label_b_degs,prefix,label_a,label_b)

    # Kruskal is a non-parametric test intended for multiple groups
    print('flattened_percentiles[label_a]',len(flattened_percentiles[label_a]))
    print('flattened_percentiles[label_b]',len(flattened_percentiles[label_b]))
    flattened_percentiles[label_a].sort()
    flattened_percentiles[label_b].sort()
    print(flattened_percentiles[label_a])
    print(flattened_percentiles[label_b])
    U1, p = kruskal(flattened_percentiles[label_a], flattened_percentiles[label_b])
    print('p-value', p)

    # what percent of the genes are in the top 10%?
    print('Percent of genes in top 10%')
    print(label_a, len([x for x in flattened_percentiles[label_a] if x > 0.8])/len(flattened_percentiles[label_a]))
    print(label_b, len([x for x in flattened_percentiles[label_b] if x > 0.8])/len(flattened_percentiles[label_b]))

    U1_8, p_8 = kruskal([x for x in flattened_percentiles[label_a] if x >= .8], [x for x in flattened_percentiles[label_b] if x >= .8])
    print('p-value for scores > 0.8', p_8)

    return p
    
def split_sort_return_str(line):
    # split a line on tabs, sort the elements and return a string
    return '\t'.join(sorted(line.split('\t')))