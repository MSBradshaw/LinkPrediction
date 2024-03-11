import argparse
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import torch
from pykeen.datasets.base import PathDataset
from scipy.stats import kruskal
from typing import List
from pykeen import predict
import os
import pickle
import warnings

# repress the SettingWithCopyWarning warning
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

# add the current dir to the path
import sys
sys.path.append('Scripts/')
from compare_groups_test_omatic import get_scores_for_edges, read_terms_from_file, load_degs, get_triples, annotate

def get_args():
    parser = argparse.ArgumentParser(description='Run the test-omatic experiment')
    parser.add_argument('--terms', type=str, required=True, help='list of terms, terms must be in the training set')
    parser.add_argument('--label', type=str, required=True, help='label for the terms group')
    parser.add_argument('--relation', type=str, required=True, help='relation to use for predictions')
    parser.add_argument('--prediction_target', type=str, required=True, help='"head" or "tail", which entity are you trying to predict? Edges are assumed to be directed.')
    parser.add_argument('--prediction_prefix', type=str, required=True, help='the type of entity to predict, the prefix of the entity id eg "MONDO:" "HGNC:"')
    parser.add_argument('--train_triples', type=str, required=True, help='training triple list')
    parser.add_argument('--validation_triples', type=str, required=True, help='validation triple list')
    parser.add_argument('--test_triples', type=str, required=True, help='test triple list')
    parser.add_argument('--model', type=str, required=True, help='trained pykeen model in .pkl format')
    parser.add_argument('--output', type=str, required=True, help='output file name')
    parser.add_argument('--progress_bar', dest='progress_bar', action='store_true',default=False, help='show progress bar')
    args = parser.parse_args()
    return args

def check_terms_in_training(terms: List[str], train_triples: List[List[str]], outputfile:str) -> List[str]:
    """
    Ensure that all query terms are in the training set and of type HGNC or MONDO, we dont care about other types
    """
    terms_in_training = []
    training_nodes = set([triple[0] for triple in train_triples] + [triple[2] for triple in train_triples])
    missing_nodes = []
    for term in terms:
        if term in training_nodes and ('HGNC:' in term or 'MONDO:' in term):
            terms_in_training.append(term)
        else:
            missing_nodes.append(term)
    if len(terms_in_training) == 0:
        raise ValueError('None of the terms are in the training set')
    
    with open(outputfile + '.warning','w') as f:
        f.write('Nodes not in training edges\n')
        for node in missing_nodes:
            f.write(node + '\n')
    return terms_in_training

def check_content(row, triples):
    return (row['head_label'],row['relation_label'],row['tail_label']) in triples or (row['tail_label'],row['relation_label'],row['head_label']) in triples

def label_scores_df(df: pd.DataFrame, term: str, test_triples: List[tuple]) -> pd.DataFrame:
    # get the test triples that have a term in them
    test_triples_with_term = [triple for triple in test_triples if term in triple]
    # label the test triples in the df
    df['in_test_set'] = [check_content(row, test_triples_with_term) for i,row in df.iterrows()]
    return df

def get_scores_for_edges(terms: List[str],
                         relation: str,
                         model,
                         degs: dict,
                         train_triples: List[List[str]],
                         test_triples: List[List[str]],
                         validation_triples: List[List[str]],
                         population: str,
                         prefix_of_predicted: str,
                         prediction_target: str,
                         data: PathDataset,
                         progress_bar: bool = False) -> (List[float], pd.DataFrame):
    all_results_df = None
    for i,term in enumerate(terms):
        print(term)
        try:
            if i % (len(terms) // 50) == 0 and progress_bar:
                print(f'{i} / {len(terms)}')
        except ZeroDivisionError:
            pass

        if prediction_target == 'head':
            pred_df = predict.predict_target(model=model, head=None, relation=relation, tail=term, triples_factory=data).df            
            og_col = 'tail_label'
            new_col = 'head_label'
            target_edge_index = 0
        elif prediction_target == 'tail':
            pred_df = predict.predict_target(model=model, head=term, relation=relation, tail=None, triples_factory=data).df
            og_col = 'head_label'
            new_col = 'tail_label'
            target_edge_index = 2
        else:
            raise ValueError('prediction_target must be "head" or "tail"')
        pred_df = annotate(df=pred_df,data=data,og_label=term,og_col=og_col,new_col=new_col,relation=relation, prefix_of_predicted=prefix_of_predicted,degs=degs)
        pred_df = label_scores_df(pred_df, term, test_triples)
        pred_df['query_term'] = term
        pred_df['population'] = population
        # sort the predictions by score
        pred_df = pred_df.sort_values(by=['score'],ascending=False)
        # assign a percentile to each prediction
        pred_df['rank'] = pred_df['score'].rank(pct=True)
        # remove the row where head and tail are the same
        pred_df = pred_df[pred_df['head_label'] != pred_df['tail_label']]

        if all_results_df is None:
            all_results_df = pred_df
        else:
            all_results_df = pd.concat([all_results_df,pred_df])
    return all_results_df

def main():
    args = get_args()
    model = torch.load(
        args.model,
        map_location=torch.device("cpu"),
    )

    # load terms
    terms = read_terms_from_file(args.terms)
    
    # load network to get degrees
    degs = load_degs(args.train_triples)

    # load triples
    train_triples, test_triples, validation_triples = get_triples(args.train_triples,args.validation_triples,args.test_triples)

    terms = check_terms_in_training(terms, train_triples, args.output)

    # write test triples to del.tsv
    with open('del.tsv','w') as outfile:
        for edge in test_triples:
            outfile.write('\t'.join(edge) + '\n')

    # load data as factory
    class TheKG(PathDataset):
        def __init__(self, **kwargs):
            super().__init__(
                training_path=args.train_triples,
                testing_path=args.test_triples,
                validation_path=args.validation_triples,
                **kwargs,
            )
    data = TheKG()


    all_results_df = get_scores_for_edges(terms = terms,
                            relation=args.relation,
                            model=model,
                            degs=degs,
                            train_triples=train_triples,
                            test_triples=test_triples,
                            validation_triples=validation_triples,
                            population=args.label,
                            prefix_of_predicted=args.prediction_prefix,
                            prediction_target=args.prediction_target,
                            data=data,
                            progress_bar=args.progress_bar)

    all_results_df.to_csv(args.output,sep='\t',index=False)

# if main run it
if __name__ == '__main__':
    main()