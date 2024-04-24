import pandas as pd
import matplotlib.pyplot as plt
import typing
import networkx as nx
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Generate data for plotting kurve')
    parser.add_argument('--train_path', type=str, help='Path to training set')
    parser.add_argument('--valid_path', type=str, help='Path to validation set')
    parser.add_argument('--save_dir', type=str, help='Path to save the data')
    parser.add_argument('--df_path', type=str, help='Path to dataframe with scores')
    parser.add_argument('--filter', action='store_true', help='Filter on query term test set presence',default=False)
    return parser.parse_args()


def generate_plotting_data(df:pd.DataFrame, G:nx.Graph, save_dir:str=None, filter_on_query_test_prescense=True) -> pd.DataFrame:
    extra_suffix = 'normal'
    if filter_on_query_test_prescense:
        extra_suffix = 'filtered'
    # remove rows that were in the training set
    df['in_train'] = df.apply(lambda x: (x['head_label'], x['tail_label']) in G.edges, axis=1)
    df = df[~df['in_train']]
    # remove self edges
    df = df[df['head_label'] != df['tail_label']]
    # remove duplicate rows
    df = df.drop_duplicates(subset=['head_label', 'tail_label'])
    # sort by score in descending order
    df = df.sort_values(by='score', ascending=False)
    df['rank'] = df['score'].rank(pct=True)

    # find query terms that have atleast 1 test sest edge
    queries_to_keep = []
    for q in df['query_term'].unique():
        # where query term is q and in_test_set is True
        test_set_edges = df[(df['query_term'] == q) & (df['in_test_set'] == True)]
        if len(test_set_edges) > 0:
            queries_to_keep.append(q)
    print(queries_to_keep)
    queries_to_process = None

    if filter_on_query_test_prescense:
        print('Using filter')
        if len(queries_to_keep) == 0:
            # return empty dataframes
            # throw an error and message
            # print('No query terms have test set edges')
            # make an empty df but with the same columns
            avg_k_df = pd.DataFrame(columns=['k','hits@k','percent_hits@k'])
            scored_df = pd.DataFrame(columns=['head_id','score','head_label','tail_label','relation_label','head_degree','tail_degree','in_test_set','query_term','population','rank','in_train'])
            k_df = pd.DataFrame(columns=['k','hits@k','percent_hits@k','query_term'])
            print('No query terms have test set edges, returning empty dataframes')
            avg_k_df.to_csv(save_dir + f'avg_hits_at_k.{extra_suffix}.tsv', index=False, sep='\t')
            print('Saved', save_dir + f'avg_hits_at_k.{extra_suffix}.tsv')
            scored_df.to_csv(save_dir + f'edges_scored_by_query_term.{extra_suffix}.tsv', index=False, sep='\t')
            print('Saved', save_dir + f'edges_scored_by_query_term.{extra_suffix}.tsv')
            k_df.to_csv(save_dir + f'all_hits_at_k.{extra_suffix}.tsv', index=False, sep='\t')
            print('Saved', save_dir + f'all_hits_at_k.{extra_suffix}.tsv')
            return avg_k_df, scored_df, k_df
        queries_to_process = queries_to_keep
    else:
        print('Not using filtered')
        queries_to_process = df['query_term'].unique()

    k_df = None # will contain all the hits at K data
    scored_df = None # will contain the data with ranks assigned by query term instead of globally
    for q in queries_to_process:
        sub = df[df['query_term'] == q]
        # sort by score in descending order
        df = df.sort_values(by='score', ascending=False)
        df['rank'] = df['score'].rank(pct=True)
            # create a dataframe that is K, hits@K, percent hits@K
        k_values = [1,5,10,15,25] + list(range(50, 501, 50))
        hits_at_k = []
        percent_hits_at_k = []
        for k in k_values:
            hits_at_k.append(sub.head(k)['in_test_set'].sum())
            percent_hits_at_k.append(sub.head(k)['in_test_set'].sum() / k)
        sub_k_df = pd.DataFrame({'k': k_values, 'hits@k': hits_at_k, 'percent_hits@k': percent_hits_at_k})
        sub_k_df['query_term'] = q
        if k_df is None:
            k_df = sub_k_df
        else:
            k_df = pd.concat([k_df, sub_k_df])
        if scored_df is None:
            scored_df = sub
        else:
            scored_df = pd.concat([scored_df, sub])
    
    print(k_df.head())
    print(k_df.shape)
    # create an average hits@k for each k
    # get k_df without the query term column
    k_df_no_q = k_df.drop(columns=['query_term'])
    avg_k_df = k_df_no_q.groupby('k').median().reset_index()

    if save_dir:
        if save_dir[-1] != '/':
            save_dir += '/'
        os.makedirs(save_dir, exist_ok=True)
        avg_k_df.to_csv(save_dir + f'avg_hits_at_k.{extra_suffix}.tsv', index=False, sep='\t')
        print('Saved', save_dir + f'avg_hits_at_k.{extra_suffix}.tsv')
        scored_df.to_csv(save_dir + f'edges_scored_by_query_term.{extra_suffix}.tsv', index=False, sep='\t')
        print('Saved', save_dir + f'edges_scored_by_query_term.{extra_suffix}.tsv')
        k_df.to_csv(save_dir + f'all_hits_at_k.{extra_suffix}.tsv', index=False, sep='\t')
        print('Saved', save_dir + f'all_hits_at_k.{extra_suffix}.tsv')

    return avg_k_df, scored_df, k_df

def get_el_as_G(train_path, valid_path):
    el_df = pd.read_csv(train_path, sep='\t', header=None)
    el_df_valid = pd.read_csv(valid_path, sep='\t', header=None)
    el_df.columns = ['subject', 'predicate', 'object']
    el_df_valid.columns = ['subject', 'predicate', 'object']
    el_df = pd.concat([el_df, el_df_valid])

    G = nx.from_pandas_edgelist(el_df, 'subject', 'object', 'predicate')
    return G

def main():
    args = get_args()
    df = pd.read_csv(args.df_path, sep='\t')
    G = get_el_as_G(args.train_path, args.valid_path)
    avg_df, scored_df, k_df = generate_plotting_data(df, G, args.save_dir, filter_on_query_test_prescense=args.filter)

if __name__ == '__main__':
    main()

"""
python Scripts/generate_kurve_plotting_data.py --train_path ELs_for_Rotate/Monarch_KG_Filtered/train.txt \
                                            --valid_path ELs_for_Rotate/Monarch_KG_Filtered/valid.txt \
                                            --save_dir RankedResultsIntermediate/monarch_filtered/TransE/Cancer/ \
                                            --df_path GroupRankResults/monarch_filtered/TransE/AllResults/Cancer.comparison_results.tsv \
                                            --filter
"""