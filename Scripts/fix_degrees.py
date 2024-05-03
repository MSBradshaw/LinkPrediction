import pandas as pd
import argparse
import networkx as nx

def get_args():
    parser = argparse.ArgumentParser(description='Fix degrees')
    parser.add_argument('-i', dest='input', help='Input file')
    parser.add_argument('-e', dest='edgelist', help='Input file')
    parser.add_argument('-o', dest='output', help='Output file')
    return parser.parse_args()

def main():
    args = get_args()
    df = pd.read_csv(args.input, sep='\t')
    print(df.head())
    # load the semantic triples edgelist
    triples = pd.read_csv(args.edgelist, sep='\t')
    triples.columns = ['head','relation','tail']
    G = nx.from_pandas_edgelist(triples, 'head', 'tail')
    # get the degrees
    degrees = dict(G.degree())
    #head_label tail_label head_degree tail_degree
    df['head_degree'] = df['head_label'].map(degrees)
    df['tail_degree'] = df['tail_label'].map(degrees)
    df.to_csv(args.output, sep='\t', index=False)

main()