import argparse
import networkx as nx

# arguments edgelist1, edgelist2, output
# function to parse args
def parse_args():
    parser = argparse.ArgumentParser(description="Create G2P lists")
    parser.add_argument('--edgelist1', help='tab seporated list of edges for time t')
    parser.add_argument('--edgelist2', help='tab seporated list of edges for time t+1')
    # param for a node name to number mapping for time t
    parser.add_argument('--node_map1', help='tab seporated list of node name to number mapping for time t')
    parser.add_argument('--output', help='output location for tab seporated list of g2p edges that are in el2 but not el1')
    parser.add_argument('--output_numbered', help='output location for tab seporated list of g2p edges that are in el2 but not el1')
    return parser.parse_args()

# function to read in edgelists
def read_edgelist(edgelist):
    # load el as an undirected graph
    G = nx.read_edgelist(edgelist, delimiter='\t', create_using=nx.Graph())
    return G

# function to remove non-g2p edges from one G
def remove_non_g2p(G):
    # exactly one of the edges must contain the substring 'HP:'
    for edge in G.edges():
        if sum([ 'HP:' in n for n in edge]) != 1:
            # the edge is PPI or HPO only
            G.remove_edge(*edge)
    return G

# remove any nodes in G2 that are not in G1
def remove_new_nodes(G1, G2):
    G2.remove_nodes_from([n for n in G2.nodes() if n not in G1.nodes()])
    return G2

# function to write the edges in G to a file tsv
def write_edgelist(G, output):
    with open(output, 'w') as f:
        for edge in G.edges():
            f.write('\t'.join(edge) + '\n')

# function to load a mapping file as a dictionary mapping node names to numbers
def load_node_map(node_map_file):
    node_map = {}
    with open(node_map_file) as f:
        for line in f:
            line = line.strip().split('\t')
            node_map[line[1]] = line[0]
    return node_map

# function to write the edges in G to a file tsv with node names replaced by numbers
def write_edgelist_numbered(G, node_map, output):
    with open(output, 'w') as f:
        for edge in G.edges():
            f.write('\t'.join([node_map[n] for n in edge]) + '\n')

# main function
def main():
    args = parse_args()
    G1 = read_edgelist(args.edgelist1)
    G2 = read_edgelist(args.edgelist2)
    G2.remove_edges_from(G1.edges())
    G2 = remove_non_g2p(G2)
    G2 = remove_new_nodes(G1, G2)
    write_edgelist(G2, args.output)
    node_map = load_node_map(args.node_map1)
    write_edgelist_numbered(G2, node_map, args.output_numbered)

if __name__ == '__main__':
    main()