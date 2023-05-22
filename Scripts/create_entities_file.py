import argparse

"""
This scrip will update the node number -> name mapping file to include nodes present in validation and test but not in train
"""

# function to collect params
# params: train, validate, test, g2p, output

def get_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--train', '-t', type=str, help='training triple list file')
    parser.add_argument('--validate', '-v', type=str, help='validation triple list file')
    parser.add_argument('--test', '-e', type=str, help='test triple list file')
    parser.add_argument('--mapping', '-m', type=str, help='number number to node name mapping file')
    parser.add_argument('--output', '-o', type=str, help='out entities file')
    args = parser.parse_args()
    return args

# function that reads a triple list and returns a set of nodes
def get_nodes_from_triple_list(triple_list):
    nodes = set()
    for line in open(triple_list, 'r'):
        if line.startswith('#'):
            continue
        # if the line is empty
        if line == '\n':
            continue
        row = line.strip().split('\t')
        nodes.add(row[0])
        nodes.add(row[2])
    return nodes

def ss(line):
    return line.strip().split('\t')

def main():
    args = get_params()
    train_nodes = get_nodes_from_triple_list(args.train)
    validate_nodes = get_nodes_from_triple_list(args.validate)
    test_nodes = get_nodes_from_triple_list(args.test)
    # list of nodes in validation and test but not in train
    new_nodes = validate_nodes.union(test_nodes).difference(train_nodes)
    # read the mapping file in as a dictionary
    node_map = { ss(line)[1]:int(ss(line)[0]) for line in open(args.mapping, 'r')}
    # add the new nodes to the mapping, each id should increase by 1 start from where node_map left off
    for node in new_nodes:
        node_map[node] = len(node_map)
    # write the new mapping file
    with open(args.output, 'w') as out:
        for node in node_map:
            out.write('{}\t{}\n'.format(node_map[node], node))



if __name__ == '__main__':
    main()