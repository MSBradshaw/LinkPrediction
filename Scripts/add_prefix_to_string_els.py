import argparse

# function to parse args
# params: el1, el2, nodemap, output_el1, output el2, output nodemap
def parse_args():
    parser = argparse.ArgumentParser(description="Create G2P lists")
    parser.add_argument('--edgelist1', help='tab seporated list of edges for full HPO')
    parser.add_argument('--edgelist2', help='tab seporated list of edges for just phenotypic HPO')
    # param for a node name to number mapping for time t
    parser.add_argument('--node_map1', help='tab seporated list of node name to number mapping for time t')
    parser.add_argument('--node_map2', help='tab seporated list of node name to number mapping for time t')
    parser.add_argument('--output_el1', help='output for el1 with string prefixes added')
    parser.add_argument('--output_el2', help='output for el1 with string prefixes added')
    parser.add_argument('--output_mapping1', help='output location for tab seporated list of node number to name mappings')
    parser.add_argument('--output_mapping2', help='output location for tab seporated list of node number to name mappings')
    return parser.parse_args()

# function that reads an a tab seporated el and writes it to a file with string prefixes added for every node that did not start with 'HP:'
def add_prefix_and_output(filein_name,fileout_name):
    with open(fileout_name,'w') as fileout:
        for line in open(filein_name,'r'):
            row = line.strip().split('\t')
            if 'HP:' not in row[0]:
                row[0] = 'STRING:' + row[0]
            if 'HP:' not in row[1]:
                row[1] = 'STRING:' + row[1]
            fileout.write('\t'.join(row) + '\n')

# function to add prefix in a node mapping file, only the second column is changed
def add_prefix_and_output_mapping(filein_name,fileout_name):
    with open(fileout_name,'w') as fileout:
        for line in open(filein_name,'r'):
            row = line.strip().split('\t')
            if 'HP:' not in row[1]:
                row[1] = 'STRING:' + row[1]
            fileout.write('\t'.join(row) + '\n')

def main():
    args = parse_args()
    add_prefix_and_output(args.edgelist1,args.output_el1)
    add_prefix_and_output(args.edgelist2,args.output_el2)
    add_prefix_and_output_mapping(args.node_map1,args.output_mapping1)
    add_prefix_and_output_mapping(args.node_map2,args.output_mapping2)

    
# if main then run main
if __name__ == '__main__':
    main()
