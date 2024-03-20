import argparse

# function to parse arguments
# input tsv of genes and populations
# output prefix for output files

def get_params():
    parser = argparse.ArgumentParser(description='Divides genes by population')
    parser.add_argument('-i', '--input', type=str, help='Input tsv file')
    parser.add_argument('-o', '--output', type=str, help='Output prefix')
    args = parser.parse_args()
    return args

def main():
    args = get_params()
    pops = {}
    for line in open(args.input,'r'):
        row = line.strip().split('\t')
        g,p = row[0],row[1]
        if p not in pops:
            pops[p] = []
        pops[p].append(g)
    for p in pops:
        with open(args.output + '.' + p + '.txt','w') as o:
            for g in pops[p]:
                o.write(g + '\n')
    with open(args.output + '.all.txt','w') as o:
        for p in pops:
            o.write(p + '\t' + str(len(pops[p])) + '\n')

if __name__ == '__main__':
    main()