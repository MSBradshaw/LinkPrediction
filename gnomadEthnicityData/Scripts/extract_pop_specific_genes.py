import argparse

# get params input a tsv file input and outputs a list of genes
def get_params():
    parser = argparse.ArgumentParser(description='Extracts genes from a tsv file')
    parser.add_argument('-i', '--input', type=str, help='Input tsv file')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    args = parser.parse_args()
    return args

def main():
    args = get_params()
    genes = {}
    columns = []
    pop_columns = []
    with open(args.input, 'r') as f:
        with open(args.output, 'w') as o:
            for line in f:
                row = line.strip().split('\t')
                if row[0] == 'CHROM':
                    columns = row
                    pop_columns = [(i,n) for i,n in enumerate(columns) if 'AF_' in n]
                    continue
                pop = None
                # the last column is the gene
                gene = row[-1]
                for i,n in pop_columns:
                    if row[i] != '?' and float(row[i]) > 0.00:
                        pop = n
                        break
                if pop is None:
                    print('No pop found for {}'.format(line))
                else:
                    if gene not in genes:
                        genes[gene] = set()
                    genes[gene].add(pop)
            for gene in genes:
                if len(genes[gene]) == 1:
                    o.write(gene + '\t' + list(genes[gene])[0] + '\n')
                        
                

if __name__ == '__main__':
    main()

