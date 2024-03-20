import argparse

# get params input a tsv file input and outputs a list of genes
def get_params():
    parser = argparse.ArgumentParser(description='Extracts genes from a tsv file')
    parser.add_argument('-i', '--input', type=str, help='Input tsv file')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('-m', '--mode', type=str, help='Mode to run in, strict "s", majority "m", plurality "m"', default='s')
    args = parser.parse_args()
    return args

def main():
    args = get_params()
    genes = {}
    columns = []
    pop_columns = []
    if args.mode == 's':
        """
        If all of the populations associated with ancestroy specific variants in the gene are the same, then the gene is considered to be ancestroy specific
        """
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
    elif args.mode == 'm':
        """
        If a majority (>50%) of the populations associated with ancestroy specific variants in the gene are the same, then the gene is considered to be ancestroy specific
        """
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
                            genes[gene] = []
                        genes[gene].append(pop)
                for gene in genes:
                    # find the most common pop
                    pop_counts = {}
                    for pop in genes[gene]:
                        if pop not in pop_counts:
                            pop_counts[pop] = 0
                        pop_counts[pop] += 1
                    max_pop = None
                    max_count = 0
                    for pop in pop_counts:
                        if pop_counts[pop] > max_count:
                            max_count = pop_counts[pop]
                            max_pop = pop
                    if max_pop is None:
                        print('No pop found for {}'.format(line))
                    elif max_count > len(genes[gene])/2:
                        o.write(gene + '\t' + list(genes[gene])[0] + '\n')
    elif args.mode == 'p':
        """
        If plurality of the populations associated with ancestroy specific variants in the gene are the same, then the gene is considered to be ancestroy specific
        """
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
                            genes[gene] = []
                        genes[gene].append(pop)
                for gene in genes:
                    # find the most common pop
                    pop_counts = {}
                    for pop in genes[gene]:
                        if pop not in pop_counts:
                            pop_counts[pop] = 0
                        pop_counts[pop] += 1
                    max_pop = None
                    max_count = 0
                    tie = False
                    for pop in pop_counts:
                        if pop_counts[pop] > max_count:
                            max_count = pop_counts[pop]
                            max_pop = pop
                            tie = False
                        elif pop_counts[pop] == max_count:
                            tie = True
                            
                    if max_pop is None:
                        print('No pop found for {}'.format(line))
                    elif max_count > len(genes[gene])/2 and not tie:
                        o.write(gene + '\t' + list(genes[gene])[0] + '\n')
                        
                

if __name__ == '__main__':
    main()

