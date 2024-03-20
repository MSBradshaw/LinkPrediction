import argparse

# function to parse args

def parse_args():
    # input a tsv file with the output of the previously urn VCFtools command
    # output a filtered version of the input where only variants with AF > 0 in a single population are kept
    parser = argparse.ArgumentParser(description='Filter variants with AF > 0 in a single population')
    parser.add_argument('-i', '--input', type=str, metavar='', required=True, help='Input file')
    parser.add_argument('-o', '--output', type=str, metavar='', required=True, help='Output file')
    parser.add_argument('-m', '--mode', type=str, metavar='', required=False, help='Mode default all "a", pathogenic "p"')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    columns = None
    pop_columns = None
    all_56s = set()
    count_56s = {}
    with open(output_file, 'w') as out:
        for line in open(input_file, 'r'):
            line = line.strip()
            if line.startswith('CHROM'):
                columns = line.split('\t')
                # get the indexs of the population columns, those that are not CHROM, POS, REF, ALT, vep
                pop_columns = [i for i in range(len(columns)) if columns[i] not in ['CHROM', 'POS', 'REF', 'ALT', 'vep']]
                # write the header to the output file
                outline_list = ['CHROM', 'POS', 'REF', 'ALT']
                for i in pop_columns:
                    outline_list.append(columns[i])
                outline_list.append('gene')
                outline = '\t'.join(outline_list)
                out.write(outline + '\n')
                continue
            row = line.split('\t')
            # count the number of populations with AF > 0
            count = 0
            for i in pop_columns:
                try:
                    if row[i] != '?' and float(row[i]) > 0:
                        count += 1
                except ValueError:
                    print(i)
                    print(row[i])
                    print(line)
                    exit()
            # if only one population has AF > 0, write the line to the output file
            if count == 1:
                outline_list = [row[columns.index(i)] for i in ['CHROM', 'POS', 'REF', 'ALT']]
                for i in pop_columns:
                    outline_list.append(row[i])
                # get the gene symbol from the vep column
                vep = row[columns.index('vep')]
                gene = vep.split('|')[3] # index of the gene symbol
                if args.mode == 'p':
                    # if the mode is pathogenic, only keep genes with a pathogenic variant
                    all_56s.add(vep.split('|')[56])
                    if vep.split('|')[56] not in count_56s:
                        count_56s[vep.split('|')[56]] = 0
                    count_56s[vep.split('|')[56]] += 1
                    if 'pathogenic' != vep.split('|')[56]: # 56 is the CLIN_SIG column index
                        # if in pathogenic mode and the variant is not pathogenic, skip
                        continue
                outline_list.append(gene)
                outline = '\t'.join(outline_list)
                out.write(outline + '\n')
    print('All CLIN_SIG observed values', all_56s)
    for key in count_56s:
        print(key, count_56s[key])
            

        
            

if __name__ == '__main__':
    main()