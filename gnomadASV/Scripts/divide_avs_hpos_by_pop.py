import argparse

# function for argparse
"""
params
input: tsv
output prefix

"""

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help="tsv file with pathogenic ASVs")
    parser.add_argument('-o', '--output', type=str, required=True, help="output prefix")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    pops = {}
    for line in open(args.input,'r'):
        row = line.strip().split('\t')
        if row[5] not in pops: # 5 is the population column
            pops[row[5]] = []
        pops[row[5]].append(line)

    for pop in pops:
        with open(args.output + '_' + pop + '.tsv','w') as out:
            out.write('\n'.join(pops[pop]))
    
    # report number of variants per pop and the mean and median number of HPO terms associated with each variant
    with open(args.output + '_' + 'summary_stats' + '.tsv','w') as out:
        out.write('{pop}\t{num_variants}\t{mean_hpos}\t{median_hpos}\t{total}\n'.format(pop='pop', num_variants='num_variants', mean_hpos='mean_hpos', median_hpos='median_hpos',total='total_hpos'))
        for pop in pops:
            num_variants = len(pops[pop])
            num_hpos = []
            for line in pops[pop]:
                num_hpos.append(len(line.strip().split('\t')[6].split(',')))
            mean_hpos = sum(num_hpos)/len(num_hpos)
            median_hpos = sorted(num_hpos)[len(num_hpos)//2]
            # total HPO terms
            total = sum(num_hpos)
            out.write('{pop}\t{num_variants}\t{mean_hpos}\t{median_hpos}\t{total}\n'.format(pop=pop, num_variants=num_variants, mean_hpos=mean_hpos, median_hpos=median_hpos,total=total))


if __name__ == '__main__':
    main()