import argparse
import pysam


# function for argparse
"""
input
1. tsv
2. ClinVar vcf
output tsv
"""

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help="tsv file with pathogenic ASVs")
    parser.add_argument('-c', '--clinvar', type=str, required=True, help="ClinVar vcf.gz file with accompanying tbi index")
    parser.add_argument('-o', '--output', type=str, required=True)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    vcf_file = pysam.VariantFile(args.clinvar, "r")
    columns = None
    with open(args.output,'w') as out:
        for line in open(args.input,'r'):
            row = line.strip().split('\t')
            if 'CHROM' in line:
                columns = row
                AF_col_idxs = [i for i,x in enumerate(columns) if 'AF_' in x]
                continue
            # if the column's value is not ? and > 0, it is the only population with any AF
            pop = [columns[i] for i in AF_col_idxs if row[i] != '?' and float(row[i]) > 0.00][0]
            # search vcf
            # get HPO terms
            # write to output
            target_pos = int(row[1])
            target_ref = row[2]
            target_alt = row[3]
            target_chrom = row[0]
            hpos = []
            for record in vcf_file.fetch(target_chrom, target_pos - 1, target_pos):
                if  record.pos == target_pos and record.ref == target_ref and target_alt in record.alts:
                    try:
                        terms = record.info['CLNDISDB']
                    except KeyError:
                        print('No CLNDISDB for record: {record}'.format(record=record.id))
                        print(str(record))
                        continue
                    # print(terms)
                    for section in terms:
                        if 'Human_Phenotype_Ontology' in section:
                            # use regex to get HPO terms that match this patter: 'HP:\d\d\d\d\d\d\d' in section
                            tmp_hpos = [term for term in section.split('|') if 'HP:' in term]
                            tmp_hpos = [ x.replace('Human_Phenotype_Ontology:','') for x in tmp_hpos]
                            if len(tmp_hpos) != 0:
                                hpos += tmp_hpos
                    if len(hpos) == 0:
                        continue
                    print(hpos)
                    template='{chr}\t{start}\t{ref}\t{alt}\t{cv_id}\t{pop}\t{hpos}\n'
                    new_line = template.format(chr=target_chrom, start=target_pos, ref=target_ref, alt=target_alt, cv_id=record.id, hpos=','.join(hpos), pop=pop)
                    out.write(new_line)
                    print()
                
                
            
                            

if __name__ == '__main__':
    main()