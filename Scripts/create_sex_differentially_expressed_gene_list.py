import pandas as pd
import argparse
import networkx as nx
# get args function, output female, output male
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--output_female',help='output file name',required=True)
    parser.add_argument('-m','--output_male',help='output file name',required=True)
    args = parser.parse_args()
    return args


# Load the data
xls = pd.ExcelFile('Resources/aba3066-table-s2.xlsx')

gene_tissue_counts_female = {}
gene_tissue_counts_male = {}
t2g2s = {'tissue':[],'gene':[],'sex':[]}
tissue_2_gene_by_sex = {}
for sheet_name in xls.sheet_names[2:]:
    xdf = xls.parse(sheet_name)
    # sort xdf by "MASH LFSR" in ascending order
    xdf = xdf.sort_values(by='MASH LFSR',ascending=True)
    # ensure the first row is smaller than the last row in the "MASH LSFR" column
    assert xdf['MASH LFSR'].iloc[0] < xdf['MASH LFSR'].iloc[-1]
    # get the top 100 genes
    xdf = xdf.head(100)
    for i,row in xdf.iterrows():
        g = row['HUGO_gene_id']
        is_female = row['MASH beta'] > 0
        if is_female:
            if g not in gene_tissue_counts_female:
                gene_tissue_counts_female[g] = 0
            gene_tissue_counts_female[g] += 1
            t2g2s['tissue'].append(sheet_name)
            t2g2s['gene'].append(g)
            t2g2s['sex'].append('female')
        else:
            if g not in gene_tissue_counts_male:
                gene_tissue_counts_male[g] = 0
            gene_tissue_counts_male[g] += 1
            t2g2s['tissue'].append(sheet_name)
            t2g2s['gene'].append(g)
            t2g2s['sex'].append('male')

# build a DF for these that are differentially expressed and with greater expression in females
genes_female = gene_tissue_counts_female.keys()
data_female = {'gene':genes_female,'tissue_count':[gene_tissue_counts_female[g] for g in genes_female]}
df_female = pd.DataFrame(data_female)
df_female['tissue_count'].max()
print('females',df_female[df_female['tissue_count'] == df_female['tissue_count'].max()])
df_female = df_female[df_female['tissue_count'] < 45]
print('females',df_female.shape)

# build a DF for these that are differentially expressed and with greater expression in males
genes_male = gene_tissue_counts_male.keys()
data_male = {'gene':genes_male,'tissue_count':[gene_tissue_counts_male[g] for g in genes_male]}
df_male = pd.DataFrame(data_male)
df_male['tissue_count'].max()
print('males',df_male[df_male['tissue_count'] == df_male['tissue_count'].max()])
df_male = df_male[df_male['tissue_count'] < 45]
print('males',df_male.shape)

genes_female = list(df_female['gene'])
genes_male = list(df_male['gene'])

# convert gene symbols to HGNC ids
# read Resources/monarch-kg_nodes.sep_22.tsv line by line and make a dictionary mapping symbol (5) to HGNC id (0)
gene_symbol_2_hgnc = {}
with open('Resources/monarch-kg_nodes.sep_22.tsv','r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene_symbol_2_hgnc[line[5]] = line[0]

# convert gene symbols to HGNC ids
genes_female = [ gene_symbol_2_hgnc[x] for x in genes_female if x in gene_symbol_2_hgnc] 
genes_male = [ gene_symbol_2_hgnc[x] for x in genes_male if x in gene_symbol_2_hgnc] 

# read edge list as a networkx graph
# ELs_for_Rotate/Monarch_KG/train.txt

df: pd.DataFrame = pd.read_csv(
    'ELs_for_Rotate/Monarch_KG/train.txt', sep="\t", header=None, names=['source', 'relation', 'target']
)
# Load these edges into a NX graph and compute the degree for each entity
G: nx.MultiGraph = nx.from_pandas_edgelist(
    df, "source", "target", create_using=nx.MultiGraph()
)

# remove genes that are not in the graph
genes_female = [x for x in genes_female if x in G.nodes()]
genes_male = [x for x in genes_male if x in G.nodes()]

args = get_args()
# write genes to files
with open(args.output_female,'w') as outfile:
    for g in genes_female:
        outfile.write(g+'\n')

# write genes to files
with open(args.output_male,'w') as outfile:
    for g in genes_male:
        outfile.write(g+'\n')
