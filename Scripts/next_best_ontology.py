import pandas as pd

train = pd.read_csv('ELs_for_Rotate/Monarch_KG/train.txt', sep='\t', header=None)
valid = pd.read_csv('ELs_for_Rotate/Monarch_KG/valid.txt', sep='\t', header=None)
test = pd.read_csv('ELs_for_Rotate/Monarch_KG/test.txt', sep='\t', header=None)
df = pd.concat([train, valid, test], ignore_index=True)
df.columns = ['head', 'relation', 'tail']

GROUPS_2_path = {'Cancer':'work_comparison/cancer_genes.txt', 
                'PedCancer':'work_comparison/pediatric_cancer_genes.txt', 
                'European':'work_comparison/mondo_terms.nfe_onf.txt', 
                'EastAsian':'work_comparison/mondo_terms.eas.txt', 
                'Latino':'work_comparison/mondo_terms.amr.txt', 
                'African':'work_comparison/mondo_terms.afr.txt', 
                'Female':'work_comparison/female_differentially_expressed_genes.txt', 
                'Male':'work_comparison/male_differentially_expressed_genes.txt', 
                'RandomGenes':'work_comparison/random_500_genes.txt',
                'UltraRareDisease':'OrphanetEpidemiology/300_mondo_ids.txt',
                'RareDisease':'OrphanetEpidemiology/290_non_ultra_rare.mondo_ids.tsv',
                'RandomDiseases':'work_comparison/random_300_diseases.txt'}

group_2_terms = {}
for group, path in GROUPS_2_path.items():
    with open(path, 'r') as f:
        group_2_terms[group] = set([line.strip() for line in f if 'MONDO' in line or 'HGNC' in line])

# print each group and the number of terms
for group, terms in group_2_terms.items():
    print(group, len(terms))

# get a list of all node types
node_types = set([node.strip().split(':')[0] for node in df['head']])
node_types = node_types.union(set([node.strip().split(':')[0] for node in df['tail']]))

print('meta nodes',len(node_types))
print(node_types)

# remove MONDO and HGNC from metanodes
node_types = [node for node in node_types if node not in ['MONDO', 'HGNC']]

# for each group create a dataframe of all edges involving a group term

group_2_df = {}
for group, terms in group_2_terms.items():
    group_df = df[df['head'].isin(terms) | df['tail'].isin(terms)]
    group_2_df[group] = group_df

# for each meta node
metanode_counts_by_group = {'meta_node':[]}
for node_type in node_types:
    print(node_type)
    # for each group
    metanode_counts_by_group['meta_node'].append(node_type)
    for group, group_df in group_2_df.items():
        # get the number of edges involving the meta node
        num_edges = group_df[group_df['head'].str.startswith(node_type) | group_df['tail'].str.startswith(node_type)].shape[0]
        if group not in metanode_counts_by_group:
            metanode_counts_by_group[group] = []
        metanode_counts_by_group[group].append(num_edges)

metanode_counts_by_group_df = pd.DataFrame(metanode_counts_by_group)

# remove metanodes that are 0 for all groups, but report how many are moved
print('before', metanode_counts_by_group_df.shape)
# sum across all columns excet meta_node
tmp_df = metanode_counts_by_group_df.drop('meta_node', axis=1)
metanode_counts_by_group_df = metanode_counts_by_group_df[tmp_df.sum(axis=1) > 0]
print('after', metanode_counts_by_group_df.shape)

# rank metanodes by number of cancer edges added, high is bad
metanode_counts_by_group_df['inverse_cancer_rank'] = metanode_counts_by_group_df['Cancer'].rank(ascending=False)

# rank based on the number of Female Male, and RandomGenes edges added, low is bad
metanode_counts_by_group_df['non_cancer_genes_count'] = metanode_counts_by_group_df['Female'] + metanode_counts_by_group_df['Male'] + metanode_counts_by_group_df['RandomGenes']
metanode_counts_by_group_df['non_cancer_genes_rank'] = metanode_counts_by_group_df['non_cancer_genes_count'].rank(ascending=True)

# rank based on the number of UltraRareDisease, RareDisease, and RandomDiseases edges added, low is bad
metanode_counts_by_group_df['rare_disease_count'] = metanode_counts_by_group_df['UltraRareDisease'] + metanode_counts_by_group_df['RareDisease'] + metanode_counts_by_group_df['RandomDiseases']
metanode_counts_by_group_df['rare_disease_rank'] = metanode_counts_by_group_df['rare_disease_count'].rank(ascending=True)

# rank based on the number of European, EastAsian, Latino, and African edges added and random disease, low is bad
metanode_counts_by_group_df['asvd_count'] = metanode_counts_by_group_df['European'] + metanode_counts_by_group_df['EastAsian'] + metanode_counts_by_group_df['Latino'] + metanode_counts_by_group_df['African'] + metanode_counts_by_group_df['RandomDiseases']
metanode_counts_by_group_df['asvd_rank'] = metanode_counts_by_group_df['asvd_count'].rank(ascending=True)

# sum of ranks
metanode_counts_by_group_df['sum_rank'] = metanode_counts_by_group_df['inverse_cancer_rank'] + metanode_counts_by_group_df['non_cancer_genes_rank'] + metanode_counts_by_group_df['rare_disease_rank'] + metanode_counts_by_group_df['asvd_rank']

metanode_counts_by_group_df = metanode_counts_by_group_df.sort_values('sum_rank')
metanode_counts_by_group_df.to_csv('Data/metanode_counts_by_group_ranked.csv', index=False)

# print the top 10 by sum of rank
print(metanode_counts_by_group_df.head(10))

# print the bottom 10 by sum of rank
print(metanode_counts_by_group_df.tail(10))

# print the best for each group
for group in group_2_terms.keys():
    print(group, metanode_counts_by_group_df.sort_values(group).head(1))


