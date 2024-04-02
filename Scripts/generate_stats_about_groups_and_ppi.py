import pandas as pd

def count_group_occurances(group_terms,ppi_edge_list):
    count = 0
    group_set = set(group_terms)
    for line in open(ppi_edge_list,'r'):
        line = line.strip().split('\t')
        if line[0] in group_set or line[2] in group_set:
            count += 1
    return count

def get_test(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/test.txt'
    elif ppi == 'HuRI_filtered':
        return 'ELs_for_Rotate/Monarch_HuRI_Filtered/test.txt'
    elif ppi == 'monarch_filtered':
        return 'ELs_for_Rotate/Monarch_KG_Filtered/test.txt'
    elif ppi == 'original_monarch':
        return 'ELs_for_Rotate/Monarch_KG/test.txt'
    elif ppi == 'string_filtered_t25':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t25_new/test.txt'
    elif ppi == 'string_filtered_t50':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t50_new/test.txt'
    elif ppi == 'string_filtered_t100':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t100_new/test.txt'
    elif ppi == 'string_t25':
        return 'ELs_for_Rotate/Monarch_STRING_t25_new/test.txt'
    elif ppi == 'string_t50':
        return 'ELs_for_Rotate/Monarch_STRING_t50_new/test.txt'
    elif ppi == 'string_t100':
        return 'ELs_for_Rotate/Monarch_STRING_t100_new/test.txt'
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_train(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/train.txt'
    elif ppi == 'HuRI_filtered':
        return 'ELs_for_Rotate/Monarch_HuRI_Filtered/train.txt'
    elif ppi == 'monarch_filtered':
        return 'ELs_for_Rotate/Monarch_KG_Filtered/train.txt'
    elif ppi == 'original_monarch':
        return 'ELs_for_Rotate/Monarch_KG/train.txt'
    elif ppi == 'string_filtered_t25':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t25_new/train.txt'
    elif ppi == 'string_filtered_t50':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t50_new/train.txt'
    elif ppi == 'string_filtered_t100':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t100_new/train.txt'
    elif ppi == 'string_t25':
        return 'ELs_for_Rotate/Monarch_STRING_t25_new/train.txt'
    elif ppi == 'string_t50':
        return 'ELs_for_Rotate/Monarch_STRING_t50_new/train.txt'
    elif ppi == 'string_t100':
        return 'ELs_for_Rotate/Monarch_STRING_t100_new/train.txt'
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_valid(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/valid.txt'
    elif ppi == 'HuRI_filtered':
        return 'ELs_for_Rotate/Monarch_HuRI_Filtered/valid.txt'
    elif ppi == 'monarch_filtered':
        return 'ELs_for_Rotate/Monarch_KG_Filtered/valid.txt'
    elif ppi == 'original_monarch':
        return 'ELs_for_Rotate/Monarch_KG/valid.txt'
    elif ppi == 'string_filtered_t25':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t25_new/valid.txt'
    elif ppi == 'string_filtered_t50':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t50_new/valid.txt'
    elif ppi == 'string_filtered_t100':
        return 'ELs_for_Rotate/Monarch_STRING_Filtered_t100_new/valid.txt'
    elif ppi == 'string_t25':
        return 'ELs_for_Rotate/Monarch_STRING_t25_new/valid.txt'
    elif ppi == 'string_t50':
        return 'ELs_for_Rotate/Monarch_STRING_t50_new/valid.txt'
    elif ppi == 'string_t100':
        return 'ELs_for_Rotate/Monarch_STRING_t100_new/valid.txt'
    else:
        raise ValueError('Unknown PPI: ' + ppi)

PPIs = ['original_monarch','HuRI','HuRI_filtered','monarch_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100', 'string_t25', 'string_t50', 'string_t100'] # all

ppi_names = {'original_monarch':'Monarch',
            'HuRI':'HuRI',
            'HuRI_filtered':'Filtered HuRI',
            'monarch_filtered':'Filtered Monarch',
            'string_filtered_t25':'Filtered STRING top 25%',
            'string_filtered_t50':'Filtered STRING top 50%',
            'string_filtered_t100':'Filtered STRING top 100%',
            'string_t25':'STRING top 25%',
            'string_t50':'STRING top 50%',
            'string_t100':'STRING top 100%'}
groups2names = {
'work_comparison/cancer_genes.txt':'Cancer',
'work_comparison/random_500_genes.txt':'Random Genes',
'work_comparison/random_300_diseases.txt':'Random Diseases',
'work_comparison/pediatric_cancer_genes.txt':'Pediatric Cancer',
'work_comparison/female_differentially_expressed_genes.txt':'Female',
'work_comparison/male_differentially_expressed_genes.txt':'Male',
'work_comparison/mondo_terms.afr.txt':'African',
'work_comparison/mondo_terms.amr.txt':'Latino',
'work_comparison/mondo_terms.eas.txt':'East Asian',
'work_comparison/mondo_terms.nfe_onf.txt':'European',
'OrphanetEpidemiology/300_mondo_ids.txt':'Ultra Rare Diseases',
'OrphanetEpidemiology/290_non_ultra_rare.mondo_ids.tsv':'Rare Diseases'
}
# for each input keep only the HGNC and MONDO terms in them
groups = {}
for input_file in groups2names.keys():
    groups[input_file] = []
    for line in open(input_file,'r'):
        if 'HGNC' in line or 'MONDO' in line:
            groups[input_file].append(line.strip())
data = {'group':[],'KG':[],'total':[],'train':[],'valid':[],'test':[]}
for group in groups:
    print(group)
    for ppi in PPIs:
        data['group'].append(groups2names[group])
        data['KG'].append(ppi_names[ppi])
        data['train'].append(count_group_occurances(groups[group],get_train(ppi)))
        data['valid'].append(count_group_occurances(groups[group],get_valid(ppi)))
        data['test'].append(count_group_occurances(groups[group],get_test(ppi)))
        data['total'].append(data['train'][-1]+data['valid'][-1]+data['test'][-1])
df = pd.DataFrame(data)
df.to_csv('work_comparison/stats_about_groups_and_ppi.tsv',sep='\t',index=False)
    