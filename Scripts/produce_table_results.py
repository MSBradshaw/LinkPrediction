import pandas as pd
import matplotlib.pyplot as plt
import typing
import networkx as nx
import os
import pickle
import random
from pykeen.evaluation import evaluator_resolver
from pykeen.triples import TriplesFactory
import torch
from pykeen.metrics.ranking import HitsAtK
from pykeen.evaluation import RankBasedEvaluator
import json




def generate_results(model_path, test_path, grouping_path, output):
    loaded_model = torch.load(model_path,map_location=torch.device('cuda'))
    # load the grouping
    group = set()
    for line in open(grouping_path,'r'):
        group.add(line.strip())
    
    # get a random number 10000-99999
    rand = str(random.randint(10000,99999)) 
    tmp_test_edges = output + rand + ".txt"
    with open(tmp_test_edges,'w') as f:
        for line in open(test_path,'r'):
            line = line.strip().split('\t')
            if line[0] in group or line[2] in group:
                f.write('\t'.join(line)+'\n')


    test_triples_factory = TriplesFactory.from_path(tmp_test_edges)
    # test_triples_factory = TriplesFactory.from_path(test_path)
    
    # remove the tmp file
    os.remove(tmp_test_edges)

    # create a rankbased ev
    evaluator = evaluator_resolver.make(RankBasedEvaluator, metrics=[HitsAtK(k=25),HitsAtK(k=50),HitsAtK(k=100)], clear_on_finalize=False)
    evaluator.evaluate(model=loaded_model, mapped_triples=test_triples_factory.mapped_triples)
    res = evaluator.finalize().to_dict()
    with open(output, 'w') as f:
        for key in res.keys():
            res[key] = res[key]
            f.write(f"{key}\t{res[key]}\n")


def get_model(ppi,model):
    # fix string for model name
    model = model.lower()
    if model not in ['transe','rotate','complex']:
        raise ValueError('Unknown model: '+model)
        
    if ppi == 'HuRI':
        '{}_monarch_huri_filtered_3.trained_model.pkl'
        return 'Models/{}_monarch_huri.trained_model.pkl'.format(model)
    elif ppi == 'HuRI_filtered':
        return 'Models/{}_monarch_huri_filtered.trained_model.pkl'.format(model)
    elif ppi == 'monarch_filtered':
        return 'Models/{}_monarch_filtered.trained_model.pkl'.format(model)
    elif ppi == 'original_monarch':
        return 'Models/{}_monarch.trained_model.pkl'.format(model)
    elif ppi == 'string_filtered_t25':
        return 'Models/{}_monarch_string_filtered_t25_optimized.trained_model.pkl'.format(model)
    elif ppi == 'string_filtered_t50':
        return 'Models/{}_monarch_string_filtered_t50_optimized.trained_model.pkl'.format(model)
    elif ppi == 'string_filtered_t100':
        return 'Models/{}_monarch_string_filtered_t100_optimized.trained_model.pkl'.format(model)
    elif ppi == 'string_t25':
        return 'Models/{}_monarch_kg_string_t25.trained_model.pkl'.format(model)
    elif ppi == 'string_t50':
        return 'Models/{}_monarch_kg_string_t50.trained_model.pkl'.format(model)
    elif ppi == 'string_t100':
        return 'Models/{}_monarch_kg_string_t100.trained_model.pkl'.format(model)
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_test(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.test.tsv'
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
        return 'ELs_for_Rotate/Monarch_STRING_t25/test.txt'
    elif ppi == 'string_t50':
        return 'ELs_for_Rotate/Monarch_STRING_t50/test.txt'
    elif ppi == 'string_t100':
        return 'ELs_for_Rotate/Monarch_STRING_t100/test.txt'
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_group(group):
    mappings = {'Cancer':'work_comparison/cancer_genes.txt', 
     'PedCancer':'work_comparison/pediatric_cancer_genes.txt', 
     'European':'work_comparison/mondo_terms.nfe_onf.txt', 
     'EastAsian':'work_comparison/mondo_terms.eas.txt', 
     'Latino':'work_comparison/mondo_terms.amr.txt',
     'African':'work_comparison/mondo_terms.afr.txt',
     'Female':'work_comparison/female_differentially_expressed_genes.txt',
     'Male':'work_comparison/male_differentially_expressed_genes.txt',
     'Random':'work_comparison/random_500_genes.txt',
     'UltraRareDisease':'OrphanetEpidemiology/300_mondo_ids.txt',
     'RareDisease':'OrphanetEpidemiology/290_non_ultra_rare.mondo_ids.tsv',
     'RandomDiseases':'work_comparison/random_300_diseases.txt'}
    return mappings[group]

def load_result_file(file:str,head_or_tail:str,ppi:str,group:str,model:str) -> pd.DataFrame:
    print(file)
    results = {'hits_at_3':[],'hits_at_5':[],'hits_at_10':[],'hits_at_25':[],'hits_at_50':[],'hits_at_100':[],'inverse_harmonic_mean_rank':[],'group':[],'ppi':[],'model':[]}
    categories = ['hits_at_3','hits_at_5','hits_at_10','hits_at_25','hits_at_50','hits_at_100','inverse_harmonic_mean_rank']
    for line in open(file,'r'):
        line = line.strip().split('\t')
        if head_or_tail not in line[0]:
            continue
        tmp = json.loads(line[1].replace("'",'"'))
        for c in categories:
            results[c].append(tmp['realistic'][c])
        results['group'].append(group)
        results['ppi'].append(ppi)
        results['model'].append(model)
    return pd.DataFrame(results)

# PPIs = ['original_monarch','HuRI','HuRI_filtered','monarch_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100', 'string_t25', 'string_t50', 'string_t100'] # all
# PPIs = ['original_monarch','HuRI','HuRI_filtered','monarch_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100'] # all but t25 because rotate t25 is not ready yet
PPIs = ['original_monarch','monarch_filtered','HuRI_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100'] # all but t25 because rotate t25 is not ready yet
Models = ['TransE','RotatE','ComplEx']
Models = ['TransE','RotatE'] # no comple because I dont have/wont ever have string results for it
GROUPS = ['Cancer', 'PedCancer', 'European', 'EastAsian', 'Latino', 'African', 'Female', 'Male', 'Random','UltraRareDisease','RareDisease','RandomDiseases']

metrics = ['{head_or_tail}.realistic.inverse_harmonic_mean_rank']

# disease predict head, genes predict tail
group_to_head_tail  = {'Cancer':'tail', 'PedCancer':'tail', 'European':'head', 'EastAsian':'head', 'Latino':'head', 'African':'head', 'Female':'tail', 'Male':'tail', 'Random':'tail','UltraRareDisease':'head','RareDisease':'head','RandomDiseases':'head'}
ppi_to_pretty_name = {'original_monarch':'Monarch','HuRI':'HuRI','HuRI_filtered':'HuRI Filtered','monarch_filtered':'Monarch Filtered','string_filtered_t25':'String 25% Filtered', 'string_filtered_t50':'String 50% Filtered','string_filtered_t100':'String Filtered', 'string_t25':'String 25%', 'string_t50':'String 50%', 'string_t100':'String'}

all_results = None
for ppi in PPIs:
    test_el = get_test(ppi)
    for model in Models:
        model_path = get_model(ppi,model)
        for group in GROUPS:
            group_path = get_group(group)
            print('Doing combo:',ppi,model,group)
            output = f'TableResults/{ppi}.{model}.{group}.testing_results.txt'
            # check if output exists, if it does, skip
            if os.path.exists(output):
                print('\tAlready done, loading ')
                tmp_df = load_result_file(output,group_to_head_tail[group],ppi,group,model)
                if all_results is None:
                    all_results = tmp_df
                else:
                    all_results = pd.concat([all_results, tmp_df])
            else:
                print('\tProducing results')
                generate_results(model_path,test_el,group_path,output)

# write all results to file
all_results.to_csv('TableResults/all_results.tsv', sep='\t', index=False)

# seporate df by group
for group in GROUPS:
    tmp_df = all_results[all_results['group'] == group]
    tmp_df['ppi'] = tmp_df['ppi'].apply(lambda x: ppi_to_pretty_name[x])
    tmp_df.to_latex(f'TableResults/Latex/{group}.latex_table.txt',index=False)
    