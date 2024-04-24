import pandas as pd
import networkx as nx
import random
import sys

POPs = ['nfe_onf','afr','amr','eas']

# PPIs = ['original_monarch','HuRI']
PPIs = ['original_monarch','HuRI','monarch_filtered','string_filtered_t25','string_filtered_t50','string_filtered_t100','HuRI_filtered']

PPIs2ELDirName = {'original_monarch':'Monarch_KG',
                    'HuRI':'Monarch_HuRI',
                    'HuRI_filtered':'Monarch_HuRI_Filtered',
                    'monarch_filtered':'Monarch_KG_Filtered',
                    'string_filtered_t25':'Monarch_STRING_Filtered_t25_new',
                    'string_filtered_t50':'Monarch_STRING_Filtered_t50_new',
                    'string_filtered_t100':'Monarch_STRING_Filtered_t100_new'}

PPIs2pretty = {'original_monarch':'Monarch',
                    'HuRI':'HuRI',
                    'HuRI_filtered':'HuRI Filtered',
                    'monarch_filtered':'Monarch Filtered',
                    'string_filtered_t25':'String Filtered 25',
                    'string_filtered_t50':'String Filtered 50',
                    'string_filtered_t100':'String Filtered 100'}

Models = ['TransE','RotatE','ComplEx']

GROUPS = ['Cancer', 'PedCancer', 'European', 'EastAsian', 'Latino', 'African', 'Female', 'Male', 'RandomGenes','UltraRareDisease','RareDisease','RandomDiseases']
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

# diseases predict head, genes predict tail
GROUP_2_predictor_info = {'Cancer':['tail','biolink:interacts_with'], 
                        'PedCancer':['tail','biolink:interacts_with'], 
                        'European':['head','biolink:causes'], 
                        'EastAsian':['head','biolink:causes'], 
                        'Latino':['head','biolink:causes'], 
                        'African':['head','biolink:causes'], 
                        'Female':['tail','biolink:interacts_with'], 
                        'Male':['tail','biolink:interacts_with'], 
                        'RandomGenes':['tail','biolink:interacts_with'],
                        'UltraRareDisease':['head','biolink:causes'],
                        'RareDisease':['head','biolink:causes'],
                        'RandomDiseases':['head','biolink:causes']}

Filtered = ['normal','filtered']

rule all:
    input:
        expand('RankResults/{PPI}/{model}/{group}_ranking_data.tsv', PPI=PPIs, model=Models, group=GROUPS),
        expand('RankResultsHistograms/{PPI}.{model}.hist.png', PPI=PPIs, model=Models),
        expand('PercentTestAtK/{PPI}.{model}.percent_of_test_edges_at_k.png', PPI=PPIs, model=Models)

def load_ancestry_hpo_file(_f):
    _hpos = []
    for _line in open(_f,'r'):
        _row = _line.strip().split('\t')
        _tmp = _row[-1].split(',')
        _hpos += _tmp
    return list(set(_hpos))

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

european_terms = load_ancestry_hpo_file('AncestrySpecificMONDO/all_pathogenic_ancestry_groups__AF_nfe_onf.tsv')
african_terms = load_ancestry_hpo_file('AncestrySpecificMONDO/all_pathogenic_ancestry_groups__AF_afr.tsv')    
latino_terms = load_ancestry_hpo_file('AncestrySpecificMONDO/all_pathogenic_ancestry_groups__AF_amr.tsv')
east_asian_terms = load_ancestry_hpo_file('AncestrySpecificMONDO/all_pathogenic_ancestry_groups__AF_eas.tsv')

rule make_asvd_group_lists:
    input:
        'AncestrySpecificMONDO/all_pathogenic_ancestry_groups__AF_{population}.tsv',
        training_data = 'ELs_for_Rotate/Monarch_KG/train.txt',
    output:
        'work_comparison/mondo_terms.{population}.txt'
    run:
        df: pd.DataFrame = pd.read_csv(
            input[1], sep="\t", header=None, names=['source', 'relation', 'target']
        )
        # Load these edges into a NX graph and compute the degree for each entity
        H: nx.MultiGraph = nx.from_pandas_edgelist(
            df, "source", "target", create_using=nx.MultiGraph()
        )
        _hpos = load_ancestry_hpo_file(input[0])
        with open(output[0],'w') as _out:
            for _hpo in _hpos:
                if _hpo in H.nodes():
                    _out.write(_hpo+'\n')

rule make_sex_differentially_expressed_genes:
    input:
        'Resources/aba3066-table-s2.xlsx',
        'Resources/monarch-kg_nodes.sep_22.tsv',
        'ELs_for_Rotate/Monarch_KG/train.txt',
    output:
        female='work_comparison/female_differentially_expressed_genes.txt',
        male='work_comparison/male_differentially_expressed_genes.txt'
    shell:
        """ 
        python Scripts/create_sex_differentially_expressed_gene_list.py -m {output.male} -f {output.female}
        """

rule make_cancer_gene_list:
    output:
        'work_comparison/cancer_genes.txt'
    run:
        # Cancer genes
        cancer_genes = []
        #Breast cancer in women	
        cancer_genes = cancer_genes + ['ATM', 'BARD1', 'BRCA1', 'BRCA2', 'BRIP1', 'CHEK2',  'CDH1', 'NF1', 'PALB2', 'PTEN', 'RAD51C', 'RAD51D', 'STK11', 'TP53']

        # Breast cancer in men	
        cancer_genes = cancer_genes + ['BRCA1', 'BRCA2', 'CHEK2', 'PALB2']

        # Colorectal cancer	
        cancer_genes = cancer_genes + ['APC', 'EPCAM', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'CHEK2', 'PTEN', 'STK11', 'TP53', 'MUTYH']

        # Endometrial cancer	
        cancer_genes = cancer_genes + ['BRCA1', 'EPCAM', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'PTEN', 'STK11']

        # Fallopian tube, ovarian, primary peritoneal cancer	
        cancer_genes = cancer_genes + ['ATM', 'BRCA1', 'BRCA2', 'BRIP1', 'EPCAM', 'MLH1',  'MSH2',  'MSH6', 'NBN', 'PALB2', 'RAD51C', 'RAD51D']

        # Gastric cancer	
        cancer_genes = cancer_genes + ['APC', 'CDH1', 'STK11', 'EPCAM', 'MLH1', 'MSH2', 'MSH6', 'PMS2']

        # Melanoma	
        cancer_genes = cancer_genes + ['BAP1', 'BRCA2' 'CDK4', 'CDKN2A', 'PTEN', 'TP53']

        # Pancreatic cancer	
        cancer_genes = cancer_genes + ['ATM', 'BRCA1', 'BRCA2', 'CDKN2A', 'EPCAM', 'MLH1',  'MSH2', 'MSH6', 'PALB2', 'STK11', 'TP53']

        # Prostate cancer	
        cancer_genes = cancer_genes + ['ATM', 'BRCA1', 'BRCA2', 'CHEK2', 'HOXB13', 'PALB2', 'EPCAM', 'MLH1', 'MSH2', 'MSH6',  'PMS2']

        cancer_genes = list(set(cancer_genes))
        cancer_genes.sort()

        # convert symbol to hgnc
        gene_symbol_2_hgnc = {}
        with open('Resources/monarch-kg_nodes.sep_22.tsv','r') as f:
            for line in f:
                line = line.strip().split('\t')
                gene_symbol_2_hgnc[line[5]] = line[0]
        
        cancer_genes = [gene_symbol_2_hgnc[gene] for gene in cancer_genes if gene in gene_symbol_2_hgnc]

        with open(output[0],'w') as _out:
            for _gene in cancer_genes:
                _out.write(_gene+'\n')
  
rule make_random_gene_list:
    output:
        'work_comparison/random_500_genes.txt'
    run:
        # Random genes
        genes = []
        with open('Resources/monarch-kg_nodes.sep_22.tsv','r') as f:
            for line in f:
                line = line.strip().split('\t')
                if 'HGNC' in line[0]: 
                    genes.append(line[0])
                if 'HGNC' in line[2]: 
                    genes.append(line[2])
        genes = list(set(genes))
        genes.sort()
        # choose 500 genes at random from genes
        random.seed(42)
        random_genes = random.sample(genes, 500)

        with open(output[0],'w') as _out:
            for _gene in random_genes:
                _out.write(_gene+'\n')

rule make_random_disease_list:
    output:
        'work_comparison/random_300_diseases.txt'
    run:
        # Random diseases from all diseases with gene connections
        diseases = []
        for x in ['train','valid','test']:
            with open(f'ELs_for_Rotate/Monarch_KG/{x}.txt','r') as f:
                for line in f:
                    if 'MONDO:' in line and 'HGNC:' in line:
                        line = line.strip().split('\t')
                        if 'MONDO' in line[0]: 
                            diseases.append(line[0])
                        if 'MONDO' in line[2]: 
                            diseases.append(line[2])
        diseases = list(set(diseases))
        print(len(diseases))
        diseases.sort()
        # choose 500 genes at random from genes
        random.seed(42)
        random_diseases = random.sample(diseases, 300)

        with open(output[0],'w') as _out:
            for d in random_diseases:
                _out.write(d+'\n')

# ----------------------------------------- Generate Ranking Results -----------------------------------------




# cancer vs random
rule do_all:
    input:
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt'
    output:
        'RankResults/{PPI}/{model}/{group}_ranking_data.tsv'
    threads: 2
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        terms_list = lambda wildcards: GROUPS_2_path[wildcards.group],
        target = lambda wildcards: GROUP_2_predictor_info[wildcards.group][0],
        relation = lambda wildcards: GROUP_2_predictor_info[wildcards.group][1],
        train = lambda wildcards: get_el_split(wildcards.PPI,'train'),
        validation = lambda wildcards: get_el_split(wildcards.PPI,'valid'),
        test = lambda wildcards: get_el_split(wildcards.PPI,'test')
    shell:
        """
        mkdir -p RankResults/{wildcards.PPI}/{wildcards.model}/
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {params.terms_list} \
            --a_label {wildcards.group} \
            --relation {params.relation} \
            --prediction_target {params.target} \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.validation} \
            --test_triples {params.test} \
            --model {params.model} \
            --output {output} \
            --progress_bar
        """

# -------------------------------------------------------- Result Plotting --------------------------------------------------------

def count_group_occurances(group_terms,ppi_edge_list):
    count = 0
    group_set = set(group_terms)
    for line in open(ppi_edge_list,'r'):
        line = line.strip().split('\t')
        if line[0] in group_set or line[2] in group_set:
            count += 1
    return count

rule generate_stats_about_groups_and_ppi:
    input:
        cancer = 'work_comparison/cancer_genes.txt',
        random = 'work_comparison/random_500_genes.txt',
        random_diseases = 'work_comparison/random_300_diseases.txt',
        pedcancer = 'work_comparison/pediatric_cancer_genes.txt',
        female='work_comparison/female_differentially_expressed_genes.txt',
        male='work_comparison/male_differentially_expressed_genes.txt',
        african='work_comparison/mondo_terms.afr.txt',
        latino='work_comparison/mondo_terms.amr.txt',
        east_asian='work_comparison/mondo_terms.eas.txt',
        european='work_comparison/mondo_terms.nfe_onf.txt',
        ultra_rare='OrphanetEpidemiology/300_mondo_ids.txt',
        rare_disease='OrphanetEpidemiology/290_non_ultra_rare.mondo_ids.tsv'
    output:
        'work_comparison/stats_about_groups_and_ppi.tsv'
    params:
        ppis = PPIs
    run:
        groups2names = {
        'work_comparison/cancer_genes.txt':'Cancer',
        'work_comparison/random_500_genes.txt':'Random Genes',
        'work_comparison/random_300_diseases':'Random Diseases',
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
        for input_file in input:
            groups[input_file] = []
            for line in open(input_file,'r'):
                if 'HGNC' in line or 'MONDO' in line:
                    groups[input_file].append(line.strip())
        data = {'group':[],'KG':[],'total':[],'train':[],'valid':[],'test':[]}
        for group in groups:
            print(group)
            for ppi in params.ppis:
                data['group'].append(groups2names[group])
                data['KG'].append(ppi)
                data['train'].append(count_group_occurances(groups[group],get_train(ppi)))
                data['valid'].append(count_group_occurances(groups[group],get_valid(ppi)))
                data['test'].append(count_group_occurances(groups[group],get_test(ppi)))
                data['total'].append(data['train'][-1]+data['valid'][-1]+data['test'][-1])
        df = pd.DataFrame(data)
        df.to_csv(output[0],sep='\t',index=False)



rule plot_hist_by_model_and_kg:
    input: 
        expand('RankResults/{PPI}/{model}/{group}_ranking_data.tsv', PPI='{PPI}', model='{model}', group=GROUPS)
    output:
        'RankResultsHistograms/{PPI}.{model}.hist.png'
    params:
        title = lambda wc: f'{PPIs2pretty[wc.PPI]} {wc.model}'
    shell:
        """
        mkdir -p RankResultsHistograms
        python Scripts/plot_hist_within_model_ppi.py -i RankResults/{wildcards.PPI}/{wildcards.model}/ -o {output} -t "{params.title}"
        """

rule plot_percent_of_test_edge_at_k:
    input:
        expand('RankResults/{PPI}/{model}/{group}_ranking_data.tsv', PPI='{PPI}', model='{model}', group=GROUPS),
    output:
        'PercentTestAtK/{PPI}.{model}.percent_of_test_edges_at_k.png'
    params:
        title = lambda wc: f'{PPIs2pretty[wc.PPI]} {wc.model}'
    shell:
        """
        mkdir -p PercentTestAtK
        python Scripts/plot_percent_of_test_edge_at_k.py -i RankResults/{wildcards.PPI}/{wildcards.model}/ -o {output} -t "{params.title}"
        """
"""
rj -n 'alt_results' -T 64 -m 200G -c "snakemake -s run_group_comparison_experiments.smk --cores 64 --rerun-incomplete"
"""