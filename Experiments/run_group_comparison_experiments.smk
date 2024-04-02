import pandas as pd
import networkx as nx
import random
import sys
import matplotlib.pyplot as plt
sys.path.append('Scripts/')
from compare_groups_test_omatic import plot_two_groups_hists, kruskal_test

POPs = ['nfe_onf','afr','amr','eas']

# '00', '01' ... '99'
SEX_SHARDS = [ '0'+str(i) if i < 10 else str(i) for i in range(0,100)]

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

# PPIs = ['original_monarch','HuRI','HuRI_filtered','monarch_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100', 'string_t25', 'string_t50', 'string_t100'] # all
PPIs = ['original_monarch','HuRI','HuRI_filtered','monarch_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100'] # all but t25 because rotate t25 is not ready yet
PPIs = ['original_monarch','HuRI_filtered','monarch_filtered','string_filtered_t25', 'string_filtered_t50','string_filtered_t100'] # all but t25 because rotate t25 is not ready yet
Models = ['TransE','RotatE','ComplEx']
# Models = ['TransE','RotatE'] # no comple because I dont have/wont ever have string results for it
GROUPS = ['Cancer', 'PedCancer', 'European', 'EastAsian', 'Latino', 'African', 'Female', 'Male', 'Random','UltraRareDisease','RareDisease','RandomDiseases']

def generate_plotting_data(df:pd.DataFrame, G:nx.Graph, save_dir:str=None, filter_on_query_test_prescense=True) -> pd.DataFrame:
    # remove rows that were in the training set
    df['in_train'] = df.apply(lambda x: (x['head_label'], x['tail_label']) in G.edges, axis=1)
    df = df[~df['in_train']]
    # remove self edges
    df = df[df['head_label'] != df['tail_label']]
    # remove duplicate rows
    df = df.drop_duplicates(subset=['head_label', 'tail_label'])
    # sort by score in descending order
    df = df.sort_values(by='score', ascending=False)
    df['rank'] = df['score'].rank(pct=True)

    # find query terms that have atleast 1 test sest edge
    queries_to_keep = []
    for q in df['query_term'].unique():
        # where query term is q and in_test_set is True
        test_set_edges = df[(df['query_term'] == q) & (df['in_test_set'] == True)]
        if len(test_set_edges) > 0:
            queries_to_keep.append(q)
    print(queries_to_keep)

    queries_to_process = None

    if filter_on_query_test_prescense:
        if len(queries_to_keep) == 0:
            # return empty dataframes
            # throw an error and message
            # print('No query terms have test set edges')
            exit(1, 'No query terms have test set edges')
        queries_to_process = queries_to_keep
    else:
        queries_to_process = df['query_term'].unique()



    k_df = None # will contain all the hits at K data
    scored_df = None # will contain the data with ranks assigned by query term instead of globally
    for q in queries_to_process:
        sub = df[df['query_term'] == q]
        # sort by score in descending order
        df = df.sort_values(by='score', ascending=False)
        df['rank'] = df['score'].rank(pct=True)
            # create a dataframe that is K, hits@K, percent hits@K
        k_values = [1,5,10,15,25] + list(range(50, 501, 50))
        hits_at_k = []
        percent_hits_at_k = []
        for k in k_values:
            hits_at_k.append(sub.head(k)['in_test_set'].sum())
            percent_hits_at_k.append(sub.head(k)['in_test_set'].sum() / k)
        sub_k_df = pd.DataFrame({'k': k_values, 'hits@k': hits_at_k, 'percent_hits@k': percent_hits_at_k})
        sub_k_df['query_term'] = q
        if k_df is None:
            k_df = sub_k_df
        else:
            k_df = pd.concat([k_df, sub_k_df])
        if scored_df is None:
            scored_df = sub
        else:
            scored_df = pd.concat([scored_df, sub])
    
    print(k_df.head())
    print(k_df.shape)
    # create an average hits@k for each k
    # get k_df without the query term column
    k_df_no_q = k_df.drop(columns=['query_term'])
    avg_k_df = k_df_no_q.groupby('k').mean().reset_index()

    if save_dir:
        if save_dir[-1] != '/':
            save_dir += '/'
        os.makedirs(save_dir, exist_ok=True)
        avg_k_df.to_csv(save_dir + 'avg_hits_at_k.tsv', index=False, sep='\t')
        scored_df.to_csv(save_dir + 'edges_scored_by_query_term.tsv', index=False, sep='\t')
        k_df.to_csv(save_dir + 'all_hits_at_k.tsv', index=False, sep='\t')

    return avg_k_df, scored_df, k_df

rule all:
    input:
        expand('GroupRankResults/{PPI}/{model}/AllResults/{group}.comparison_results.tsv', PPI=PPIs, model=Models, group=GROUPS),
        expand('RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k.tsv', PPI=PPIs, model=Models, group=GROUPS),
        expand('HitsAtKurvePlots/kurve.{PPI}.{model}.png', PPI=PPIs, model=Models),
        'RankedResultsIntermediate/combined_avg_hits_at_k.tsv',
        'RankedResultsIntermediate/combined_avg_hits_at_k.filtered.tsv',
        'work_comparison/stats_about_groups_and_ppi.tsv'

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
        # Random genes
        diseases = []
        with open('Resources/monarch-kg_nodes.sep_22.tsv','r') as f:
            for line in f:
                line = line.strip().split('\t')
                if 'MONDO' in line[0]: 
                    diseases.append(line[0])
                if 'MONDO' in line[2]: 
                    diseases.append(line[2])
        diseases = list(set(diseases))
        diseases.sort()
        # choose 500 genes at random from genes
        random.seed(42)
        random_diseases = random.sample(diseases, 300)

        with open(output[0],'w') as _out:
            for d in random_diseases:
                _out.write(d+'\n')

# ----------------------------------------- Generate Ranking Results -----------------------------------------

rule cancer_rank_results:
    input:
        terms = 'work_comparison/cancer_genes.txt'
    output:
        'GroupRankResults/{PPI}/{model}/AllResults/Cancer.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/Cancer/
        python Scripts/generate_ranked_results.py \
                --terms {input.terms} \
                --label Cancer \
                --relation "biolink:interacts_with" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output} \
                --progress_bar
        """

rule random_rank_results:
    input:
        terms = 'work_comparison/random_500_genes.txt'
    output:
        'GroupRankResults/{PPI}/{model}/AllResults/Random.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/Random/
        python Scripts/generate_ranked_results.py \
                --terms {input.terms} \
                --label Random_500_42 \
                --relation "biolink:interacts_with" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output} \
                --progress_bar
        """

rule random_disease_rank_results:
    input:
        terms = 'work_comparison/random_300_diseases.txt'
    output:
        'GroupRankResults/{PPI}/{model}/AllResults/RandomDiseases.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/Random/
        python Scripts/generate_ranked_results.py \
                --terms {input.terms} \
                --label Random_Diseases_500_42 \
                --relation "biolink:interacts_with" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output} \
                --progress_bar
        """

rule pediatric_cancer_rank_results:
    input:
        terms = 'work_comparison/pediatric_cancer_genes.txt'
    output:
        'GroupRankResults/{PPI}/{model}/AllResults/PedCancer.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/PediatricCancer/
        python Scripts/generate_ranked_results.py \
                --terms {input.terms} \
                --label Pediatric_Cancer \
                --relation "biolink:interacts_with" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output} \
                --progress_bar
        """

rule female_sex_differentially_expressed_ranked_results:
    input:
        female='work_comparison/female_differentially_expressed_genes.txt'
    output:
        female='GroupRankResults/{PPI}/{model}/AllResults/Female.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/SexDiffExpShards/
        # remove non human terms if they made it in somehow
        grep "HGNC:" {input.female} > {input.female}.filtered

        python Scripts/generate_ranked_results.py \
                --terms {input.female} \
                --label Female \
                --relation "biolink:interacts_with" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output.female} \
                --progress_bar
        """

rule male_sex_differentially_expressed_ranked_results:
    input:
        male='work_comparison/male_differentially_expressed_genes.txt'
    output:
        male='GroupRankResults/{PPI}/{model}/AllResults/Male.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/SexDiffExpShards/
        # remove non human terms if they made it in somehow
        grep "HGNC:" {input.male} > {input.male}.filtered

        python Scripts/generate_ranked_results.py \
                --terms {input.male} \
                --label Male \
                --relation "biolink:interacts_with" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output.male} \
                --progress_bar
        """

                
rule asvd_ranked_results:
    input:
        euro = 'work_comparison/mondo_terms.nfe_onf.txt',
        african = 'work_comparison/mondo_terms.afr.txt',
        latino = 'work_comparison/mondo_terms.amr.txt',
        east_asian = 'work_comparison/mondo_terms.eas.txt'
    output:
        euro='GroupRankResults/{PPI}/{model}/AllResults/European.comparison_results.tsv',
        african='GroupRankResults/{PPI}/{model}/AllResults/African.comparison_results.tsv',
        latino='GroupRankResults/{PPI}/{model}/AllResults/Latino.comparison_results.tsv',
        east_asian='GroupRankResults/{PPI}/{model}/AllResults/EastAsian.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)        
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/AllResults/
        # EURO
        python Scripts/generate_ranked_results.py \
                --terms {input.euro} \
                --label European \
                --relation "biolink:causes" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output.euro} \
                --progress_bar
        
        # AFRICAN
        python Scripts/generate_ranked_results.py \
                --terms {input.african} \
                --label African \
                --relation "biolink:causes" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output.african} \
                --progress_bar
        
        # LATINO
        python Scripts/generate_ranked_results.py \
                --terms {input.latino} \
                --label Latino \
                --relation "biolink:causes" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output.latino} \
                --progress_bar
        
        # EAST ASIAN
        python Scripts/generate_ranked_results.py \
                --terms {input.east_asian} \
                --label East_Asian \
                --relation "biolink:causes" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output.east_asian} \
                --progress_bar
        """

rule ultra_rare_disease_rank_results:
    input:
        'OrphanetEpidemiology/300_mondo_ids.txt'
    output:
        'GroupRankResults/{PPI}/{model}/AllResults/UltraRareDisease.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/AllResults/
        python Scripts/generate_ranked_results.py \
                --terms {input} \
                --label Rare_Diseases \
                --relation "biolink:causes" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output} \
                --progress_bar
        """

rule rare_disease_rank_results:
    input:
        'OrphanetEpidemiology/290_non_ultra_rare.mondo_ids.tsv'
    output:
        'GroupRankResults/{PPI}/{model}/AllResults/RareDisease.comparison_results.tsv'
    threads: 1
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/AllResults/
        python Scripts/generate_ranked_results.py \
                --terms {input} \
                --label Rare_Diseases \
                --relation "biolink:causes" \
                --prediction_target head \
                --prediction_prefix "HGNC:" \
                --train_triples {params.train} \
                --validation_triples {params.valid} \
                --test_triples {params.test} \
                --model {params.model} \
                --output {output} \
                --progress_bar
        """

def get_el_as_G(ppi):
    train_path = get_train(ppi)
    valid_path = get_valid(ppi)
    el_df = pd.read_csv(train_path, sep='\t', header=None)
    el_df_valid = pd.read_csv(valid_path, sep='\t', header=None)
    el_df.columns = ['subject', 'predicate', 'object']
    el_df_valid.columns = ['subject', 'predicate', 'object']
    el_df = pd.concat([el_df, el_df_valid])

    G = nx.from_pandas_edgelist(el_df, 'subject', 'object', 'predicate')
    return G

rule generate_hits_at_k_plotting_data:
    input:
        'GroupRankResults/{PPI}/{model}/AllResults/{group}.comparison_results.tsv'
    output:
        avg_hits_at_k = 'RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k.tsv',
        edges_scored_by_query_term = 'RankedResultsIntermediate/{PPI}/{model}/{group}/edges_scored_by_query_term.tsv',
        all_hits_at_k = 'RankedResultsIntermediate/{PPI}/{model}/{group}/all_hits_at_k.tsv',
        f_avg_hits_at_k = 'RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k.filtered.tsv',
        f_edges_scored_by_query_term = 'RankedResultsIntermediate/{PPI}/{model}/{group}/edges_scored_by_query_term.filtered.tsv',
        f_all_hits_at_k = 'RankedResultsIntermediate/{PPI}/{model}/{group}/all_hits_at_k.filtered.tsv'
    params:
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI)
    shell:
        """
        python Scripts/generate_kurve_plotting_data.py --train_path {params.train} \
                                                    --valid_path {params.valid} \
                                                    --save_dir RankedResultsIntermediate/{wildcards.PPI}/{wildcards.model}/{wildcards.group}/ \
                                                    --df_path {input} \
                                                    --filter

        python Scripts/generate_kurve_plotting_data.py --train_path {params.train} \
                                                    --valid_path {params.valid} \
                                                    --save_dir RankedResultsIntermediate/{wildcards.PPI}/{wildcards.model}/{wildcards.group}/ \
                                                    --df_path {input}
        """

rule plot_hits_at_kurve:
    input:
        avg_hits_at_k = expand('RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k{filtered}.tsv', PPI='{PPI}', model='{model}', group=GROUPS, filtered=['', '.filtered']),
    output:
        'HitsAtKurvePlots/kurve.{PPI}.{model}{filtered}.png'
    params:
        ppi_name = None
    shell:
        """
        mkdir -p HitsAtKurvePlots
        python Scripts/plot_hits_at_kurve.py --ppi_model_dir RankedResultsIntermediate/{wildcards.PPI}/{wildcards.model}/ \
                                    --title {wildcards.PPI} \
                                    --outfile {output}
        """

rule annotate_hits_at_kurve:
    input:
        avg_hits_at_k = 'RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k{filtered}.tsv'
    output:
        'RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k.annotated{filtered}.tsv'
    run:
        df = pd.read_csv(input[0], sep='\t')
        df['group'] = wildcards.group
        df['model'] = wildcards.model
        df['ppi'] = wildcards.PPI
        df.to_csv(output[0], sep='\t', index=False)

rule combine_hits_at_kurve:
    input:
        expand('RankedResultsIntermediate/{PPI}/{model}/{group}/avg_hits_at_k.annotated{filtered}.tsv', PPI=PPIs, model=Models, group=GROUPS, filtered='{filtered}')
    output:
        'RankedResultsIntermediate/combined_avg_hits_at_k{filtered}.tsv'
    run:
        df = pd.concat([pd.read_csv(f, sep='\t') for f in input])
        df.to_csv(output[0], sep='\t', index=False)

rule combine_edges_scored_by_query_term:
    input:
        expand('RankedResultsIntermediate/{PPI}/{model}/{group}/edges_scored_by_query_term.tsv', PPI=PPIs, model=Models, group=GROUPS)
    output:
        'RankedResultsIntermediate/combined_edges_scored_by_query_term.tsv'
    run:
        combined_df = None
        for f in input:
            df = pd.read_csv(f, sep='\t')
            df['model'] = f.split('/')[2]
            df['ppi'] = f.split('/')[1]
            df['group'] = f.split('/')[3]
            if combined_df == None:
                combined_df = df
            else:
                combined_df = pd.concat([combined_df, df])
        combined_df.to_csv(output[0], sep='\t', index=False)


            

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
        for input_file in inputs:
            groups[input_file] = []
            for line in open(input_file,'r'):
                if 'HGNC' in line or 'MONDO' in line:
                    groups[input_file].append(line.strip())
        data = {'group':[],'KG':[],'total':[],'train':[],'valid':[],'test':[]}
        for group in groups:
            print(group)
            for ppi in params.ppis:
                data['group'].append(groups2names[group])
                data['KG'].append(groups2names[ppi])
                data['train'].append(count_group_occurances(groups[group],get_train(ppi)))
                data['valid'].append(count_group_occurances(groups[group],get_valid(ppi)))
                data['test'].append(count_group_occurances(groups[group],get_test(ppi)))
                data['total'].append(data['train'][-1]+data['valid'][-1]+data['test'][-1])
        df = pd.DataFrame(data)
        df.to_csv(output[0],sep='\t',index=False)
            