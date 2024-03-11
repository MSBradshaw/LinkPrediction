import pandas as pd
import networkx as nx
import random
import sys
sys.path.append('Scripts/')
from compare_groups_test_omatic import plot_two_groups_hists, kruskal_test

POPs = ['nfe_onf','afr','amr','eas']

# '00', '01' ... '99'
SEX_SHARDS = [ '0'+str(i) if i < 10 else str(i) for i in range(0,100)]

# PPIs = ['original_monarch','HuRI','HuRI_filtered','monarch_filtered','string_filtered_t25','string_filtered_t50','string_filtered_t100']
# PPIs = ['HuRI']
# PPIs = ['HuRI_filtered','monarch_filtered']
# PPIs = ['monarch_filtered']
# PPIs = ['original_monarch']
PPIs = ['original_monarch','string_filtered_t25','string_filtered_t50','string_filtered_t100']
# Models = ['TransE','RotatE','ComplEx']
Models = ['TransE'] # so I can force run only TransE HuRI while I wait for the other two to finish running
# Models = ['RotatE']
# Models = ['ComplEx']

rule all:
    input:
        # expand('work_comparison/mondo_terms.{population}.txt',population=POPs),
        # expand('GroupComparisonResults/{PPI}/{model}/ASVD/euro_afr_gene_causes_mondo_European_v_African_g2p_rankings_hist.csv', PPI=PPIs, model=Models),
        # expand('GroupComparisonResults/{PPI}/{model}/ASVD/euro_latino_gene_causes_mondo_European_v_Latino_g2p_rankings_hist.csv', PPI=PPIs, model=Models),
        # expand('GroupComparisonResults/{PPI}/{model}/ASVD/euro_eas_gene_causes_mondo_European_v_East_Asian_g2p_rankings_hist.csv', PPI=PPIs, model=Models),
        # 'work_comparison/female_differentially_expressed_genes.txt',
        # 'work_comparison/male_differentially_expressed_genes.txt',
        # expand('GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.csv',PPI=PPIs, model=Models),
        # expand('GroupComparisonResults/{PPI}/{model}/SexDiffExpShards/sex_diff_genes_mondo_{i}_Female_v_Male_g2p_rankings_hist.csv',i=SEX_SHARDS, PPI=PPIs, model=Models),
        # 'work_comparison/cancer_genes.txt',
        # expand('GroupComparisonResults/{PPI}/{model}/CancerVsRandom/monarch_transE_Cancer_v_Random_500_42_g2p_rankings_hist.csv', PPI=PPIs, model=Models),
        # expand('GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.png', PPI=PPIs, model=Models),
        # expand('GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.txt', PPI=PPIs, model=Models),
        # expand('GroupComparisonResults/{PPI}/{model}/PediatricCancerVsRandom/pediatric_Pediatric_Cancer_v_Random_500_42_g2p_rankings_hist.csv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/Cancer/cancer.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/SexDiffExp/female_differentially_expressed_genes.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/SexDiffExp/male_differentially_expressed_genes.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/PediatricCancer/pediatric_cancer.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/ASVD/european.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/ASVD/african.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/ASVD/latino.tsv', PPI=PPIs, model=Models),
        expand('GroupRankResults/{PPI}/{model}/ASVD/east_asian.tsv', PPI=PPIs, model=Models)

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
        return 'Models/{}_monarch_string_filtered_t25_new.trained_model.pkl'.format(model)
    elif ppi == 'string_filtered_t50':
        return 'Models/{}_monarch_string_filtered_t50_new.trained_model.pkl'.format(model)
    elif ppi == 'string_filtered_t100':
        return 'Models/{}_monarch_string_filtered_t100_new.trained_model.pkl'.format(model)
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_datasplits(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.train.tsv', 'ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.valid.tsv', 'ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.test.tsv'
    elif ppi == 'original_monarch':
        return 'ELs_for_Rotate/Monarch_KG/train.txt', 'ELs_for_Rotate/Monarch_KG/valid.txt', 'ELs_for_Rotate/Monarch_KG/test.txt'
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
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_train(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.train.tsv'
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
    else:
        raise ValueError('Unknown PPI: ' + ppi)

def get_valid(ppi):
    if ppi == 'HuRI':
        return 'ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.valid.tsv'
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

rule do_asvd_comparisons:
    input:
        euro = 'work_comparison/mondo_terms.nfe_onf.txt',
        african = 'work_comparison/mondo_terms.afr.txt',
        latino = 'work_comparison/mondo_terms.amr.txt',
        east_asian = 'work_comparison/mondo_terms.eas.txt',
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt'
    output:
        'GroupComparisonResults/{PPI}/{model}/ASVD/euro_afr_gene_causes_mondo_European_v_African_g2p_rankings_hist.csv',
        'GroupComparisonResults/{PPI}/{model}/ASVD/euro_latino_gene_causes_mondo_European_v_Latino_g2p_rankings_hist.csv',
        'GroupComparisonResults/{PPI}/{model}/ASVD/euro_eas_gene_causes_mondo_European_v_East_Asian_g2p_rankings_hist.csv'
    threads: 2
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
        
    shell:
        """
        mkdir -p GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/ASVD/
        # EURO vs AFR
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.euro} \
            --b_terms {input.african} \
            --a_label European \
            --b_label African \
            --relation "biolink:causes" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.valid} \
            --test_triples {params.test} \
            --model {params.model} \
            --output_prefix GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/ASVD/euro_afr_gene_causes_mondo_

        # EURO vs AMR
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.euro} \
            --b_terms {input.latino} \
            --a_label European \
            --b_label Latino \
            --relation "biolink:causes" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.valid} \
            --test_triples {params.test} \
            --model {params.model} \
            --output_prefix GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/ASVD/euro_latino_gene_causes_mondo_
        
        # EURO vs EAS
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.euro} \
            --b_terms {input.east_asian} \
            --a_label European \
            --b_label East_Asian \
            --relation "biolink:causes" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.valid} \
            --test_triples {params.test} \
            --model {params.model} \
            --output_prefix GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/ASVD/euro_eas_gene_causes_mondo_
        """

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

rule shard_sex_lists:
    input:
        female='work_comparison/female_differentially_expressed_genes.txt',
        male='work_comparison/male_differentially_expressed_genes.txt',
    output:
        expand('work_comparison/SexDifChunks/female_differentially_expressed_genes_chunk_{i}.txt',i=SEX_SHARDS),
        expand('work_comparison/SexDifChunks/male_differentially_expressed_genes_chunk_{i}.txt',i=SEX_SHARDS),
    shell:
        """
        mkdir -p work_comparison/SexDifChunks/
        split -n l/100 -d --additional-suffix .txt {input.female} work_comparison/SexDifChunks/female_differentially_expressed_genes_chunk_
        split -n l/100 -d --additional-suffix .txt {input.male} work_comparison/SexDifChunks/male_differentially_expressed_genes_chunk_
        """

rule shards_sex_differentially_expressed_experiment:
    input:
        female='work_comparison/SexDifChunks/female_differentially_expressed_genes_chunk_{i}.txt',
        male='work_comparison/SexDifChunks/male_differentially_expressed_genes_chunk_{i}.txt',
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt'
    output:
        'GroupComparisonResults/{PPI}/{model}/SexDiffExpShards/sex_diff_genes_mondo_{i}_Female_v_Male_g2p_rankings_hist.csv'
    threads: 2
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
        grep "HGNC:" {input.male} > {input.male}.filtered

        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.female} \
            --b_terms {input.male} \
            --a_label Female \
            --b_label Male \
            --relation "biolink:interacts_with" \
            --prediction_target tail \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.valid} \
            --test_triples {params.test} \
            --model {params.model} \
            --output_prefix GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/SexDiffExpShards/sex_diff_genes_mondo_{wildcards.i}_ \
            --progress_bar
        """

rule combine_sex_shards:
    input:
        expand('GroupComparisonResults/{PPI}/{model}/SexDiffExpShards/sex_diff_genes_mondo_{i}_Female_v_Male_g2p_rankings_hist.csv',i=SEX_SHARDS, PPI='{PPI}', model='{model}'),
    output:
        'GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.csv'
    shell:
        """
        mkdir -p GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/SexDiffExp
        echo ',head_label,relation_label,tail_label,rank,head_degree,tail_degree,population' > {output}
        cat {input} | grep -v 'head_label' >> {output}
        """

rule sex_hist_and_kruskal:
    input:
        'GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.csv'
    output:
        'GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.png',
        'GroupComparisonResults/{PPI}/{model}/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.txt'
    run:
        df = pd.read_csv(input[0])
        female_ranks = list(df[df['population'] =='Female']['rank'])
        male_ranks = list(df[df['population'] =='Male']['rank'])
        plot_two_groups_hists(female_ranks,male_ranks,'Female','Male',f'GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/SexDiffExp/sex_diff_genes_mondo_')
        kruskal_test(female_ranks,male_ranks,'Female','Male',f'GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/SexDiffExp/sex_diff_genes_mondo_')

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

# cancer vs random
rule do_cancer_vs_random:
    input:
        cancer = 'work_comparison/cancer_genes.txt',
        random = 'work_comparison/random_500_genes.txt',
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt'
    output:
        'GroupComparisonResults/{PPI}/{model}/CancerVsRandom/monarch_transE_Cancer_v_Random_500_42_g2p_rankings_hist.csv'
    threads: 2
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/CancerVsRandom/
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.cancer} \
            --b_terms {input.random} \
            --a_label Cancer \
            --b_label Random_500_42 \
            --relation "biolink:interacts_with" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.valid} \
            --test_triples {params.test} \
            --model {params.model} \
            --output_prefix GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/CancerVsRandom/monarch_transE_ \
            --progress_bar
        """

rule make_peds_cancer_list:
    input:
        chco = 'Resources/chco_rna_fusions.csv',
        h2s = 'Resources/HGNC_to_symbol.tsv'
    output:
        'work_comparison/pediatric_cancer_genes.txt'
    run:
        h2s = {line.strip().split('\t')[1]:line.strip().split('\t')[0] for line in open(input.h2s,'r')}
        fiveprime = 2
        threeprime = 3
        written = set()
        with open(output[0],'w') as outfile:
            for line in open(input.chco,'r'):
                row = line.strip().split(',')
                if len(row) <= threeprime:
                    print('skipping',row)
                    continue
                if row[fiveprime] in h2s and row[fiveprime] not in written:
                    outfile.write(h2s[row[fiveprime]]+'\n')
                    written.add(h2s[row[fiveprime]])
                if row[threeprime] in h2s and row[threeprime] not in written:
                    outfile.write(h2s[row[threeprime]]+'\n')
                    written.add(h2s[row[threeprime]])

rule do_pediatric_cancer_vs_random:
    input:
        ped_cancer = 'work_comparison/pediatric_cancer_genes.txt',
        random = 'work_comparison/random_500_genes.txt',
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt'
    output:
        'GroupComparisonResults/{PPI}/{model}/PediatricCancerVsRandom/pediatric_Pediatric_Cancer_v_Random_500_42_g2p_rankings_hist.csv'
    threads: 2
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
        mkdir -p GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/PediatricCancerVsRandom/
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.ped_cancer} \
            --b_terms {input.random} \
            --a_label Pediatric_Cancer \
            --b_label Random_500_42 \
            --relation "biolink:interacts_with" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {params.train} \
            --validation_triples {params.valid} \
            --test_triples {params.test} \
            --model {params.model} \
            --output_prefix GroupComparisonResults/{wildcards.PPI}/{wildcards.model}/PediatricCancerVsRandom/pediatric_ \
            --progress_bar
        """

# ----------------------------------------- Generate Ranking Results -----------------------------------------

rule cancer_rank_results:
    input:
        terms = 'work_comparison/cancer_genes.txt'
    output:
        'GroupRankResults/{PPI}/{model}/Cancer/cancer.tsv'
    threads: 2
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
        'GroupRankResults/{PPI}/{model}/Random/random_500_42.tsv'
    threads: 2
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

rule pediatric_cancer_rank_results:
    input:
        terms = 'work_comparison/pediatric_cancer_genes.txt'
    output:
        'GroupRankResults/{PPI}/{model}/PediatricCancer/pediatric_cancer.tsv'
    threads: 2
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
        female='GroupRankResults/{PPI}/{model}/SexDiffExp/female_differentially_expressed_genes.tsv'
    threads: 2
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
        male='GroupRankResults/{PPI}/{model}/SexDiffExp/male_differentially_expressed_genes.tsv'
    threads: 2
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
        euro='GroupRankResults/{PPI}/{model}/ASVD/european.tsv',
        african='GroupRankResults/{PPI}/{model}/ASVD/african.tsv',
        latino='GroupRankResults/{PPI}/{model}/ASVD/latino.tsv',
        east_asian='GroupRankResults/{PPI}/{model}/ASVD/east_asian.tsv'
    threads: 2
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)        
    shell:
        """
        mkdir -p GroupRankResults/{wildcards.PPI}/{wildcards.model}/ASVD/
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

rule rare_disease_rank_results:
    input:
        'OrphanetEpidemiology/300_mondo_ids.txt'
    output:
        'GroupRankResults/{PPI}/{model}/RareDisease/rare_disease.tsv'
    threads: 2
    params:
        model = lambda wildcards: get_model(wildcards.PPI, wildcards.model),
        train = lambda wildcards: get_train(wildcards.PPI),
        valid = lambda wildcards: get_valid(wildcards.PPI),
        test = lambda wildcards: get_test(wildcards.PPI)
    shell:
        """
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