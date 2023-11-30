import pandas as pd
import networkx as nx
import random

POPs = ['nfe_onf','afr','amr','eas']

# '00', '01' ... '99'
SEX_SHARDS = [ '0'+str(i) if i < 10 else str(i) for i in range(0,100)]

rule all:
    input:
        expand('work_comparison/mondo_terms.{population}.txt',population=POPs),
        'GroupComparisonResults/ASVD/euro_afr_gene_causes_mondo_European_v_African_g2p_rankings_hist.csv',
        'GroupComparisonResults/ASVD/euro_latino_gene_causes_mondo_European_v_Latino_g2p_rankings_hist.csv',
        'GroupComparisonResults/ASVD/euro_eas_gene_causes_mondo_European_v_East_Asian_g2p_rankings_hist.csv',
        'work_comparison/female_differentially_expressed_genes.txt',
        'work_comparison/male_differentially_expressed_genes.txt',
        'GroupComparisonResults/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.csv',
        expand('GroupComparisonResults/SexDiffExpShards/sex_diff_genes_mondo_{i}_Female_v_Male_g2p_rankings_hist.csv',i=SEX_SHARDS),
        'work_comparison/cancer_genes.txt',
        'GroupComparisonResults/CancerVsRandom/monarch_transE_Cancer_v_Random_500_42_g2p_rankings_hist.csv'

def load_ancestry_hpo_file(_f):
    _hpos = []
    for _line in open(_f,'r'):
        _row = _line.strip().split('\t')
        _tmp = _row[-1].split(',')
        _hpos += _tmp
    return list(set(_hpos))

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
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt',
        model = 'Models/transE_monarch.pkl'
    output:
        'GroupComparisonResults/ASVD/euro_afr_gene_causes_mondo_European_v_African_g2p_rankings_hist.csv',
        'GroupComparisonResults/ASVD/euro_latino_gene_causes_mondo_European_v_Latino_g2p_rankings_hist.csv',
        'GroupComparisonResults/ASVD/euro_eas_gene_causes_mondo_European_v_East_Asian_g2p_rankings_hist.csv'
    threads: 2
    shell:
        """
        mkdir -p GroupComparisonResults/ASVD/
        # EURO vs AFR
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.euro} \
            --b_terms {input.african} \
            --a_label European \
            --b_label African \
            --relation "biolink:causes" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {input.train} \
            --validation_triples {input.validation} \
            --test_triples {input.test} \
            --model {input.model} \
            --output_prefix GroupComparisonResults/ASVD/euro_afr_gene_causes_mondo_

        # EURO vs AMR
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.euro} \
            --b_terms {input.latino} \
            --a_label European \
            --b_label Latino \
            --relation "biolink:causes" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {input.train} \
            --validation_triples {input.validation} \
            --test_triples {input.test} \
            --model {input.model} \
            --output_prefix GroupComparisonResults/ASVD/euro_latino_gene_causes_mondo_
        
        # EURO vs EAS
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.euro} \
            --b_terms {input.east_asian} \
            --a_label European \
            --b_label East_Asian \
            --relation "biolink:causes" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {input.train} \
            --validation_triples {input.validation} \
            --test_triples {input.test} \
            --model {input.model} \
            --output_prefix GroupComparisonResults/ASVD/euro_eas_gene_causes_mondo_
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

# rule sex_differentially_expressed_experiment:
#     input:
#         female='work_comparison/female_differentially_expressed_genes.txt',
#         male='work_comparison/male_differentially_expressed_genes.txt',
#         test =  'ELs_for_Rotate/Monarch_KG/test.txt',
#         train = 'ELs_for_Rotate/Monarch_KG/train.txt',
#         validation = 'ELs_for_Rotate/Monarch_KG/valid.txt',
#         model = 'Models/transE_monarch.pkl'
#     output:
#         'GroupComparisonResults/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.csv'
#     threads: 2
#     shell:
#         """
#         mkdir -p GroupComparisonResults/SexDiffExp/
#         python Scripts/compare_groups_test_omatic.py \
#             --a_terms {input.female} \
#             --b_terms {input.male} \
#             --a_label Female \
#             --b_label Male \
#             --relation "biolink:interacts_with" \
#             --prediction_target tail \
#             --prediction_prefix "HGNC:" \
#             --train_triples {input.train} \
#             --validation_triples {input.validation} \
#             --test_triples {input.test} \
#             --model {input.model} \
#             --output_prefix GroupComparisonResults/SexDiffExp/sex_diff_genes_mondo_ \
#             --progress_bar
#         """

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
        male='work_comparison/SexDifChunks/female_differentially_expressed_genes_chunk_{i}.txt',
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt',
        model = 'Models/transE_monarch.pkl'
    output:
        'GroupComparisonResults/SexDiffExpShards/sex_diff_genes_mondo_{i}_Female_v_Male_g2p_rankings_hist.csv'
    threads: 2
    shell:
        """
        mkdir -p GroupComparisonResults/SexDiffExpShards/
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.female} \
            --b_terms {input.male} \
            --a_label Female \
            --b_label Male \
            --relation "biolink:interacts_with" \
            --prediction_target tail \
            --prediction_prefix "HGNC:" \
            --train_triples {input.train} \
            --validation_triples {input.validation} \
            --test_triples {input.test} \
            --model {input.model} \
            --output_prefix GroupComparisonResults/SexDiffExpShards/sex_diff_genes_mondo_{wildcards.i}_ \
            --progress_bar
        """

rule combine_sex_shards:
    input:
        expand('GroupComparisonResults/SexDiffExpShards/sex_diff_genes_mondo_{i}_Female_v_Male_g2p_rankings_hist.csv',i=SEX_SHARDS)
    output:
        'GroupComparisonResults/SexDiffExp/sex_diff_genes_mondo_Female_v_Male_g2p_rankings_hist.csv'
    shell:
        """
        mkdir -p GroupComparisonResults/SexDiffExp
        echo ',head_label,relation_label,tail_label,rank,head_degree,tail_degree,population' > {output}
        cat {input} | grep -v 'head_label' >> {output}
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

# cancer vs random
rule do_cancer_vs_random:
    input:
        cancer = 'work_comparison/cancer_genes.txt',
        random = 'work_comparison/random_500_genes.txt',
        test =  'ELs_for_Rotate/Monarch_KG/test.txt',
        train = 'ELs_for_Rotate/Monarch_KG/train.txt',
        validation = 'ELs_for_Rotate/Monarch_KG/valid.txt',
        model = 'Models/transE_monarch.pkl'
    output:
        'GroupComparisonResults/CancerVsRandom/monarch_transE_Cancer_v_Random_500_42_g2p_rankings_hist.csv'
    threads: 2
    shell:
        """
        mkdir -p GroupComparisonResults/CancerVsRandom/
        python Scripts/compare_groups_test_omatic.py \
            --a_terms {input.cancer} \
            --b_terms {input.random} \
            --a_label Cancer \
            --b_label Random_500_42 \
            --relation "biolink:interacts_with" \
            --prediction_target head \
            --prediction_prefix "HGNC:" \
            --train_triples {input.train} \
            --validation_triples {input.validation} \
            --test_triples {input.test} \
            --model {input.model} \
            --output_prefix GroupComparisonResults/CancerVsRandom/monarch_transE_ \
            --progress_bar
        """