import pandas as pd
import networkx as nx

POPs = ['nfe_onf','afr','amr','eas']

rule all:
    input:
        expand('work_comparison/mondo_terms.{population}.txt',population=POPs),
        'GroupComparisonResults/ASVD/euro_afr_gene_causes_mondo_European_v_African_g2p_rankings_hist.csv',
        'GroupComparisonResults/ASVD/euro_latino_gene_causes_mondo_European_v_Latino_g2p_rankings_hist.csv',
        'GroupComparisonResults/ASVD/euro_eas_gene_causes_mondo_European_v_East_Asian_g2p_rankings_hist.csv'

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