import random
import os

"""
The purpose of this pipeline is to download and process the Monarch KG data so that it is in the triple format needed for training KGE models with PyKeen
"""

dates = ["2022-09-27", "2023-03-16", "2023-09-28"]
types = ["edges", "nodes"]

meta_types = ['dictyBase','NCBIGene','PomBase','Xenbase','RGD','Reactome','MGI','SGD','FB','WB','ZFIN','NCIT','HSAPDV','NCBITaxon','SO','PATO','CL','CHR','CHEBI','UBERON','GO','HP']

rule all:
    input:
        "ELs_for_Rotate/Monarch_KG/test.txt",
        "ELs_for_Rotate/Monarch_KG/train.txt",
        "ELs_for_Rotate/Monarch_KG/valid.txt",
        "ELs_for_Rotate/Monarch_HuRI/train.txt",
        "ELs_for_Rotate/Monarch_HuRI/valid.txt",
        "ELs_for_Rotate/Monarch_HuRI/test.txt",
        "ELs_for_Rotate/Monarch_KG_Filtered/train.txt",
        "ELs_for_Rotate/Monarch_KG_Filtered/valid.txt",
        "ELs_for_Rotate/Monarch_KG_Filtered/test.txt",
        "ELs_for_Rotate/Monarch_HuRI_Filtered/train.txt",
        "ELs_for_Rotate/Monarch_HuRI_Filtered/valid.txt",
        "ELs_for_Rotate/Monarch_HuRI_Filtered/test.txt",
        expand('Metanode_ELs/{meta_type}/train.txt', meta_type=meta_types),
        expand('Metanode_ELs/{meta_type}/valid.txt', meta_type=meta_types),
        expand('Metanode_ELs/{meta_type}/test.txt', meta_type=meta_types),
        "Metanode_ELs/CC/train.txt",
        "Metanode_ELs/CC/valid.txt",
        "Metanode_ELs/CC/test.txt",
        "Metanode_ELs/CC/all.txt"


# # Commented out to avoid rerunning the creation of the Monarch KG, which creates a version slightly different that what I trained on and changes test set results.
# rule download:
#     output:
#         nodes="Monarch_KG/monarch-kg_nodes.2023-12-16.tsv",
#         edges="Monarch_KG/monarch-kg_edges.2023-12-16.tsv"
#     shell:
#         """
#         mkdir -p Monarch_KG/
#         wget https://data.monarchinitiative.org/monarch-kg-dev/2023-12-16/monarch-kg.tar.gz
#         tar -xf monarch-kg.tar.gz
#         mv monarch-kg_edges.tsv Monarch_KG/monarch-kg_edges.2023-12-16.tsv
#         mv monarch-kg_nodes.tsv Monarch_KG/monarch-kg_nodes.2023-12-16.tsv
#         rm monarch-kg.tar.gz
#         """

# rule create_triples:
#     input:
#         nodes="Monarch_KG/monarch-kg_nodes.2023-12-16.tsv",
#         edges="Monarch_KG/monarch-kg_edges.2023-12-16.tsv"
#     output:
#         'Monarch_KG/monarch-kg_triples.2023-12-16.tsv'
#     shell:
#         """
#         mkdir -p ELs_for_Rotate
#         mkdir -p ELs_for_Rotate/Monarch_KG

#         cat Monarch_KG/monarch-kg_edges.2023-12-16.tsv | cut -f3,18,19 | awk '{{print $2"\t"$1"\t"$3}}' | grep -v subject > {output}
#         """

# rule split_triples:
#     input:
#         "Monarch_KG/monarch-kg_triples.2023-12-16.tsv"
#     output:
#         train="ELs_for_Rotate/Monarch_KG/train.txt",
#         valid="ELs_for_Rotate/Monarch_KG/valid.txt",
#         test="ELs_for_Rotate/Monarch_KG/test.txt"
#     run:
#         # split the triples into train, valid, and test with .8, .1, .1 split
#         random.seed(42)
#         with open(output.train,'w') as train, open(output.valid,'w') as valid, open(output.test,'w') as test:
#             for line in open(input[0],'r'):
#                 rand = random.random()
#                 if rand < .8:
#                     train.write(line)
#                 elif rand < .9:
#                     valid.write(line)
#                 else:
#                     test.write(line)
    
rule map_HGNC_to_ENSG:
    input:
        "Monarch_KG/monarch-kg_nodes.2023-12-16.tsv",
    output:
        "Resources/HGNC_to_ENSG.tsv",
    run:
        hgnc_i = 0
        ensg_i = 3
        h2e = {}
        for line in open(input[0],'r'):
            row = line.strip().split('\t')
            ensgs = [x for x in row[ensg_i].split('|') if x.startswith('ENSEMBL')]
            if len(ensgs) == 0:
                continue
            h2e[row[hgnc_i]] = ensgs[0]
        with open(output[0],'w') as out:
            for hgnc in h2e:
                out.write(f"{hgnc}\t{h2e[hgnc]}\n")

rule convert_huri_to_hgnc:
    input:
        mapping="Resources/HGNC_to_ENSG.tsv",
        huri="Resources/HuRI.tsv"
    output:
        "Resources/huri.hgnc.triples.tsv"
    run:
        # load the mapping files as a dictionary but in reverse
        hgnc_to_ensg = {}
        for line in open(input.mapping,'r'):
            hgnc, ensg = line.strip().split('\t')
            hgnc_to_ensg[ensg] = hgnc
        # convert all huri triples to hgnc ids
        skipped_mappings = 0
        non_skipped = 0 
        with open(output[0],'w') as out:
            for line in open(input.huri,'r'):
                row = line.strip().split('\t')
                row[0] = 'ENSEMBL:' + row[0]
                row[1] = 'ENSEMBL:' + row[1]
                if row[0].startswith('ENSEMBL') and row[1].startswith('ENSEMBL'):
                    if row[0] in hgnc_to_ensg:
                        hgnc1 = hgnc_to_ensg[row[0]]
                        non_skipped += 1
                    else:
                        skipped_mappings += 1
                        continue
                    if row[1] in hgnc_to_ensg:
                        hgnc2 = hgnc_to_ensg[row[1]]
                        non_skipped += 1
                    else:
                        skipped_mappings += 1
                        continue
                    out.write(f"{hgnc1}\tbiolink:interacts_with\t{hgnc2}\n")
        print(f"------------------------------Skipped {skipped_mappings} mappings------------------------------")
        print(f"------------------------------Converted {non_skipped} mappings------------------------------")

rule make_monarch_HuRI:
    input:
        mapping="Resources/HGNC_to_ENSG.tsv",
        huri="Resources/huri.hgnc.triples.tsv",
        monarch="Monarch_KG/monarch-kg_triples.2023-12-16.tsv"
    output:
        monarch_huri="Resources/monarch_HuRI.triples.tsv",
    shell:
        """
        mkdir -p ELs_for_Rotate/Monarch_HuRI
        python Scripts/create_monarch_huri.py -a {input.monarch} -b {input.huri} -p 'HGNC' -o {output.monarch_huri}.tmp
        # remove duplicates
        cat {output.monarch_huri}.tmp | sort | uniq > {output.monarch_huri}
        rm {output.monarch_huri}.tmp
        """

rule split_monarch_HuRI:
    input:
        monarch_huri="Resources/monarch_HuRI.triples.tsv"
    output:
        train="ELs_for_Rotate/Monarch_HuRI/train.txt",
        valid="ELs_for_Rotate/Monarch_HuRI/valid.txt",
        test= "ELs_for_Rotate/Monarch_HuRI/test.txt"
    run:
        os.makedirs('ELs_for_Rotate/Monarch_HuRI', exist_ok=True)
        # split the triples into train, valid, and test with .8, .1, .1 split
        random.seed(42)
        with open(output.train,'w') as train, open(output.valid,'w') as valid, open(output.test,'w') as test:
            for line in open(input[0],'r'):
                rand = random.random()
                if rand < .8:
                    train.write(line)
                elif rand < .9:
                    valid.write(line)
                else:
                    test.write(line)

rule keep_only_hgnc_mondo:
    input:
        "Monarch_KG/monarch-kg_triples.2023-12-16.tsv"
    output:
        out="Monarch_KG/monarch-kg_triples.2023-12-16.filtered.tsv"
    run:
        # split the triples into train, valid, and test with .8, .1, .1 split
        random.seed(42)
        with open(output.out,'w') as out:
            for line in open(input[0],'r'):
                row = line.strip().split('\t')
                type1 = row[0].split(':')[0]
                type2 = row[2].split(':')[0]
                if type1 in ['HGNC','MONDO'] and type2 in ['HGNC','MONDO']:
                    out.write(line)
                    

rule make_monarch_HuRI_filtered:
    input:
        mapping="Resources/HGNC_to_ENSG.tsv",
        huri="Resources/huri.hgnc.triples.tsv",
        monarch="Monarch_KG/monarch-kg_triples.2023-12-16.filtered.tsv"
    output:
        out="Resources/monarch_HuRI_filtered.tsv"
    shell:
        """
        mkdir -p ELs_for_Rotate/Monarch_HuRI_Filtered
        python Scripts/create_monarch_huri.py -a {input.monarch} -b {input.huri} -p 'HGNC' -o {output.out}.tmp
        # remove duplicates
        cat {output.out}.tmp | sort | uniq > {output.out}
        rm {output.out}.tmp
        """

rule split_monarch_filtered_triples:
    input:
        "Monarch_KG/monarch-kg_triples.2023-12-16.filtered.tsv"
    output:
        train="ELs_for_Rotate/Monarch_KG_Filtered/train.txt",
        valid="ELs_for_Rotate/Monarch_KG_Filtered/valid.txt",
        test="ELs_for_Rotate/Monarch_KG_Filtered/test.txt"
    run:
        os.makedirs('ELs_for_Rotate/Monarch_KG_Filtered', exist_ok=True)
        # split the triples into train, valid, and test with .8, .1, .1 split
        random.seed(42)
        with open(output.train,'w') as train, open(output.valid,'w') as valid, open(output.test,'w') as test:
            for line in open(input[0],'r'):
                rand = random.random()
                if rand < .8:
                    train.write(line)
                elif rand < .9:
                    valid.write(line)
                else:
                    test.write(line)

rule split_huri_filtered_triples:
    input:
        "Resources/monarch_HuRI_filtered.tsv"
    output:
        train="ELs_for_Rotate/Monarch_HuRI_Filtered/train.txt",
        valid="ELs_for_Rotate/Monarch_HuRI_Filtered/valid.txt",
        test="ELs_for_Rotate/Monarch_HuRI_Filtered/test.txt"
    run:
        os.makedirs('ELs_for_Rotate/Monarch_HuRI_Filtered', exist_ok=True)
        # split the triples into train, valid, and test with .8, .1, .1 split
        random.seed(42)
        with open(output.train,'w') as train, open(output.valid,'w') as valid, open(output.test,'w') as test:
            for line in open(input[0],'r'):
                rand = random.random()
                if rand < .8:
                    train.write(line)
                elif rand < .9:
                    valid.write(line)
                else:
                    test.write(line)

# ------------------ Make metanode ELs ------------------
rule make_meta_els:
    input:
        "Monarch_KG/monarch-kg_triples.2023-12-16.tsv"
    output:
        train="Metanode_ELs/{meta_type}/train.txt",
        valid="Metanode_ELs/{meta_type}/valid.txt",
        test="Metanode_ELs/{meta_type}/test.txt",
        el_all="Metanode_ELs/{meta_type}/all.txt"
    run:
        # split the triples into train, valid, and test with .8, .1, .1 split
        os.makedirs(f'Metanode_ELs/{wildcards.meta_type}', exist_ok=True)
        random.seed(42)
        with open(output.el_all,'w') as out:
            for line in open(input[0],'r'):
                row = line.strip().split('\t')
                type1 = row[0].split(':')[0]
                type2 = row[2].split(':')[0]
                if type1 in ['HGNC','MONDO',wildcards.meta_type] and type2 in ['HGNC','MONDO',wildcards.meta_type]:
                    out.write(line)
        # remove duplicates from all
        os.system(f"cat {output.el_all} | sort | uniq > {output.el_all}.tmp")
        os.system(f"mv {output.el_all}.tmp {output.el_all}")
        with open(output.train,'w') as train, open(output.valid,'w') as valid, open(output.test,'w') as test:
            for line in open(output.el_all,'r'):
                rand = random.random()
                if rand < .8:
                    train.write(line)
                elif rand < .9:
                    valid.write(line)
                else:
                    test.write(line)

rule make_connected_components_graph:
    input:
        "Monarch_KG/monarch-kg_triples.2023-12-16.tsv"
    output:
        train="Metanode_ELs/CC/train.txt",
        valid="Metanode_ELs/CC/valid.txt",
        test="Metanode_ELs/CC/test.txt",
        el_all="Metanode_ELs/CC/all.txt"
    params:
        metanodes = ['HGNC','MONDO'] + meta_types
    run:
        # split the triples into train, valid, and test with .8, .1, .1 split
        os.makedirs(f'Metanode_ELs/CC', exist_ok=True)
        random.seed(42)
        with open(output.el_all,'w') as out:
            for line in open(input[0],'r'):
                row = line.strip().split('\t')
                type1 = row[0].split(':')[0]
                type2 = row[2].split(':')[0]
                if type1 in params.metanodes and type2 in params.metanodes:
                    out.write(line)
        # remove duplicates from all
        os.system(f"cat {output.el_all} | sort | uniq > {output.el_all}.tmp")
        os.system(f"mv {output.el_all}.tmp {output.el_all}")
        with open(output.train,'w') as train, open(output.valid,'w') as valid, open(output.test,'w') as test:
            for line in open(output.el_all,'r'):
                rand = random.random()
                if rand < .8:
                    train.write(line)
                elif rand < .9:
                    valid.write(line)
                else:
                    test.write(line)
        