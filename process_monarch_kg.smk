"""
The purpose of this pipeline is to download and process the Monarch KG data so that it is in the triple format needed for training KGE models with PyKeen
"""

dates = ["2022-09-27", "2023-03-16", "2023-09-28"]
types = ["edges", "nodes"]


rule all:
    input:
        "ELs_for_Rotate/Monarch_KG/test.txt",
        "ELs_for_Rotate/Monarch_KG/train.txt",
        "ELs_for_Rotate/Monarch_KG/valid.txt",
        "ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.train.tsv",
        "ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.valid.tsv",
        "ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.test.tsv"

rule download:
    output:
        expand("Monarch_KG/monarch-kg_{type}.{date}.tsv", type=types, date=dates)
    shell:
        """
        mkdir -p Monarch_KG/
        wget https://data.monarchinitiative.org/monarch-kg-dev/2022-09-27/monarch-kg.tar.gz
        tar -xf monarch-kg.tar.gz
        mv monarch-kg_edges.tsv Monarch_KG/monarch-kg_edges.2022-09-27.tsv
        mv monarch-kg_nodes.tsv Monarch_KG/monarch-kg_nodes.2022-09-27.tsv
        rm monarch-kg.tar.gz
        
        wget https://data.monarchinitiative.org/monarch-kg-dev/2023-03-16/monarch-kg.tar.gz
        tar -xf monarch-kg.tar.gz
        mv monarch-kg_edges.tsv Monarch_KG/monarch-kg_edges.2023-03-16.tsv
        mv monarch-kg_nodes.tsv Monarch_KG/monarch-kg_nodes.2023-03-16.tsv
        rm monarch-kg.tar.gz
        

        wget https://data.monarchinitiative.org/monarch-kg-dev/2023-09-28/monarch-kg.tar.gz
        tar -xf monarch-kg.tar.gz
        mv monarch-kg_edges.tsv Monarch_KG/monarch-kg_edges.2023-09-28.tsv
        mv monarch-kg_nodes.tsv Monarch_KG/monarch-kg_nodes.2023-09-28.tsv
        rm monarch-kg.tar.gz

        wget https://data.monarchinitiative.org/monarch-kg-dev/2022-03-07/monarch-kg.tar.gz
        tar -xf monarch-kg.tar.gz
        mv monarch-kg_edges.tsv Monarch_KG/monarch-kg_edges.2022-03-07.tsv
        mv monarch-kg_nodes.tsv Monarch_KG/monarch-kg_nodes.2022-03-07.tsv
        rm monarch-kg.tar.gz
        """

rule remove_nodes_and_edges_not_in_first:
    input:
        train_nodes="Monarch_KG/monarch-kg_nodes.2022-09-27.tsv",
        train="Monarch_KG/monarch-kg_edges.2022-09-27.tsv",
        valid="Monarch_KG/monarch-kg_edges.2023-03-16.tsv",
        test="Monarch_KG/monarch-kg_edges.2023-09-28.tsv"
    output:
        train="ELs_for_Rotate/Monarch_KG/train.txt",
        test="ELs_for_Rotate/Monarch_KG/test.txt",
        valid="ELs_for_Rotate/Monarch_KG/valid.txt"
    shell:
        """
        mkdir -p ELs_for_Rotate
        mkdir -p ELs_for_Rotate/Monarch_KG

        python Scripts/filter_edges_based_on_node_list.py -n {input.train_nodes} -e {input.train} -o {output.train}
        python Scripts/filter_edges_based_on_node_list.py -n {input.train_nodes} -e {input.valid} -o {output.valid}.tmp
        python Scripts/filter_edges_based_on_node_list.py -n {input.train_nodes} -e {input.test} -o {output.test}.tmp -s 17 -p 2 -b 18

        # remove semantic trips from valid that are in train
        python Scripts/remove_preexisting_triples.py -i {input.train} -n {output.valid}.tmp -o {output.valid}
        # remove triples from test that are in train or valid
        python Scripts/remove_preexisting_triples.py -i {input.train} -n {output.test}.tmp -o {output.test}.tmp2
        python Scripts/remove_preexisting_triples.py -i {output.valid} -n {output.test}.tmp2 -o {output.test}

        rm ELs_for_Rotate/Monarch_KG/*.tmp
        rm ELs_for_Rotate/Monarch_KG/*.tmp2
        """

rule map_HGNC_to_ENSG:
    input:
        "Resources/monarch-kg_nodes.sep_22.tsv",
    output:
        "Resources/HGNC_to_ENSG.tsv",
    run:
        hgnc_i = 0
        ensg_i = 7
        h2e = {}
        for line in open(input[0],'r'):
            row = line.strip().split('\t')
            ensgs = [x for x in row[ensg_i].split('|') if x.startswith('ENSEMBL')]
            if len(ensgs) == 0:
                continue
            h2e[row[hgnc_i]] = row[ensg_i][0]
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
        with open(output[0],'w') as out:
            for line in open(input.huri,'r'):
                row = line.strip().split('\t')
                row[0] = 'ENSEMBL:' + row[0]
                row[1] = 'ENSEMBL:' + row[1]
                if row[0].startswith('ENSEMBL') and row[1].startswith('ENSEMBL'):
                    if row[0] in hgnc_to_ensg:
                        hgnc1 = hgnc_to_ensg[row[0]]
                    if row[1] in hgnc_to_ensg:
                        hgnc2 = hgnc_to_ensg[row[1]]
                    out.write(f"{hgnc1}\tbiolink:interacts_with\t{hgnc2}\n")

rule make_monarch_HuRI:
    input:
        mapping="Resources/HGNC_to_ENSG.tsv",
        huri="Resources/huri.hgnc.triples.tsv",
        train="ELs_for_Rotate/Monarch_KG/train.txt",
        valid="ELs_for_Rotate/Monarch_KG/valid.txt",
        test="ELs_for_Rotate/Monarch_KG/test.txt"
    output:
        train="ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.train.tsv",
        valid="ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.valid.tsv",
        test= "ELs_for_Rotate/Monarch_HuRI/monarch_HuRI.test.tsv"
    shell:
        """
        mkdir -p ELs_for_Rotate/Monarch_HuRI
        python Scripts/create_monarch_huri.py -a {input.train} -b {input.huri} -p 'HGNC' -o {output.train}
        python Scripts/create_monarch_huri.py -a {input.valid} -b {input.huri} -p 'HGNC' -o {output.valid}
        python Scripts/create_monarch_huri.py -a {input.test} -b {input.huri} -p 'HGNC' -o {output.test}
        """

