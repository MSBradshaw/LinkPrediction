# --- Parameter sweep options for node2vec ---
N2V_DIMs = [32,64,128,256]
N2V_WALK_LENGTHS = [10,40,80,160]
N2V_NUM_WALKS = [5,10,15,20]
N2V_CONTEXT_SIZE = [5,10,15]
# --------------------------------------------


YEARS = [2019,2020,2021,2022]
YEARS_short = YEARS[:-1]
short_YEARS_short = YEARS_short[1:]

assert len(YEARS) == len(YEARS_short) + 1

rule all:
    input:
        # Node2Vec embedding output
        # expand("Embeddings/String_HPO_2019.all_hpo.numbered.n2v.{dim}.{length}.{walks}.{k}.e1.p1.q1.emb",
        #     dim=N2V_DIMs,
        #     length=N2V_WALK_LENGTHS,
        #     walks=N2V_NUM_WALKS,
        #     k=N2V_CONTEXT_SIZE),
        # New G2P edge lists
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.numbered.edgelist.txt', year=short_YEARS_short),
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.edgelist.txt', year=short_YEARS_short),
        # set up edgelists for rotatE
        expand('ELs_for_Rotate/String_HPO_{year}.all_hpo/entities.dict',year=YEARS[:-2]),
        expand('ELs_for_Rotate/String_HPO_{year}.all_hpo/relations.dict',year=YEARS),
        expand('ELs_for_Rotate/String_HPO_{year}.all_hpo/train.txt',year=YEARS),
        # validation and test sets for 2019 and 2022 for rotatE
        'ELs_for_Rotate/String_HPO_2019.all_hpo/valid.txt',
        'ELs_for_Rotate/String_HPO_2019.all_hpo/test.txt',
        'ELs_for_Rotate/String_HPO_2020.all_hpo/valid.txt',
        'ELs_for_Rotate/String_HPO_2020.all_hpo/test.txt',

# generate embeddings with node2vec with a variety of parameter options
rule node2vec:
    input:
        el = 'Edgelists/String_HPO_2019.all_hpo.numbered.edgelist.txt'
    output:
        emb = 'Embeddings/String_HPO_2019.all_hpo.numbered.n2v.{dim}.{length}.{walks}.{k}.e1.p1.q1.emb'
    log:
        log = 'Logs/String_HPO_2019.all_hpo.numbered.n2v.{dim}.{length}.{walks}.{k}.e1.p1.q1.emb.log'
    
    # this should be parameterized to be the max number given to snakemake, node2vec automatically uses all available cores - you cannot change this.
    threads: 64
    
    shell:
        """
        mkdir -p Embeddings
        time node2vec -i:{input} \
        -o:{output} \
        -d:{wildcards.dim} \
        -l:{wildcards.length} \
        -r:{wildcards.walks} \
        -k:{wildcards.k} \
        -e:1 \
        -p:1 \
        -q:1 \
        -v &> {log}
        """

# rule takes each year (t) and find the g2p edges that were added in the next year (t+1) and saved it in a file labeled with t+1
rule generate_new_g2p_lists:
    input:
        els = expand('Edgelists/String_HPO_{year}.all_hpo.edgelist.txt', year=YEARS_short),
        maps = expand('Edgelists/String_HPO_{year}.all_hpo.nodenames.txt', year=YEARS_short)
    output:
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.numbered.edgelist.txt', year=short_YEARS_short),
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.edgelist.txt', year=short_YEARS_short)
    params:
        first_year = YEARS_short[0],
        last_year = YEARS_short[-1]
    shell:
        """
        mkdir -p NewG2PEdges
        # for each year first_year to last_year, print year
        for i in {{{params.first_year}..{params.last_year}}}
        do
            echo "$i"
            # add 1 to i  
            j=$((i+1))
            python Scripts/create_g2p_lists.py \
                --edgelist1 Edgelists/String_HPO_${{i}}.all_hpo.edgelist.txt \
                --edgelist2 Edgelists/String_HPO_${{j}}.all_hpo.edgelist.txt \
                --node_map1 Edgelists/String_HPO_${{i}}.all_hpo.nodenames.txt \
                --output NewG2PEdges/String_HPO_${{j}}.all_hpo.edgelist.txt \
                --output_numbered NewG2PEdges/String_HPO_${{j}}.all_hpo.numbered.edgelist.txt
        done 
        """

# rule that takes the string hpo edgelist and creates a directory with files rotate needs for embeddings
rule string_hpo_prep_for_rotatE:
    input:
        el='Edgelists/String_HPO_{year}.all_hpo.edgelist.txt',
        numbered_el='Edgelists/String_HPO_{year}.all_hpo.numbered.edgelist.txt'
    output:
        rels='ELs_for_Rotate/String_HPO_{year}.all_hpo/relations.dict', # number 2 name for edges
        # test='ELs_for_Rotate/String_HPO_{year}.all_hpo/test.txt', # triples for testing
        train='ELs_for_Rotate/String_HPO_{year}.all_hpo/train.txt', # triples for training
        # valid='ELs_for_Rotate/String_HPO_{year}.all_hpo/valid.txt', # triples for validation
    shell:
        """
        mkdir -p ELs_for_Rotate
        mkdir -p ELs_for_Rotate/String_HPO_{wildcards.year}.all_hpo/

        # create relations.dict
        echo "0\tSTRING2HPO" > {output.rels}
        echo "1\tSTRING2STRING" >> {output.rels}
        echo "2\tHPO2HPO" >> {output.rels}

        # train should just be the current year, but you need to add relationships as a middle column
        python Scripts/add_relationship_to_el.py --input {input.el} --output {output.train}
        """

rule add_validation_test_to_string_hpo_for_rotate:
    input:
        expand('ELs_for_Rotate/String_HPO_{year}.all_hpo/relations.dict',year=YEARS),
        expand('ELs_for_Rotate/String_HPO_{year}.all_hpo/train.txt',year=YEARS),
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.numbered.edgelist.txt', year=short_YEARS_short),
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.edgelist.txt', year=short_YEARS_short)
    output:
        valid19='ELs_for_Rotate/String_HPO_2019.all_hpo/valid.txt',
        test19='ELs_for_Rotate/String_HPO_2019.all_hpo/test.txt',
        valid20='ELs_for_Rotate/String_HPO_2020.all_hpo/valid.txt',
        test20='ELs_for_Rotate/String_HPO_2020.all_hpo/test.txt'

    shell:
        """
        # 2020 is validation for 2019
        python Scripts/add_relationship_to_el.py --input NewG2PEdges/String_HPO_2020.all_hpo.edgelist.txt --output {output.valid19}
        # 2021 is test for 2019
        python Scripts/add_relationship_to_el.py --input NewG2PEdges/String_HPO_2021.all_hpo.edgelist.txt --output {output.test19}

        # 2021 is validation for 2020
        python Scripts/add_relationship_to_el.py --input NewG2PEdges/String_HPO_2021.all_hpo.edgelist.txt --output {output.valid20}
        # 2022 is test for 2020
        python Scripts/add_relationship_to_el.py --input NewG2PEdges/String_HPO_2022.all_hpo.edgelist.txt --output {output.test20}
        """

rule add_entities:
    input:
        rels='ELs_for_Rotate/String_HPO_{year}.all_hpo/relations.dict', # number 2 name for edges
        test='ELs_for_Rotate/String_HPO_{year}.all_hpo/test.txt', # triples for testing
        train='ELs_for_Rotate/String_HPO_{year}.all_hpo/train.txt', # triples for training
        valid='ELs_for_Rotate/String_HPO_{year}.all_hpo/valid.txt', # triples for validation
        node_names='Edgelists/String_HPO_{year}.all_hpo.nodenames.txt'
    output:
        ents='ELs_for_Rotate/String_HPO_{year}.all_hpo/entities.dict', # number 2 name for nodes
    shell:
        """
        # create entities.dict
        python Scripts/create_entities_file.py --train {input.train} --validate {input.valid} --test {input.test} --mapping {input.node_names} --output {output.ents}
        """
