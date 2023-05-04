# --- Parameter sweep options for node2vec ---
N2V_DIMs = [32,64,128,256]
N2V_WALK_LENGTHS = [10,40,80,160]
N2V_NUM_WALKS = [5,10,15,20]
N2V_CONTEXT_SIZE = [5,10,15]
# --------------------------------------------


YEARS = [2019,2020,2021,2022]
YEARS_short = YEARS[:-1]
YEARS_short_short = YEARS_short[:-1]

assert len(YEARS) == len(YEARS_short) + 1

rule all:
    input:
        # Node2Vec embedding output
        expand("Embeddings/String_HPO_2019.all_hpo.numbered.n2v.{dim}.{length}.{walks}.{k}.e1.p1.q1.emb",
            dim=N2V_DIMs,
            length=N2V_WALK_LENGTHS,
            walks=N2V_NUM_WALKS,
            k=N2V_CONTEXT_SIZE)
        # New G2P edge lists
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.numbered.edgelist.txt', year=YEARS_short_short),
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.edgelist.txt', year=YEARS_short_short)

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
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.numbered.edgelist.txt', year=YEARS_short_short),
        expand('NewG2PEdges/String_HPO_{year}.all_hpo.edgelist.txt', year=YEARS_short_short)
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
                --edgelist1 Edgelists/String_HPO_${{$i}.all_hpo.numbered.edgelist.txt \
                --edgelist2 Edgelists/String_HPO_${{$j}.all_hpo.numbered.edgelist.txt \
                --node_map1 Edgelists/String_HPO_${{$i}.all_hpo.nodenames.txt \
                --output NewG2PEdges/String_HPO_${{$j}.all_hpo.edgelist.txt \
                --output_numbered NewG2PEdges/String_HPO_${{$j}.all_hpo.numbered.edgelist.txt
        done 
        """

# conda activate link
# rj -n 'LP-smk-embeddings' -t 72:00:00 -T 64 -m 500G -p sandbox  -c 'snakemake -s make_embeddings.smk --cores 64'
