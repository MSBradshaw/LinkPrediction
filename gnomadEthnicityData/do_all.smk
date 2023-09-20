
CHRs = [str(i) for i in range(1,23)] + ['X']

rule all:
    input:
        expand('work/INFO/chr{chr}_ancestry_groups.INFO',chr=CHRs),
        expand('work/PopSpecificAlleles/chr{chr}_ancestry_groups.tsv',chr=CHRs),
        'work/PopSpecificGenes/all_pop_specific_genes.tsv'

# rule download_gnomad:
#     input:
#         'do_all.smk' # an input file that should always exist (this file haha)
#     output:
#         gz='Data/gnomad.exomes.r2.1.1.sites.{chr}.vcf.gz',
#         tbi='Data/gnomad.exomes.r2.1.1.sites.{chr}.vcf.gz.tbi'
#     threads: 4
#     shell:
#         '''
#         wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{wildcards.chr}.vcf.bgz
#         wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{wildcards.chr}.vcf.bgz.tbi
#         mv gnomad.exomes.r2.1.1.sites.{wildcards.chr}.vcf.bgz {output.gz}
#         mv gnomad.exomes.r2.1.1.sites.{wildcards.chr}.vcf.bgz.tbi {output.tbi}
#         '''

rule get_AF_and_vep:
    input:
        'Data/gnomad.exomes.r2.1.1.sites.{chr}.vcf.gz'
    output:
        'work/INFO/chr{chr}_ancestry_groups.INFO'
    params:
        prefix=lambda wildcards, output: output[0].split('/')[-1].split('.')[0]
    threads: 4
    shell:
        '''
        mkdir -p work/INFO/

        vcftools --gzvcf {input} \
            --get-INFO AF_amr \
            --get-INFO AF_afr \
            --get-INFO AF_asj \
            --get-INFO AF_sas \
            --get-INFO AF_eas \
            --get-INFO AF_fin \
            --get-INFO AF_nfe_onf \
            --get-INFO AF_sas \
            --get-INFO vep \
            --out work/INFO/{params.prefix}
        '''

rule filter_AF_and_vep:
    input:
        'work/INFO/chr{chr}_ancestry_groups.INFO'
    output:
        'work/PopSpecificAlleles/chr{chr}_ancestry_groups.tsv'
    shell:
        '''
        mkdir -p work/PopSpecificAlleles/
        python Scripts/filter_AF_and_vep.py \
            --input {input} \
            --output {output}
        '''

rule get_pop_specific_genes:
    input:
        'work/PopSpecificAlleles/chr{chr}_ancestry_groups.tsv'
    output:
        'work/PopSpecificGenes/chr{chr}_pop_specific_genes.tsv'
    shell:
        '''
        mkdir -p work/PopSpecificGenes
        python Scripts/extract_pop_specific_genes.py --input {input} --output {output}
        '''

rule aggregate_pop_specific_genes:
    input:
        expand('work/PopSpecificGenes/chr{chr}_pop_specific_genes.tsv',chr=CHRs)
    output:
        'work/PopSpecificGenes/all_pop_specific_genes.tsv'
    shell:
        '''
        cat {input} > {output}
        python Scripts/divide_genes_by_pop.py --input {output} --output work/PopSpecificGenes/all_pop_specific_genes_
        '''