
CHRs = [str(i) for i in range(1,23)] + ['X']

rule all:
    input:
        expand('work/INFO/chr{chr}_ancestry_groups.INFO',chr=CHRs),
        expand('work/PopSpecificAlleles/chr{chr}_ancestry_groups.tsv',chr=CHRs),
        'work/PopSpecificGenes/all_pop_specific_genes.tsv',
        'work/PopMajoritySpecificGenes/all_pop_specific_genes.tsv',
        'work/PopPluralitySpecificGenes/all_pop_specific_genes.tsv',
        expand('work/PopSpecificPathogenicAllelesHPOTerms/chr{chr}_pathogenic_ancestry_groups.tsv',chr=CHRs),
        'work/PopSpecificPathogenicAllelesHPOTerms/all_pathogenic_ancestry_groups.tsv',
        'work/PopSpecificPathogenicAllelesMONDOTerms/all_pathogenic_ancestry_groups.tsv'

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
        af='work/PopSpecificAlleles/chr{chr}_ancestry_groups.tsv'
    shell:
        '''
        mkdir -p work/PopSpecificAlleles/
        python Scripts/filter_AF_and_vep.py \
            --input {input} \
            --output {output.af}
        '''

rule filter_AF_pathogenic:
    input:
        'work/INFO/chr{chr}_ancestry_groups.INFO'
    output:
        patho='work/PopSpecificPathogenicAlleles/chr{chr}_pathogenic_ancestry_groups.tsv'
    shell:
        '''
        mkdir -p work/PopSpecificPathogenicAlleles

        python Scripts/filter_AF_and_vep.py \
            --input {input} \
            --output {output.patho} -m p
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

rule get_majority_and_plurality_pop_specific_genes:
    input:
        'work/PopSpecificAlleles/chr{chr}_ancestry_groups.tsv'
    output:
        major='work/PopMajoritySpecificGenes/chr{chr}_pop_specific_genes.tsv',
        plural='work/PopPluralitySpecificGenes/chr{chr}_pop_specific_genes.tsv'
    shell:
        '''
        mkdir -p work/PopMajoritySpecificGenes
        mkdir -p work/PopPluralitySpecificGenes
        python Scripts/extract_pop_specific_genes.py --input {input} --output {output.major} -m m
        python Scripts/extract_pop_specific_genes.py --input {input} --output {output.plural} -m p
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

rule aggregate_majority_and_plurality_pop_specific_genes:
    input:
        major=expand('work/PopMajoritySpecificGenes/chr{chr}_pop_specific_genes.tsv',chr=CHRs),
        plural=expand('work/PopPluralitySpecificGenes/chr{chr}_pop_specific_genes.tsv',chr=CHRs)
    output:
        major='work/PopMajoritySpecificGenes/all_pop_specific_genes.tsv',
        plural='work/PopPluralitySpecificGenes/all_pop_specific_genes.tsv'
    shell:
        '''
        cat {input.major} > {output.major}
        cat {input.plural} > {output.plural}
        python Scripts/divide_genes_by_pop.py --input {output.major} --output work/PopMajoritySpecificGenes/all_majority_pop_specific_genes_
        python Scripts/divide_genes_by_pop.py --input {output.plural} --output work/PopPluralitySpecificGenes/all_plurality_pop_specific_genes_
        '''

rule get_clinvar_file:
    input:
        'do_all.smk' # an input file that should always exist (this file haha)
    output:
        'work/Data/clinvar.vcf.gz',
        'work/Data/clinvar.vcf.gz.tbi'
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
        wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
        mv clinvar.vcf.gz work/Data/clinvar.vcf.gz
        mv clinvar.vcf.gz.tbi work/Data/clinvar.vcf.gz.tbi
        """

rule patho_pop_specific_alleles_to_hpo:
    input:
        asv='work/PopSpecificPathogenicAlleles/chr{chr}_pathogenic_ancestry_groups.tsv',
        clinvar='work/Data/clinvar.vcf.gz',
        clinvar_tbi='work/Data/clinvar.vcf.gz.tbi' 
    output:
        'work/PopSpecificPathogenicAllelesHPOTerms/chr{chr}_pathogenic_ancestry_groups.tsv'
    shell:
        '''
        mkdir -p work/PopSpecificPathogenicAllelesHPOTerms
        python Scripts/patho_asv_to_hpo.py --input {input.asv} --output {output} --clinvar {input.clinvar}
        '''

rule aggregate_patho_pop_specific_alleles_to_hpo:
    input:
        expand('work/PopSpecificPathogenicAllelesHPOTerms/chr{chr}_pathogenic_ancestry_groups.tsv',chr=CHRs)
    output:
        'work/PopSpecificPathogenicAllelesHPOTerms/all_pathogenic_ancestry_groups.tsv'
    shell:
        '''
        cat {input} > {output}
        python Scripts/divide_avs_hpos_by_pop.py --input {output} --output work/PopSpecificPathogenicAllelesHPOTerms/all_pathogenic_ancestry_groups_
        '''

### ------------------- MONDO ------------------- ###

rule patho_pop_specific_alleles_to_mondo:
    input:
        asv='work/PopSpecificPathogenicAlleles/chr{chr}_pathogenic_ancestry_groups.tsv',
        clinvar='work/Data/clinvar.vcf.gz',
        clinvar_tbi='work/Data/clinvar.vcf.gz.tbi' 
    output:
        'work/PopSpecificPathogenicAllelesMONDOTerms/chr{chr}_pathogenic_ancestry_groups.tsv'
    shell:
        '''
        mkdir -p work/PopSpecificPathogenicAllelesMONDOTerms
        python Scripts/patho_asv_to_mondo.py --input {input.asv} --output {output} --clinvar {input.clinvar}
        '''

rule aggregate_patho_pop_specific_alleles_to_hpo:
    input:
        expand('work/PopSpecificPathogenicAPopSpecificPathogenicAllelesMONDOTermsllelesHPOTerms/chr{chr}_pathogenic_ancestry_groups.tsv',chr=CHRs)
    output:
        'work/PopSpecificPathogenicAllelesMONDOTerms/all_pathogenic_ancestry_groups.tsv'
    shell:
        '''
        cat {input} > {output}
        python Scripts/divide_avs_hpos_by_pop.py --input {output} --output work/PopSpecificPathogenicAllelesMONDOTerms/all_pathogenic_ancestry_groups_
        '''