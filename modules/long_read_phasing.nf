#!/usr/bin/env nextflow

process TumourVcfEdits {
    cpus 2
    tag "${name}"

    module params.bcftools

    input:
    tuple val(name), path(tumour_vcf)

    output:
    tuple val(name), path("regenotyped.vcf.gz"), path("regenotyped.vcf.gz.tbi")

    script:
    """
    #runnning whatshap for tumour VCF - note whatshap allows merging of all samples within a vcf 
    bcftools view -f PASS ${tumour_vcf} -Oz -o PASS.tumour.vcf.gz

    #remove GT field and re-assinginng the right genotype to fitlered (pass-only) vcf file
    bcftools annotate -x FORMAT/GT PASS.tumour.vcf.gz -Oz -o PASS.noGT.tumour.vcf.gz

    vcf-genotype-annotator PASS.noGT.tumour.vcf.gz ${name}-T 0/1 -o regenotyped.vcf
    bgzip -c regenotyped.vcf > regenotyped.vcf.gz
    tabix -p vcf regenotyped.vcf.gz
    """
    
    stub: 
    """
    touch regenotyped.vcf.gz
    touch regenotyped.vcf.gz.tbi
    """
}

process LongReadPhasing {
    cpus params.cpus
    tag "${name}"

    input:
    tuple val(name), path(tumour_vcf), path(tumour_vcf_tbi), path(hg_bam), path(hg_bam_bai)

    output:
    tuple val(name), path("${name}.phased.tumour.PASS.vcf.gz")

    script:
    """
    # Run whatshap on the minimal test VCF - singularity path has been added to this script 
    whatshap haplotagphase \
        -r ${params.ref} \
        -o ${name}.phased.tumour.PASS.vcf.gz \
        --ignore-read-groups \
        regenotyped.vcf.gz \
        ${hg_bam}

    """

    stub: 
    """
    touch ${name}.phased.tumour.PASS.vcf.gz
    """
}