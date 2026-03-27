#!/usr/bin/env nextflow

process LongReadPreprocessing {
    cpus params.cpus
    tag "${name}"

    conda params.vepenv

    /*
    module params.picard
    module params.gatk
    module params.bcftools
    */

    input:
    tuple val(name), path(phased_vcf) 

    output:
    tuple val (name), path("${name}.annotated.tumour.vcf.gz"), path("${name}.annotated.tumour.vcf.gz.tbi")

    script:

    """
    #Indexing the phased VCF file 
    echo "Indexing the phased VCF..."
    tabix -p vcf ${phased_vcf}

    #Not adding genotype sample information to the .vcf using vatools (included in my mamba) - already added and phased using whatshap haplotagphase based on haplotagged bam files from lumos

    #Annotating VCF with VEP
    echo "Running VEP annotation on the merged file..."
    vep \
        --input_file ${phased_vcf} \
        --output_file ${name}.vep.vcf \
        --format vcf --vcf --symbol --terms SO --mane_select --canonical --tsl --biotype \
        --hgvs --fasta ${params.ref} \
        --offline --cache --dir_cache ${params.vepcache} \
        --plugin Frameshift --plugin Wildtype \
        --dir_plugins ${params.program} \
        --fork ${params.cpus} \
        --transcript_version
        #--pick 
        
    echo "VEP annotation completed."

    # Running expression annotation if expression matrix is provided, otherwise copying the VEP-annotated VCF to the next step
    ${params.expression ? """
    echo "Running expression annotation on the VEP-annotated VCF..."
    vcf-expression-annotator \
        ${name}.vep.vcf \
        ${params.expression} \
        custom \
        gene \
        --id-column ${params.expression_geneid_column} \
        --expression-column ${name} \
        --sample-name ${name}-T \
        > ${name}.vep.gx.vcf
    """ : """
    echo "Skipping expression annotation (--expression not provided)"
    cp ${name}.vep.vcf ${name}.vep.gx.vcf
    """}
    
    #bgzip and index the pVACseq main input VCF
    echo "zipping input VCF..."
    bgzip -c ${name}.vep.gx.vcf > ${name}.annotated.tumour.vcf.gz
    tabix -p vcf ${name}.annotated.tumour.vcf.gz

    #Not adding coverage data to the vcf files from short read sequencing - there is already info about adding coverage data - contains AD, DP (and AF)

    """
    stub: 
    """
    touch ${name}.annotated.tumour.vcf.gz
    touch ${name}.annotated.tumour.vcf.gz.tbi
    """
}