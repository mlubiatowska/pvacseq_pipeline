#!/usr/bin/env nextflow

process ProximalVcf {
    cpus params.cpus
    tag "${name}"

    
    module params.bcftools
    module params.gatk
    module params.picard
    module params.mamba

    conda params.vepenv
    
    input:
    tuple val(name), path(phased_vcf), path(normal_vcf)

    output:
    tuple val(name), path("${name}_proximal_variants.vep.vcf.gz"), path("${name}_proximal_variants.vep.vcf.gz.tbi")

    script:
    """
    #Creating a phased VCF of proximal variants - load Picard GATK bgzip tabix
    
    #The sample names in the tumor.bam, the newly-created tumor_only.vcf, and the germline.vcf need to match. If they don’t you need to edit the sample names in the VCF files to match the tumor BAM file.
    bcftools reheader \
        -s <(printf "${name}-T") \
        -o ${name}.normal.renamed.vcf.gz \
        ${normal_vcf}
        

    #Combine somatic and germline variants using GATK's CombineVariants - equivalent of 
    echo "Merging VCFs into a single sample for proximate variant phasing..."
    gatk MergeVcfs \
        -I ${phased_vcf} \
        -I ${name}.normal.renamed.vcf.gz \
        -O combined_somatic_plus_germline.vcf \
        -R ${params.ref}

    #Sort combined VCF
    echo "Sorting and annotating merged VCF file..."
    gatk SortVcf \
        -I combined_somatic_plus_germline.vcf \
        -O combined_somatic_plus_germline.sorted.vcf

    #Not performing the phasing step using ReadBackedPhasing as we already have a phased VCF
    #Annotate VCF with VEP
    vep \
        --input_file combined_somatic_plus_germline.sorted.vcf --output_file ${name}_proximal_variants.vep.vcf \
        --format vcf --vcf --symbol --terms SO --tsl \
        --hgvs --fasta ${params.ref} \
        --offline --cache --dir_cache ${params.vepcache} \
        --plugin Downstream --plugin Wildtype \
        --dir_plugins ${params.program} \
        --fork ${params.cpus}
        #--pick

    #bgzip and index the phased VCF
    echo "zipping and indexing phased, merged VCF..."
    bgzip -c ${name}_proximal_variants.vep.vcf > ${name}_proximal_variants.vep.vcf.gz
    tabix -p vcf ${name}_proximal_variants.vep.vcf.gz

    """
    stub: 
    """
    touch ${name}_proximal_variants.vep.vcf.gz
    """
}
