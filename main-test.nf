#!/usr/bin/env nextflow
nextflow.enable.dsl     = 2
/*
* Pipeline parameters
*/

include { TumourVcfEdits } from './modules/long_read_phasing.nf'
include { LongReadPhasing } from './modules/long_read_phasing.nf'
include { LongReadPreprocessing } from './modules/long_read_preprocessing.nf'
include { ProximalVcf } from './modules/proximal_vcf.nf'
include { PvacseqLongread } from './modules/pvacseq_longread.nf'

workflow {
    main:
    // Parse input CSV and create two separate channels 
    if (params.pvac_only) {
    // Pvac-only mode: parse pre-computed files directly
    pvac_data = Channel
        .fromPath(params.pvac_only)
        .splitCsv(header: ['name', 'annotated_tumour_vcf', 'annotated_tumour_vcf_tbi', 'proximal_vcf', 'proximal_vcf_tbi', 'hla_alleles'], sep: ",")
        .map{ row -> tuple(
            row.name,
            file(row.annotated_tumour_vcf, checkIfExists: true),
            file(row.annotated_tumour_vcf_tbi, checkIfExists: true),
            file(row.proximal_vcf, checkIfExists: true),
            file(row.proximal_vcf_tbi, checkIfExists: true),
            file(row.hla_alleles, checkIfExists: true)
        ) }
        
    annotated_prox_alleles_ch = pvac_data
        
    } else {

    input_data = Channel
        .fromPath(params.input)
        .splitCsv(header: ['name', 'normal_vcf', 'tumour_vcf', 'hg_bam', 'hla_alleles'], sep: ",")
        .map{ row -> tuple(
            row.name, 
            file(row.normal_vcf, checkIfExists: true),      // Convert to Path object
            file(row.tumour_vcf, checkIfExists: true),      // Convert to Path object
            file(row.hg_bam, checkIfExists: true),          // Convert to Path object
            file(row.hla_alleles, checkIfExists: true)      // Convert to Path object
        ) }
        //.map{ row -> tuple(row.name, row.normal_vcf, row.tumour_vcf, row.hg_bam, row.hla_alleles) }


    // Channel 1: name, normal_vcf
    normal_channel = input_data
        .map{ name, normal_vcf, tumour_vcf, hg_bam, hla_alleles -> tuple(name, normal_vcf) }


    // Channel 2: name, tumour_vcf, hg_bam
    tumour_vcf_channel = input_data
        .map{ name, normal_vcf, tumour_vcf, hg_bam, hla_alleles -> tuple(name, tumour_vcf)}


    tumour_bam_channel = input_data
        .map{ name, normal_vcf, tumour_vcf, hg_bam, hla_alleles -> tuple(name, hg_bam, "${hg_bam}.bai")}
   

    //Channel 3: name, hla_alleles
    alleles_channel = input_data
        .map{ name, normal_vcf, tumour_vcf, hg_bam, hla_alleles -> tuple(name, hla_alleles)}
    

    //Workflow logic
    TumourVcfEdits(tumour_vcf_channel)
    tumour_channel = TumourVcfEdits.out.join(tumour_bam_channel)
    //tumour_channel.view()

    LongReadPhasing(tumour_channel)
    phased_ch = LongReadPhasing.out

    joined_ch = phased_ch.join(normal_channel)
    ProximalVcf(joined_ch)
    proximal_ch = ProximalVcf.out

    LongReadPreprocessing(phased_ch)
    annotated_ch = LongReadPreprocessing.out

    annotated_prox_ch = annotated_ch.join(proximal_ch)

    annotated_prox_alleles_ch = annotated_prox_ch.join(alleles_channel)

    }

    PvacseqLongread(annotated_prox_alleles_ch)

    publish:
    pvacseq_neoag         = PvacseqLongread.out
}

output {    
    pvacseq_neoag          { path { name, pvacseq_neoag -> "pvac_out/${name}" } }
}
