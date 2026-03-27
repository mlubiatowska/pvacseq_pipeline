The input is a csv (comma-seperated) file, including:
sample_name,phased_normal_vcf,UNphased_tumour_vcf,haplotagged_tumour_bam,hla_la_output

Expression matrix is an optional input - expression matrix must be a tsv/txt file for the annotation to work.
NOTE:inside the expression matrix there must be --id-column ${sample_name}_B \ --expression-column "ENSEMBL_Gene_ID" and in the VCF --sample-name ${name}-T --> change in the config file if necessary 

check the config file to provide paths to programmes/packages/modules as necessary

reference directory must inclue all references for all softwares, as recommended by pvactools input file preperation at https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep.html


The command for running the pipeline looks like this:
nextflow run ${WORKFLOW_DIR}/main.nf \
    --input ${INPUT_FILE} \
    --outdir ${OUTPUT_DIR} \
    --ref ${REF_DIR} \
    --expression ${EXPRESSION_MATRIX}
