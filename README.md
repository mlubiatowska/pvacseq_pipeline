**This is Nextflow pipeline I run for pvacseq. This pipeline in run in Nextflow (v25).**

The pipeline loads picard/3.3.0-Java-17, GATK/4.4.0.0-GCCcore-12.3.0-Java-17, BCFtools/1.11 (specified in config file), and load conda env bio-perl-env.yml (saved in this directory.
 
It used 2 singularity images - one is /data/scratch/DMP/DUDMP/MYEGRP/ahyde01/lumos/singularity/whatshap_2.3.sif and another is /data/scratch/shared/SINGULARITY-DOWNLOAD/tools/.singularity/pvactools_latest.sif.

The singularity includes v6.0.5 of pvactools.

Check the config file to provide paths to programmes/packages/modules as necessary

**nextflow run https://github.com/mlubiatowska/pvacseq_pipeline.git -r main \
    --input ${INPUT_FILE} \
    --outdir ${OUTPUT_DIR} \
    --ref ${REF_DIR} \
    --expression ${EXPRESSION_MATRIX}**

It takes arguments:
--input is a comma-separeted table (csv) without a header, with sample name, path to phased normal vcf, somatic vcf (unphased), haplotagged tumour bam, normal bam, and path to HLA-LA output (R1 bestguesses.txt)
    format: sample_name,phased_normal_vcf,UNphased_tumour_vcf,haplotagged_tumour_bam,hla_la_output
--outdir is an output directory
--ref is path to reference genome .fa reference directory
--expression is Expression matrix - it is an optional input - expression matrix must be a tsv/txt file for the annotation to work.

**IMPORTANT NOTES:** 
1. *NAMING INSIDE FILES.* inside the expression matrix there must be --id-column ${sample_name}_B \ --expression-column "ENSEMBL_Gene_ID" and in the VCF --sample-name ${name}-T --> change in the config file if necessary
2. *INSTALLATION* all instalations specified by pavtools for pvacseq as well as toolf for input file preparation must be installed (follow https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep.html) 


