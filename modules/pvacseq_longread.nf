#!/usr/bin/env nextflow

process PvacseqLongread {
    cpus params.cpus
    tag "${name}" 

    //scratch '/data/scratch/DGE/DUDGE/MOPOPGEN/mlubiatowska/tmp'
    //'$SCRATCH_TMPDIR'

    input:
    tuple val (name), path(annotated_tumour_vcf), path(annotated_tumour_vcf_tbi), path(proximal_vcf), path(proximal_vcf_tbi), path(hla_alleles)

    output:
    tuple val(name), path("${name}_neoag")

    script:
    """
    #extracting the column alleles HLA from the input tsv file and adding HLA before the allele names to fit pvacseq format
    awk -F'\\t' '
        NR==1 {
            for (i=1; i<=NF; i++) if (\$i=="Allele") col=i
            if (!col) {
                print "ERROR: Allele column not found" > "/dev/stderr"
                exit 1
            }
            next
        }

        \$col != "" {
            allele = \$col

            # remove whitespace
            gsub(/[[:space:]]+/, "", allele)

            # remove HLA- prefix if present
            sub(/^HLA-/, "", allele)

            # keep only 2-field resolution
            if (match(allele, /^([^*]+)\\*([0-9]+:[0-9]+)/, m)) {
                gene = m[1]
                twofield = gene "*" m[2]

                # Class I → add HLA-
                if (gene ~ /^(A|B|C|E|F|G)\$/) {
                    print "HLA-" twofield
                }
                # Class II → no prefix
                else {
                    print twofield
            }
        }
    }
    ' '${hla_alleles}' | sort -u | paste -sd "," - > "${name}_HLA_alleles.txt"


    #running singularity in an example test set, where out6 is an emplty output file and pvacseq_example_data includes exmaple dataset provided by the pVACtools 
    
    HLA_ALLELES=\$(cat ${name}_HLA_alleles.txt)
            
    pvacseq run \
        ${annotated_tumour_vcf} \
        ${name}-T \
        \${HLA_ALLELES} \
        all \
        ${name}_neoag \
        -t ${params.cpus} \
        --phased-proximal-variants-vcf ${proximal_vcf} \
        --iedb-install-directory /opt/iedb


    """

    stub: 
    """
    mkdir ${name}
    touch ${name}/neoantigens.tsv
    """
}
