process PvacseqLongread {
    cpus params.cpus
    tag "${name}" 
    container '/data/scratch/DGE/DUDGE/MOPOPGEN/mlubiatowska/pvacseq/pvactools/pvactools.sif'

    input:
    tuple val(name), path(annotated_tumour_vcf), path(annotated_tumour_vcf_tbi), path(proximal_vcf), path(proximal_vcf_tbi), path(hla_alleles)

    output:
    tuple val(name), path("${name}_neoag")

    //containerOptions "--bind /mnt:/mnt --bind /data:/data --bind ${task.workDir}/tmp:/tmp"

    script:
    """
    # Nextflow working directory is already mounted in container
    # CREATE SUBDIRECTORY - Keeps work dir clean
    mkdir -p ./tmp
    
    # sUPDATE TMPDIR - Point to subdirectory
    export TMPDIR=\${PWD}/tmp
    export TEMP=\${PWD}/tmp
    export TMP=\${PWD}/tmp
    export SINGULARITY_TMPDIR=\${PWD}/tmp  # Also add this
    export APPTAINERENV_TMPDIR=\${PWD}/tmp
    export APPTAINER_TMPDIR=\${PWD}/tmp     # And this

    # Extract HLA alleles and format for pVACseq
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
