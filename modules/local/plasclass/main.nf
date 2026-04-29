process PLASCLASS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plasclass=0.1.1"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/plasclass:0.1.1--pyhdfd78af_0' :
        'quay.io/biocontainers/plasclass:0.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    val threshold

    output:
    tuple val(meta), path("*.plasclass.txt"),     emit: report
    tuple val(meta), path("*.plasclass_classified.txt"), emit: classified
    tuple val("${task.process}"), val('plasclass'), eval("classify_fasta.py --version 2>&1 | grep -oP 'PlasClass \\K[0-9]+\\.[0-9]+\\.[0-9]+' || echo '0.1.1'"), topic: versions, emit: versions_plasclass
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plasclass_threshold = threshold ?: 0.9
    """
    if [[ ${fasta} == *.gz ]]; then
        gunzip -c ${fasta} > uncompressed.fasta
        input_fasta=uncompressed.fasta
    else
        input_fasta=${fasta}
    fi

    classify_fasta.py \\
        $args \\
        -f \$input_fasta \\
        -o ${prefix}.plasclass.txt

    awk -v threshold=${plasclass_threshold} '
    BEGIN { print "Contig_ID\\tClassification" }
    { if (\$2 >= threshold) print \$1"\\tplasmid";
      else print \$1"\\tchromosome" }
    ' ${prefix}.plasclass.txt > ${prefix}.plasclass_classified.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasclass: \$(classify_fasta.py --version 2>&1 | grep -oP 'PlasClass \\K[0-9]+\\.[0-9]+\\.[0-9]+' || echo "0.1.1")
    END_VERSIONS
    """
}