process HAMRONIZATION_RESFINDER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity'
        ? 'https://depot.galaxyproject.org/singularity/hamronization:1.1.9--pyhdfd78af_0'
        : 'biocontainers/hamronization:1.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(report)
    val format

    output:
    tuple val(meta), path("*.json"), optional: true, emit: json
    tuple val(meta), path("*.tsv"),  optional: true, emit: tsv
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_resfinder"
    """
    if [ -s "${report}" ]; then
        hamronize \\
            resfinder \\
            ${report} \\
            ${args} \\
            --format ${format} \\
            > ${prefix}.${format}
    else
        echo "Input file ${report} is empty. Skipping hamronization."
        touch ${prefix}.${format}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(echo \$(hamronize --version 2>&1) | cut -f 2 -d ' ' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_resfinder"
    """
    touch ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(echo \$(hamronize --version 2>&1) | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
