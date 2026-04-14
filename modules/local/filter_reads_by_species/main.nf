process FILTER_READS_BY_SPECIES {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(centrifuge_report), path(centrifuge_results)
    val target_species

    output:
    tuple val(meta), path("*.filtered_reads.txt"),                                                      emit: filtered_read_ids
    tuple val(meta), path("*.species_summary.txt"),                                                     emit: species_summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),         topic: versions, emit: versions_filter
    path "versions.yml",                                                                                emit: versions
    path "*.log",                                                                                       emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def species_arg = target_species ? "'${target_species}'" : ''
    """
    filter_reads_by_species.py \\
        ${centrifuge_report} \\
        ${centrifuge_results} \\
        ${prefix}.filtered_reads.txt \\
        ${prefix}.species_summary.txt \\
        ${species_arg} 2> ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}