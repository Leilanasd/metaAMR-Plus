process GENERATE_REPORT {
    tag "all_samples"
    label 'process_single'

    conda "conda-forge::python=3.8"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.8--1' :
        'quay.io/biocontainers/python:3.8' }"

    publishDir "${params.outdir}/report", mode: params.publish_dir_mode, overwrite: true

    input:
    val  results_dir
    val  run_name

    output:
    path "metaamr_report.html", emit: report
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def run_label = run_name ? "--run_name \"${run_name}\"" : ""
    """
    python ${projectDir}/bin/generate_report.py \\
        --results_dir ${results_dir} \\
        --outdir      . \\
        ${run_label}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        generate_report: 1.0.0
    END_VERSIONS
    """
}
