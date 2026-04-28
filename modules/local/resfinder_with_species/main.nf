process RESFINDER_WITH_SPECIES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::resfinder=4.1.11"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.1.11--hdfd78af_0' :
        'quay.io/biocontainers/resfinder:4.1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(species_info), path(db)
    val target_species

    output:
    tuple val(meta), path("${meta.id}.resfinder_results.txt"),                                                          emit: amr_results
    tuple val("${task.process}"), val('resfinder'), eval("run_resfinder.py --version 2>&1 | sed 's/run_resfinder.py //'"), topic: versions, emit: versions_resfinder
    path "versions.yml",                                                                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_cmd = reads.name.toLowerCase().endsWith(".fq") || reads.name.toLowerCase().endsWith(".fastq") || reads.name.toLowerCase().endsWith(".fastq.gz")
        ? "-ifq ${reads}"
        : "-ifa ${reads}"
    """
    run_resfinder.py \\
        -acq \\
        ${input_cmd} \\
        -db_res ${db} \\
        -db_point ${db} \\
        -o . \\
        -l 0.6 \\
        -t 0.8

    resfinder_output=\$(ls *ResFinder_results_tab.txt 2>/dev/null)
    if [ -z "\$resfinder_output" ]; then
        echo "Error: ResFinder output file not found" >&2
        exit 1
    fi

    associate_amr_with_species.py \$resfinder_output ${species_info} "${target_species}" "${meta.id}" > ${meta.id}.resfinder_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$(run_resfinder.py --version 2>&1 | sed 's/run_resfinder.py //')
    END_VERSIONS
    """
}