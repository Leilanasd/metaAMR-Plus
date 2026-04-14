process EXTRACT_FILTERED_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqtk=1.3"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(read_ids), path(original_reads)

    output:
    tuple val(meta), path("*.filtered.fastq.gz"),                                                          emit: filtered_reads
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | grep Version | sed 's/.*Version: //'"), topic: versions, emit: versions_seqtk
    path "versions.yml",                                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk subseq ${original_reads} <(cut -f1 ${read_ids}) | gzip > ${prefix}.filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}