process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.7--hf69ffd2_0':
        'biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0' }"

    publishDir "${params.outdir}/databases/amrfinderplus", mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    output:
    path "amrfinderdb.tar.gz", emit: db
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    amrfinder_update -d amrfinderdb
    tar czvf amrfinderdb.tar.gz -C amrfinderdb/\$(readlink amrfinderdb/latest) ./
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """

    stub:
    """
    touch amrfinderdb.tar
    gzip amrfinderdb.tar
    """
}
