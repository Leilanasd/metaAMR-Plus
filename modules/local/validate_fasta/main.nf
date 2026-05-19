process VALIDATE_FASTA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/validated_assemblies/${meta.id}", mode: params.publish_dir_mode, saveAs: { filename -> filename }
    conda "conda-forge::bash=5.2"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_validated.fasta"), emit: validated_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input_file = "${meta.id}_input.fasta"
    """
    # Validate that the input FASTA file exists
    if [ ! -f "${fasta}" ]; then
        echo "WARNING: Sample ${meta.id} has a missing or inaccessible FASTA assembly. PlasmidFinder and PlasClass will be skipped for this sample." >&2
        echo "Error: FASTA input file does not exist or is missing" >&2
        exit 1
    fi

    # Decompress if gzipped
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > "${input_file}"
    else
        cp "${fasta}" "${input_file}"
    fi

    # Sanitize FASTA headers to simple sequential contig IDs
    awk '/^>/ {print $1} !/^>/ {print}' "${input_file}" > "${meta.id}_validated.fasta"

    # Validate that the sanitized FASTA is non-empty and starts with ">"
    if [[ ! -s "${meta.id}_validated.fasta" ]] || [[ "\$(head -c 1 ${meta.id}_validated.fasta)" != ">" ]]; then
        echo "WARNING: Sample ${meta.id} has an invalid FASTA assembly. PlasmidFinder and PlasClass will be skipped for this sample." >&2
        echo "Error: Invalid FASTA file after validation" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | awk '{print \$4}')
    END_VERSIONS
    """
}
