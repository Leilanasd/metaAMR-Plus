process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.12.8--h283d18e_0':
        'biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}_amrfinder.tsv")            , emit: report
    tuple val(meta), path("${prefix}_amrfinder_mutations.tsv"), emit: mutation_report, optional: true
    tuple val(meta), val("amrfinderplus")                       , emit: format
    path "versions.yml"                             , emit: versions
    env VER                                         , emit: tool_version
    env DBVER                                       , emit: db_version

    when:
    task.ext.when == null || task.ext.when

    script:
def args = task.ext.args ?: ''
def is_compressed_fasta = fasta.getName().endsWith(".gz") ? true : false
def is_compressed_db = db.getName().endsWith(".gz") ? true : false
prefix = task.ext.prefix ?: "${meta.id}"
organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
fasta_name = fasta.getName().replace(".gz", "")
fasta_param = "-n"
if (meta.containsKey("is_proteins")) {
    if (meta.is_proteins) {
        fasta_param = "-p"
    }
}
"""
set -e  # Exit on any error
if [ "$is_compressed_fasta" == "true" ]; then
    gzip -c -d $fasta > $fasta_name
fi

DB_PATH=""

if [ "$is_compressed_db" == "true" ]; then
    mkdir amrfinderdb
    tar xzvf $db -C amrfinderdb
    DB_PATH="amrfinderdb"
else
    DB_PATH="$db"
fi

echo "Using AMRFinder database at: \$DB_PATH"
ls -lah "\$DB_PATH" || true

amrfinder \
    $fasta_param $fasta_name \
    $organism_param \
    $args \
    --database "\$DB_PATH" \
    --threads $task.cpus > ${prefix}_amrfinder.tsv

VER=\$(amrfinder --version)
DBVER=\$(echo \$(amrfinder --database "\$DB_PATH" --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    amrfinderplus: \$(amrfinder --version)
    amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database "\$DB_PATH" --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
END_VERSIONS
"""

}