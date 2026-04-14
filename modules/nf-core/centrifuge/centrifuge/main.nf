
process CENTRIFUGE_CENTRIFUGE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4.2--hdcf5f25_0' :
        'biocontainers/centrifuge:1.0.4.2--hdcf5f25_0' }"

    input:
    tuple val(meta), path(reads)
    path db
    val save_unaligned
    val save_aligned

    output:
    tuple val(meta), path('*_centrifuge_report.txt'),              emit: report
    tuple val(meta), path('*_centrifuge_results.txt'),             emit: results
    tuple val(meta), path('*.{sam,tab}'),          optional: true, emit: sam
    tuple val(meta), path('*.mapped.fastq{,.1,.2}.gz'),   optional: true, emit: fastq_mapped
    tuple val(meta), path('*.unmapped.fastq{,.1,.2}.gz'), optional: true, emit: fastq_unmapped
    tuple val("${task.process}"), val('centrifuge'), eval("centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //'"), topic: versions, emit: versions_centrifuge
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_read = reads instanceof List ? reads[0] : reads
    def is_fasta = single_read.name =~ /\.(fasta|fa|fna)(\.gz)?$/
    def is_gzipped = single_read.name.endsWith('.gz')
    def input_type = is_fasta ? "-f" : (meta.single_end ? "-U" : "-1")
    def unaligned = ''
    def aligned = ''
    if (!is_fasta) {
        if (meta.single_end) {
            unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
            aligned   = save_aligned   ? "--al-gz ${prefix}.mapped.fastq.gz"   : ''
        } else {
            unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
            aligned   = save_aligned   ? "--al-conc-gz ${prefix}.mapped.fastq.gz"   : ''
        }
    }
    // For FASTA: pipe through zcat if gzipped, pass via stdin
    // For FASTQ: pass directly as file arguments
    def input_command = ''
    def reads_arg = ''
    if (is_fasta) {
        input_command = is_gzipped ? "zcat ${single_read} |" : "cat ${single_read} |"
        reads_arg = "-f -"
    } else if (meta.single_end) {
        input_command = ''
        reads_arg = "-U ${reads}"
    } else {
        input_command = ''
        reads_arg = "-1 ${reads[0]} -2 ${reads[1]}"
    }
    """
    db_name=`find -L ${db} -name "*.1.cf" -not -name "._*" | sed 's/\\.1.cf\$//'`

    mkdir ./temp

    ${input_command} centrifuge \\
        -x \$db_name \\
        --temp-directory ./temp \\
        -p ${task.cpus} \\
        ${reads_arg} \\
        --report-file ${prefix}_centrifuge_report.txt \\
        -S ${prefix}_centrifuge_results.txt \\
        ${unaligned} \\
        ${aligned} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_centrifuge_report.txt
    touch ${prefix}_centrifuge_results.txt
    touch ${prefix}.sam
    echo | gzip -n > ${prefix}.unmapped.fastq.gz
    echo | gzip -n > ${prefix}.mapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}