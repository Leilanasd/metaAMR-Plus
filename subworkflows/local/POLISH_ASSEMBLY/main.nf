include { MINIMAP2_ALIGN as MINIMAP2_POLISH } from '../../../modules/nf-core/minimap2/align/main'
include { RACON } from '../../../modules/nf-core/racon/main'

//
// Subworkflow: Polish assembly with Minimap2 + Racon
//

workflow POLISH_ASSEMBLY {
    take:
    ch_input // [ val(meta), path(reads), path(assembly) ]

    main:
    ch_versions = Channel.empty()

    // Fork ch_input into three named sub-channels to avoid queue channel consumption issues
    ch_input
        .multiMap { meta, reads, assembly ->
            reads_ch:    [meta, reads]
            ref_ch:      [[id: assembly.baseName], assembly]
            join_ch:     [meta, reads, assembly]
        }
        .set { ch_split }

    // Minimap2 alignment for polishing — produces PAF for Racon
    MINIMAP2_POLISH(
        ch_split.reads_ch,
        ch_split.ref_ch,
        false, "bai", true, false
    )

    // Prepare input for Racon: [ meta, reads, assembly, paf ]
    ch_racon_input = ch_split.join_ch
        .join(MINIMAP2_POLISH.out.paf)
        .map { meta, reads, assembly, paf ->
            [meta, reads instanceof List ? reads[0] : reads, assembly, paf]
        }

    // Racon polishing — outputs *.fasta.gz
    RACON(ch_racon_input)
    ch_versions = ch_versions.mix(RACON.out.versions)
    ch_versions = ch_versions.mix(Channel.topic('versions'))

    emit:
    polished_assembly = RACON.out.improved_assembly
    versions            = ch_versions
}
