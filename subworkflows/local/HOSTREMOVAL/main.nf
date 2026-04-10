//
// Remove host reads via alignment and export off-target reads
//

include { MINIMAP2_INDEX             } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN             } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_VIEW              } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ             } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_INDEX             } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS             } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_SORT              } from '../../../modules/nf-core/samtools/sort/main'

workflow READS_HOSTREMOVAL {
    take:
    reads     // [ [ meta ], [ reads ] ]
    reference // /path/to/fasta
    index    // /path/to/index

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( !params.hostremoval_index ) {
        ch_minimap2_index = MINIMAP2_INDEX ( [ [], reference ] ).index.map { it[1] }
        // topic channel: MINIMAP2_INDEX versions
    } else {
        ch_minimap2_index = index
    }
     // Pass FILTLONG processed reads to the alignment step for host removal
    MINIMAP2_ALIGN ( reads , ch_minimap2_index.map { index -> [[id:"reference"], index] }, true, "bai", false, false)
    // topic channel: MINIMAP2_ALIGN versions
    ch_minimap2_mapped = MINIMAP2_ALIGN.out.bam
        .map {
            meta, reads ->
                [ meta, reads, [] ]
        }

    // Generate unmapped reads FASTQ for downstream taxprofiling
    SAMTOOLS_VIEW ( ch_minimap2_mapped, [[],[],[]], [[],[]], [[],[]], "bai" )
    // topic channel: SAMTOOLS_VIEW versions

// Filter for unmapped BAM files ending in "_unmapped.bam" and create a channel for SAMTOOLS_FASTQ
    SAMTOOLS_VIEW.out.bam
    .filter { meta, bam -> bam.name.contains(".unmapped.bam") }
    .set { ch_unmapped_bam }

    // Convert unmapped BAM to FASTQ
    SAMTOOLS_FASTQ ( ch_unmapped_bam, false )
    // topic channel: SAMTOOLS_FASTQ versions

    // Indexing whole BAM for host removal statistics
    SAMTOOLS_INDEX ( MINIMAP2_ALIGN.out.bam )
    // topic channel: SAMTOOLS_INDEX versions

    bam_bai = MINIMAP2_ALIGN.out.bam
        .join(SAMTOOLS_INDEX.out.index)

    SAMTOOLS_STATS ( bam_bai, [[],[],[]] )
    // topic channel: SAMTOOLS_STATS versions
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_STATS.out.stats )

    ch_versions = ch_versions.mix(Channel.topic('versions'))

    emit:
    stats    = SAMTOOLS_STATS.out.stats     //channel: [val(meta), [reads  ] ]
    reads    = SAMTOOLS_FASTQ.out.other
    versions = ch_versions                 // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
