include { FLYE as FLYE_META } from '../../../modules/nf-core/flye/main'
include { QUAST } from '../../../modules/nf-core/quast/main'

//
// Subworkflow: Metagenome assembly with Flye and optional QUAST quality assessment
//

workflow META_ASSEMBLY {

    take:
    reads // [ val(meta), path(reads) ]

    main:
    ch_versions = Channel.empty()

    // Flye assembly mode — default --nano-hq for modern ONT data (R10.4+, Q20+)
    // Override with --flye_mode "--nano-raw" for older ONT data (R9.4)
    def mode = params.flye_mode
    

    // meta option for metagenomic assembly
    ch_assembly = FLYE_META(
        reads,
        mode,
    ).fasta
    // Filter out failed assemblies and warn user
    ch_assembly
        .filter { meta, fasta ->
            def isEmpty = fasta.size() == 0
            if (isEmpty) log.warn "Sample ${meta.id}: Flye assembly produced no output (possibly no reads above minimum length threshold). Skipping downstream analysis for this sample."
            return !isEmpty
        }
        .set { ch_assembly }

   
    // for assembly quality evaluation
    ch_quast_results = Channel.empty()
    if (!params.skip_quast) {
        QUAST(ch_assembly, [[],[]], [[],[]])
        ch_quast_results = QUAST.out.results
    }



    ch_versions = ch_versions.mix(Channel.topic('versions'))

    emit:
    assembly      = ch_assembly
    versions      = ch_versions
    quast_results = ch_quast_results
    
    
}