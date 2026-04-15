include { FLYE as FLYE_META } from '../../../modules/nf-core/flye/main'
include { QUAST } from '../../../modules/nf-core/quast/main' 


workflow META_ASSEMBLY {

    take:
    reads

    main:
    ch_versions = Channel.empty()
    def mode = "--nano-hq" 
    

    // meta option for metagenomic assembly
    ch_assembly = FLYE_META(
        reads,
        mode,
    ).fasta
    // topic channel: FLYE_META versions

   
    // for assembly quality evaluation
    ch_quast_results = Channel.empty()
    if (!params.skip_quast) {
        QUAST(ch_assembly, [[],[]], [[],[]])
        ch_quast_results = QUAST.out.results
    }
    // topic channel: QUAST versions



    ch_versions = ch_versions.mix(Channel.topic('versions'))

    emit:
    assembly      = ch_assembly
    versions      = ch_versions
    quast_results = ch_quast_results
    
    
}