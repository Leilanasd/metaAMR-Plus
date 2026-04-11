/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CENTRIFUGE_CENTRIFUGE } from '../../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../../../modules/nf-core/centrifuge/kreport/main'
include { FILTER_READS_BY_SPECIES } from '../../../modules/local/filter_reads_by_species/main'
include { EXTRACT_FILTERED_READS } from '../../../modules/local/extract_filtered_reads/main'
include { RESFINDER_WITH_SPECIES } from '../../../modules/local/resfinder_with_species/main'

workflow TARGET_SPECIES_AMR {

    take:
    reads
    databases
    resfinder_db
    target_species

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
     * STEP 1 — Prepare reads for Centrifuge
     */
    ch_centrifuge_reads = reads.map { meta, reads_file ->
        def input_reads = reads_file instanceof List ? reads_file.flatten() : [reads_file]
        [meta, input_reads]
    }

    /*
     * STEP 2 — Select only the Centrifuge DB
     */
    ch_centrifuge_db = databases
        .filter { db_meta, db_path -> db_meta.tool == 'centrifuge' }
        .map { db_meta, db_path -> db_path }

    /*
     * STEP 3 — Combine reads with Centrifuge DB
     */
    ch_centrifuge_input = ch_centrifuge_reads.combine(ch_centrifuge_db)

    /*
     * STEP 4 — Run Centrifuge directly
     */
    CENTRIFUGE_CENTRIFUGE(
        ch_centrifuge_input.map { meta, input_reads, db -> [meta, input_reads] },
        ch_centrifuge_input.map { meta, input_reads, db -> db }.first(),
        false,
        false
    )

    ch_versions = ch_versions.mix(CENTRIFUGE_CENTRIFUGE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(CENTRIFUGE_CENTRIFUGE.out.report)

    /*
     * STEP 5 — Optional kreport generation
     */
    CENTRIFUGE_KREPORT(
        CENTRIFUGE_CENTRIFUGE.out.results,
        ch_centrifuge_db.first()
    )

    ch_versions = ch_versions.mix(CENTRIFUGE_KREPORT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(CENTRIFUGE_KREPORT.out.kreport)

    /*
     * STEP 6 — Filter reads by species
     * IMPORTANT: use CENTRIFUGE_CENTRIFUGE.out.report, not kreport
     */
    FILTER_READS_BY_SPECIES(
        CENTRIFUGE_CENTRIFUGE.out.report.join(CENTRIFUGE_CENTRIFUGE.out.results),
        target_species
    )

    ch_filtered_ids = FILTER_READS_BY_SPECIES.out.filtered_read_ids
    ch_species_summary = FILTER_READS_BY_SPECIES.out.species_summary

    /*
     * STEP 7 — Extract target-species reads
     */
    EXTRACT_FILTERED_READS(
        ch_filtered_ids.join(reads)
    )

    ch_filtered_reads = EXTRACT_FILTERED_READS.out.filtered_reads

    /*
     * STEP 8 — Run ResFinder on filtered reads
     */
    ch_resfinder_input = ch_filtered_reads
        .join(ch_species_summary)
        .combine(resfinder_db)

    RESFINDER_WITH_SPECIES(ch_resfinder_input)

    ch_versions = ch_versions.mix(RESFINDER_WITH_SPECIES.out.versions)

    emit:
    filtered_reads  = ch_filtered_reads
    species_summary = ch_species_summary
    resfinder       = RESFINDER_WITH_SPECIES.out.amr_results
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}