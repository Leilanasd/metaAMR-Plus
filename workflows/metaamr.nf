/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PORECHOP_ABI          } from '../modules/nf-core/porechop/abi/main'
include { FILTLONG               } from '../modules/nf-core/filtlong/main'
include { RESFINDER_RUN          } from '../modules/nf-core/resfinder/run/main'
include { AMRFINDERPLUS_RUN      } from '../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_RUN           } from '../modules/nf-core/abricate/run/main'
include { RGI_MAIN               } from '../modules/nf-core/rgi/main/main'
include { PLASMIDFINDER          } from '../modules/nf-core/plasmidfinder/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaamr_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_HOSTREMOVAL          } from '../subworkflows/local/HOSTREMOVAL/main'
include { META_ASSEMBLY              } from '../subworkflows/local/ASSEMBLY/main'
include { QUAST                        } from '../modules/nf-core/quast/main'
include { POLISH_ASSEMBLY            } from '../subworkflows/local/POLISH_ASSEMBLY/main'
include { PREPARE_TOOL_DBS           } from '../subworkflows/local/prepare_tool_dbs/main'
include { HAMRONIZATION              } from '../subworkflows/local/HAMRONIZATION/main'
include { VALIDATE_FASTA             } from '../modules/local/validate_fasta/main'
include { PLASCLASS                  } from '../modules/local/plasclass/main'
include { PROFILING                  } from '../subworkflows/local/PROFILING/main'
include { TARGET_SPECIES_AMR         } from '../subworkflows/local/TARGET_SPECIES_AMR/main'
include { COMBINE_CONTIGS_AND_SPECIES } from '../modules/local/combine_contigs_and_species/main'
include { GENERATE_REPORT              } from '../modules/local/generate_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [
    params.input,
    params.hostremoval_index,
    params.hostremoval_reference,
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

ch_databases = params.databases
    ? Channel.fromPath(params.databases)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [tool: row.tool, db_name: row.db_name, db_params: row.db_params]
            [ meta, file(row.db_path) ]
        }
    : Channel.empty()

ch_reference = params.hostremoval_reference ? file(params.hostremoval_reference) : []
ch_reference_index = params.hostremoval_index ? file(params.hostremoval_index) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAAMR {

    take:
    ch_samplesheet

    main:

    // Validate that assembly is enabled when assembly-dependent tools are requested
    def assembly_required = (
        params.run_rgi          ||
        params.run_amrfinderplus ||
        params.run_abricate     ||
        params.run_plasmidfinder ||
        params.run_plasclass
    )

    if (assembly_required && !params.perform_assembly) {
        error """
        Selected tool(s) require assembled contigs, but --perform_assembly was not enabled.

        Please rerun with:
            --perform_assembly

        Tools requiring assembly:
            RGI, AMRFinderPlus, Abricate, PlasmidFinder, PlasClass
        """
    }

    // Target species mode is read-based and incompatible with assembly
    if (params.target_species && params.perform_assembly) {
        error """
        --target_species is read-based and cannot be used with --perform_assembly.

        Please disable:
            --perform_assembly
            --perform_polish_assembly
        """
    }


    // Host removal requires a reference genome
    if (params.perform_hostremoval && !params.hostremoval_reference) {
        error """
        --perform_hostremoval requires a reference genome.
        Please provide:
            --hostremoval_reference /path/to/reference.fa
        Optionally also provide a pre-built index:
            --hostremoval_index /path/to/index
        """
    }

    // Profiling requires at least one tool
    if (params.run_profiling && !params.run_centrifuge && !params.run_kaiju) {
        log.warn "--run_profiling enabled but neither --run_centrifuge nor --run_kaiju specified. No profiling will be performed."
    }
    if (params.run_centrifuge && !params.run_profiling) {
        log.warn "--run_centrifuge specified but --run_profiling not enabled. Centrifuge will not run. Add --run_profiling to enable profiling."
    }
    if (params.run_kaiju && !params.run_profiling) {
        log.warn "--run_kaiju specified but --run_profiling not enabled. Kaiju will not run. Add --run_profiling to enable profiling."
    }
    if (params.target_species && !params.databases) {
        error """
        --target_species mode requires a Centrifuge database.
        Please provide a databases CSV file with a centrifuge entry:
            --databases database.csv
        Where database.csv contains a line like:
            centrifuge,centrifuge_db,,/path/to/centrifuge_db
        """
    }

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    // Note: .first() is applied to per-sample module version outputs throughout this workflow.
    // This emits exactly one versions.yml entry per tool regardless of sample count,
    // preventing duplicate rows in the MultiQC software versions table.
    // Subworkflows (READS_HOSTREMOVAL, META_ASSEMBLY, etc.) deduplicate versions internally
    // so .first() is not needed on their .out.versions.

    // Prepare tool-specific databases
    PREPARE_TOOL_DBS()
    ch_versions = ch_versions.mix(PREPARE_TOOL_DBS.out.versions)

    //
    // MODULE: FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_samplesheet)
        ch_versions      = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    }

    //
    // MODULE: Porechop + Filtlong (adapter trimming and quality filtering)
    //
    if (params.perform_trim) {
        PORECHOP_ABI(ch_samplesheet, [])

        ch_clipped_reads = PORECHOP_ABI.out.reads
            .map { meta, reads ->
                def readList = reads instanceof List ? reads : [reads]
                def trimmed = readList.find { it.name.endsWith('_trimmed.fastq.gz') } ?: readList[-1]
                [ meta + [single_end: true], trimmed ]
            }

        FILTLONG(ch_clipped_reads.map { meta, reads -> [ meta, [], reads ] })
        ch_processed_reads = FILTLONG.out.reads
            .filter { meta, reads ->
                def isEmpty = reads instanceof List ? reads.every { it.size() == 0 } : reads.size() == 0
                if (isEmpty) log.warn "Sample ${meta.id}: No reads remaining after Filtlong quality filtering. This sample will be skipped."
                return !isEmpty
            }

        ch_versions      = ch_versions.mix(PORECHOP_ABI.out.versions.first())
        ch_versions      = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_ABI.out.log.map { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FILTLONG.out.log.map { it[1] }.ifEmpty([]))
    } else {
        ch_processed_reads = ch_samplesheet
    }

    //
    // SUBWORKFLOW: Host removal
    //
    if (params.perform_hostremoval) {
        READS_HOSTREMOVAL(ch_processed_reads, ch_reference, ch_reference_index)
        ch_hostremoved = READS_HOSTREMOVAL.out.reads
            .filter { meta, reads ->
                def isEmpty = reads instanceof List ? reads.every { it.size() == 0 } : reads.size() == 0
                if (isEmpty) log.warn "Sample ${meta.id}: No reads remaining after host removal. This sample will be skipped."
                return !isEmpty
            }
        ch_versions    = ch_versions.mix(READS_HOSTREMOVAL.out.versions)
    } else {
        ch_hostremoved = ch_processed_reads
    }

    //
    // SUBWORKFLOW: Assembly
    //
    if (params.perform_assembly) {
        META_ASSEMBLY(ch_hostremoved)
        ch_assembly = META_ASSEMBLY.out.assembly
        ch_quast    = META_ASSEMBLY.out.quast_results
        ch_versions = ch_versions.mix(META_ASSEMBLY.out.versions)
    } else {
        ch_assembly = ch_hostremoved
        ch_quast    = Channel.empty()
    }

    //
    // SUBWORKFLOW: Assembly polishing
    //
    if (params.perform_polish_assembly && params.perform_assembly) {
        ch_polish_input = ch_hostremoved
            .join(ch_assembly)
            .map { meta, reads, assembly ->
                [
                    meta,
                    reads instanceof List ? reads[0] : reads,
                    assembly instanceof List ? assembly[0] : assembly
                ]
            }
        POLISH_ASSEMBLY(ch_polish_input)
        ch_final_assembly = POLISH_ASSEMBLY.out.polished_assembly
        ch_versions       = ch_versions.mix(POLISH_ASSEMBLY.out.versions)
    } else {
        if (params.perform_polish_assembly && !params.perform_assembly) {
            log.warn "Assembly polishing requested (--perform_polish_assembly) but --perform_assembly is not enabled. Skipping."
        }
        ch_final_assembly = ch_assembly
    }

    //
    // MODULE: FASTA validation — normalise contig names for all downstream tools
    // Must run before any tool that uses assembly as input so all tools
    // use identical contig names for cross-referencing in the report.
    //
    if (params.perform_assembly) {
        VALIDATE_FASTA(ch_final_assembly)
        ch_validated_assemblies = VALIDATE_FASTA.out.validated_fasta
    } else {
        ch_validated_assemblies = ch_final_assembly
    }

    //
    // MODULE: QUAST — assembly quality on validated assembly
    // Runs on polished assembly if polishing was performed, raw otherwise.
    // Stats reflect the exact contigs used by all downstream tools.
    //
    if (params.perform_assembly) {
        QUAST(ch_validated_assemblies, [[],[]], [[],[]])
        ch_quast    = QUAST.out.results
        ch_versions = ch_versions.mix(QUAST.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(ch_quast.collect { it[1] }.ifEmpty([]))
    }

    //
    // MODULE: ResFinder
    //
    if (params.run_resfinder && !params.target_species) {
        // Select best available input: polished assembly > assembly > host-removed > trimmed reads
        ch_resfinder_input = params.perform_assembly    ? ch_validated_assemblies
                           : params.perform_hostremoval ? ch_hostremoved
                           : ch_processed_reads

        // Route input to correct FASTQ or FASTA argument slot
        ch_resfinder_input = ch_resfinder_input.map { meta, files ->
            def fileList = files instanceof List ? files : [files]
            def isFastq  = fileList.any { f ->
                f.name.toLowerCase().endsWith('.fastq') || f.name.toLowerCase().endsWith('.fastq.gz')
            }
            [ meta, isFastq ? fileList : [], isFastq ? [] : fileList ]
        }

        RESFINDER_RUN(
            ch_resfinder_input,
            [],
            PREPARE_TOOL_DBS.out.resfinder_db
        )
        ch_versions          = ch_versions.mix(RESFINDER_RUN.out.versions.first())
        ch_resfinder_results = RESFINDER_RUN.out.resfinder_results_tab
    } else {
        ch_resfinder_results = Channel.empty()
    }

    //
    // MODULE: Abricate
    //
    if (params.run_abricate) {
        ABRICATE_RUN(ch_validated_assemblies, params.arg_abricate_db)
        ch_versions          = ch_versions.mix(ABRICATE_RUN.out.versions.first())
        ch_abricate_results  = ABRICATE_RUN.out.report
        ch_multiqc_files     = ch_multiqc_files.mix(ABRICATE_RUN.out.report.map { meta, report -> report }.ifEmpty([]))
    } else {
        ch_abricate_results = Channel.empty()
    }

    //
    // MODULE: AMRFinderPlus
    //
    if (params.run_amrfinderplus) {
        AMRFINDERPLUS_RUN(
            ch_validated_assemblies,
            PREPARE_TOOL_DBS.out.amrfinderplus_db
        )
        ch_versions             = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions.first())
        ch_amrfinderplus_results = AMRFINDERPLUS_RUN.out.report
    } else {
        ch_amrfinderplus_results = Channel.empty()
    }

    //
    // MODULE: RGI
    //
    if (params.run_rgi) {
        RGI_MAIN(
            ch_validated_assemblies,
            PREPARE_TOOL_DBS.out.rgi_db,
            []
        )
        ch_versions      = ch_versions.mix(RGI_MAIN.out.versions.first())
        ch_rgi_results   = RGI_MAIN.out.tsv
        ch_multiqc_files = ch_multiqc_files.mix(RGI_MAIN.out.tsv.map { meta, report -> report }.ifEmpty([]))
    } else {
        ch_rgi_results = Channel.empty()
    }



    //
    // MODULE: PlasmidFinder
    //
    if (params.run_plasmidfinder) {
        PLASMIDFINDER(ch_validated_assemblies)
        ch_versions             = ch_versions.mix(PLASMIDFINDER.out.versions.first())
        ch_plasmidfinder_results = PLASMIDFINDER.out.tsv
    } else {
        ch_plasmidfinder_results = Channel.empty()
    }

    //
    // MODULE: PlasClass
    //
    if (params.run_plasclass) {
        PLASCLASS(ch_validated_assemblies, params.plasclass_threshold)
        ch_versions          = ch_versions.mix(PLASCLASS.out.versions.first())
        ch_plasclass_results = PLASCLASS.out.classified
    } else {
        ch_plasclass_results = Channel.empty()
    }

    //
    // SUBWORKFLOW: hAMRonization (harmonise AMR results across tools)
    //
    def any_amr_enabled = params.run_abricate || params.run_amrfinderplus || params.run_rgi

    if (params.run_hamronization && any_amr_enabled) {
        HAMRONIZATION(
            ch_abricate_results,
            ch_amrfinderplus_results,
            ch_rgi_results,
            params.run_rgi && params.download_rgi_db
                ? PREPARE_TOOL_DBS.out.rgi_db_version
                : Channel.value(params.card_version)
        )
        ch_versions = ch_versions.mix(HAMRONIZATION.out.versions.ifEmpty([]))
    } else if (params.run_hamronization) {
        log.warn "hAMRonization requested but no supported AMR tools are enabled (Abricate, AMRFinderPlus, RGI)."
    }

    //
    // SUBWORKFLOW: Taxonomic profiling
    //
    if (params.run_profiling && !params.target_species) {
        // ch_validated_assemblies contains normalised contigs if --perform_assembly,
        // or reads otherwise — Centrifuge/Kaiju accept both FASTA and FASTQ
        ch_profiling_input = ch_validated_assemblies.map { meta, assembly -> [ meta, [assembly] ] }

        PROFILING(ch_profiling_input, ch_databases)
        ch_versions      = ch_versions.mix(PROFILING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PROFILING.out.multiqc_files.map { it[1] }.ifEmpty([]))
        if (!params.skip_krona) {
            ch_multiqc_files = ch_multiqc_files.mix(PROFILING.out.krona_html.map { it[1] }.ifEmpty([]))
        }

        if (params.perform_assembly && params.run_centrifuge) {
            COMBINE_CONTIGS_AND_SPECIES(
                PROFILING.out.centrifuge_results.join(PROFILING.out.centrifuge_report)
            )
            ch_versions = ch_versions.mix(COMBINE_CONTIGS_AND_SPECIES.out.versions)
        }

    }

    //
    // SUBWORKFLOW: Target species AMR mode (read-based, no assembly)
    //
    

    if (params.target_species && params.run_resfinder) {
        log.warn "--run_resfinder is ignored in --target_species mode. ResFinder runs automatically on extracted target species reads via TARGET_SPECIES_AMR."
    }

    if (params.target_species) {
        ch_reads_for_target = params.perform_hostremoval ? ch_hostremoved
                        : params.perform_trim        ? ch_processed_reads
                        : ch_samplesheet

        TARGET_SPECIES_AMR(
            ch_reads_for_target,
            ch_databases,
            PREPARE_TOOL_DBS.out.resfinder_db,
            params.target_species
        )
        ch_versions      = ch_versions.mix(TARGET_SPECIES_AMR.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TARGET_SPECIES_AMR.out.multiqc_files.map { it[1] }.ifEmpty([]))
    }

    //
    // Collate software versions
    //
    softwareVersionsToYAML(ch_versions.mix(Channel.topic('versions')))
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:     'nf_core_pipeline_software_mqc_versions.yml',
            sort:     true,
            newLine:  true
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    // Note: MULTIQC module takes a single tuple input — config/logo passed inline

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files    = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
    )

    if (params.perform_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(READS_HOSTREMOVAL.out.mqc.map { meta, stats -> stats }.ifEmpty([]))
    }
    if (params.perform_assembly) {
    }

    MULTIQC(
        ch_multiqc_files.collect().map { files ->
            [
                [id: 'multiqc'],
                files,
                params.multiqc_config
                    ? [file("$projectDir/assets/multiqc_config.yml"), file(params.multiqc_config)]
                    : [file("$projectDir/assets/multiqc_config.yml")],
                params.multiqc_logo ? [file(params.multiqc_logo)] : [],
                [],
                []
            ]
        }
    )

    //
    // MODULE: HTML Report
    //

    GENERATE_REPORT(
        MULTIQC.out.report.map { meta, report -> file(params.outdir).toAbsolutePath().toString() },
        params.run_name ?: workflow.runName
    )

    
    emit:
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }.toList()
    versions       = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
