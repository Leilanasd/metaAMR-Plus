/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PORECHOP_PORECHOP      } from '../modules/nf-core/porechop/main'
include { FILTLONG               } from '../modules/nf-core/filtlong/main'
include { RESFINDER_RUN } from '../modules/nf-core/resfinder/run/main'
include { AMRFINDERPLUS_RUN } from '../modules/nf-core/amrfinderplus/run/main' 
include { AMRFINDERPLUS_UPDATE } from '../modules/nf-core/amrfinderplus/update/main' 
include { ABRICATE_RUN } from '../modules/nf-core/abricate/run/main' 
include { RGI_CARDANNOTATION } from '../modules/nf-core/rgi/cardannotation/main' 
include { RGI_MAIN } from '../modules/nf-core/rgi/main/main' 
include { PLASMIDFINDER_RUN } from '../modules/nf-core/plasmidfinder/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaamr_pipeline'
include { CENTRIFUGE_CENTRIFUGE } from '../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../modules/nf-core/centrifuge/kreport/main'


// Check input path parameters to see if they exist

def checkPathParamList = [ params.input, 
                           params.hostremoval_index,
                           params.hostremoval_reference,
                         ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

if (params.input) {
    ch_input = Channel.fromPath(params.input)
                      .splitCsv(header:true, sep:',')
                      .map { row -> 
                          def meta = [id: row.sample, single_end: true]
                          def reads = file(row.fastq_1, checkIfExists: true)
                          return [meta, reads]
                      }
} else {
    error("Input samplesheet not specified")
}

// Check if databases file is provided and create a channel
ch_databases = params.databases ? Channel.fromPath(params.databases)
    .splitCsv(header:true, sep:',')
    .map { row -> 
        def meta = [:]
        meta.tool = row.tool
        meta.db_name = row.db_name
        meta.db_params = row.db_params
        [ meta, file(row.db_path) ]
    } : Channel.empty()

// Log information about databases
if (params.databases) {
    log.info "Using provided database file: ${params.databases}"
} else {
    log.info "No database file provided. Tools will use default databases or prepare them as needed."
}


if (params.hostremoval_reference) { 
    ch_reference = file(params.hostremoval_reference) 
}
if (params.hostremoval_index) { 
    ch_reference_index = file(params.hostremoval_index) 
} else { 
    ch_reference_index = [] 
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include {READS_HOSTREMOVAL       } from '../subworkflows/local/HOSTREMOVAL'
include {META_ASSEMBLY      } from '../subworkflows/local/ASSEMBLY'
include {POLISH_ASSEMBLY    } from '../subworkflows/local/POLISH_ASSEMBLY'
include { PREPARE_TOOL_DBS } from '../subworkflows/local/prepare_tool_dbs'
include { HAMRONIZATION } from '../subworkflows/local/HAMRONIZATION'
include { VALIDATE_FASTA } from '../modules/local/validate_fasta'
include { PLASCLASS } from '../modules/local/plasclass'
include { PLASCLASS_POSTPROCESS } from '../modules/local/plasclass_postprocess.nf'
include { PROFILING } from '../subworkflows/local/PROFILING'
include { TARGET_SPECIES_AMR } from '../subworkflows/local/TARGET_SPECIES_AMR'
include { COMBINE_CONTIGS_AND_SPECIES } from '../modules/local/combine_contigs_and_species'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAAMR {

    take:
    ch_samplesheet 
    

    main:
    /*
 * Validate assembly requirements
 * Some tools require assembled contigs
 */
    def assembly_required = (
        params.run_rgi ||
        params.run_amrfinderplus ||
        params.run_abricate ||
        params.run_plasmidfinder ||
        params.run_plasclass
    )

    if (assembly_required && !params.perform_assembly) {
        error """
        Selected tool(s) require assembled contigs, but --perform_assembly was not enabled.

        Please rerun the pipeline with:
            --perform_assembly

        Tools requiring assembly:
            RGI
            AMRFinderPlus
            Abricate
            PlasmidFinder
            PlasClass
        """
    }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
   

    
    // Prepare tool-specific databases
    PREPARE_TOOL_DBS()
    
    def ch_rgi_db_final = PREPARE_TOOL_DBS.out.rgi_db
    //
    // MODULE: Run FastQC
    //
    if (params.run_fastqc) {
        FASTQC(
            ch_samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }
    //
    // MODULE: Run PORECHOPS & FILTLONG
    //
    
    if (params.perform_trim) {
        PORECHOP_PORECHOP(
            ch_samplesheet
        )
        ch_clipped_reads = PORECHOP_PORECHOP.out.reads
            .map { meta, reads -> 
                def porechopped_reads = reads.findAll { it.name.contains('porechopped') } 
                [ meta + [single_end: true], porechopped_reads ] }
    
        ch_processed_reads = FILTLONG ( ch_clipped_reads.map { meta, reads -> [ meta, [], reads ] } ).reads

        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(
            FILTLONG.out.log
                .map { meta, log -> log }
                .ifEmpty([])
        )
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_PORECHOP.out.log.map{ it[1] } )
    } else {
        ch_processed_reads = ch_samplesheet
    }
    
    /*
        SUBWORKFLOW: HOST REMOVAL
    */
    if ( params.perform_hostremoval ) {
        ch_hostremoved = READS_HOSTREMOVAL(
            ch_processed_reads, 
            ch_reference,   
            ch_reference_index         
        ).reads
        ch_versions = ch_versions.mix(READS_HOSTREMOVAL.out.versions)
    } else {
        ch_hostremoved = ch_processed_reads
    }
  

    /*
    SUBWORKFLOW: ASSEMBLY
    */

    if ( params.perform_assembly) {
        ch_assembly = META_ASSEMBLY(ch_hostremoved).ch_assembly   
        ch_versions = ch_versions.mix(META_ASSEMBLY.out.ch_versions)
        ch_quast = META_ASSEMBLY.out.quast_results 
    } else {
        ch_assembly = ch_hostremoved
        ch_quast = Channel.empty()
    }

   // Polish assembly 
    if (params.perform_polish_assembly && params.perform_assembly) {
        ch_polish_input = ch_hostremoved.join(ch_assembly).map { meta, reads, assembly ->
            [meta, reads instanceof List ? reads[0] : reads, assembly instanceof List ? assembly[0] : assembly]
        }
        POLISH_ASSEMBLY(ch_polish_input)
        ch_final_polished_assembly = POLISH_ASSEMBLY.out.polished_assembly_1
        ch_versions = ch_versions.mix(POLISH_ASSEMBLY.out.versions)
    } else {
        ch_final_polished_assembly = ch_assembly
    }

    //  assemblies for downstream tools
    ch_assembly_for_arg = ch_final_polished_assembly.mix(ch_hostremoved)
        .groupTuple()
        .map { meta, assemblies -> [meta, assemblies.find { it != null } ?: meta.reads] }


    if (params.run_resfinder) {
        log.info "Running ResFinder"

    // Use polished assembly if available, otherwise fallback
    ch_resfinder_input = params.perform_polish_assembly ? ch_final_polished_assembly
                         : params.perform_assembly ? ch_assembly
                         : params.perform_hostremoval ? ch_hostremoved
                         : ch_processed_reads  // Use processed reads if nothing else

    // Ensure correct format handling (FASTQ vs. FASTA)
    ch_resfinder_input = ch_resfinder_input.map { meta, files -> 
        def isFastq = files.any { file ->
            file.name.toLowerCase().endsWith('.fastq') || file.name.toLowerCase().endsWith('.fastq.gz')
        }
        def fastq = isFastq ? files : []
        def fasta = isFastq ? [] : files
        return [meta, fastq, fasta]
    }

    // Combine with ResFinder database
    ch_resfinder_input = ch_resfinder_input
        .combine(PREPARE_TOOL_DBS.out.resfinder_db)
        .map { meta, fastq, fasta, db -> 
            return [meta, fastq, fasta, db, []]  // Maintain correct argument structure
        }

    // Run ResFinder
        RESFINDER_RUN(ch_resfinder_input)

    // Capture versions & outputs
        ch_versions = ch_versions.mix(RESFINDER_RUN.out.versions.first())
        ch_resfinder = RESFINDER_RUN.out.all_outputs
    }
    
    if (params.run_abricate) {
      
        log.info "Running Abricate"

        ch_abricate_input = ch_final_polished_assembly
        
        ABRICATE_RUN(
            ch_abricate_input,
            params.arg_abricate_db
        )
    
        ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ABRICATE_RUN.out.report.map { meta, report -> report })
        ABRICATE_RUN.out.report.view { meta, report -> 
            "Abricate outputs for ${meta.id}: ${report.getName()}"
        }
    }
    
    if (params.run_amrfinderplus) {
        log.info "Running AMRFinderPlus"

        ch_amrfinderplus_input = ch_final_polished_assembly
        ch_amrfinderplus_input.combine(PREPARE_TOOL_DBS.out.amrfinderplus_db)
        

        AMRFINDERPLUS_RUN(
            ch_amrfinderplus_input,
            PREPARE_TOOL_DBS.out.amrfinderplus_db,
            
        )

        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)
        AMRFINDERPLUS_RUN.out.report.view { meta, report -> 
            "AMRFinderPlus outputs for ${meta.id}: ${report.getName()}"
        }
        ch_multiqc_files = ch_multiqc_files.mix(AMRFINDERPLUS_RUN.out.report.collect{it[1]}.ifEmpty([]))
        
    }


    if (params.run_rgi) {
        log.info "Running RGI"

        ch_rgi_input = ch_final_polished_assembly
            .map { meta, assembly -> [ meta, assembly ] }
            
        ch_rgi_input
            .combine(ch_rgi_db_final)
            .map { meta, assembly, db -> [ meta, assembly, db ] }  
            .set { ch_rgi_combined_input }

        RGI_MAIN(
            ch_rgi_combined_input,
            []  
        )

        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)
    
        RGI_MAIN.out.tsv
            .view { meta, tsv -> "RGI outputs for ${meta.id}: ${tsv.getName()}" }
            .ifEmpty { log.warn "No output from RGI_MAIN process" }

        ch_multiqc_files = ch_multiqc_files.mix(
            RGI_MAIN.out.tsv.map { meta, tsv -> tsv }.ifEmpty([])
        )
    }    
    
    //  if FASTA validation is needed
    def run_validate_fasta = params.run_plasmidfinder || params.run_plasclass

    // Run FASTA validation only if PlasmidFinder or PlasClass is enabled
    if (run_validate_fasta) {
        log.info "Validating FASTA files"
        VALIDATE_FASTA(ch_final_polished_assembly)
        ch_validated_assemblies = VALIDATE_FASTA.out.validated_fasta
    } else {
        log.info "Skipping FASTA validation"
        ch_validated_assemblies = ch_final_polished_assembly
    }


    if (params.run_plasmidfinder) {
        log.info "Running PlasmidFinder"

    // Run PlasmidFinder
        PLASMIDFINDER_RUN(ch_validated_assemblies, PREPARE_TOOL_DBS.out.plasmidfinder_db.first())

        ch_versions = ch_versions.mix(PLASMIDFINDER_RUN.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(
            PLASMIDFINDER_RUN.out.tsv.collect { it[1] }.ifEmpty([])
        )
        
    }

    // Run PlasClass
    if (params.run_plasclass) {
        log.info "Running PlasClass"

        PLASCLASS(ch_validated_assemblies)

        PLASCLASS_POSTPROCESS(PLASCLASS.out.report)

        ch_versions = ch_versions.mix(PLASCLASS.out.versions)
        ch_versions = ch_versions.mix(PLASCLASS_POSTPROCESS.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(
            PLASCLASS_POSTPROCESS.out.classified.collect { it[1] }.ifEmpty([])
        )
    }
    
  
    // Collect results from AMR/plasmid tools
    ch_abricate_results       = params.run_abricate ? ABRICATE_RUN.out.report : Channel.empty()
    ch_amrfinderplus_results  = params.run_amrfinderplus ? AMRFINDERPLUS_RUN.out.report : Channel.empty()
    ch_rgi_results            = params.run_rgi ? RGI_MAIN.out.tsv : Channel.empty()
    ch_resfinder_results = params.run_resfinder ? RESFINDER_RUN.out.txt : Channel.empty()
    ch_plasclass_results      = params.run_plasclass ? PLASCLASS.out.report : Channel.empty()
   

    // Run HAMRONIZATION
    def any_amr_enabled = params.run_abricate || params.run_amrfinderplus || params.run_rgi

    if (params.run_hamronization && any_amr_enabled) {
        HAMRONIZATION (
            ch_abricate_results,
            ch_amrfinderplus_results,
            ch_rgi_results
        )
        ch_versions = ch_versions.mix(HAMRONIZATION.out.versions.ifEmpty([]))
        
    } else if (params.run_hamronization) {
        log.warn "HAMRONIZATION requested but no AMR tools are enabled."
    } else {
        log.info "Skipping HAMRONIZATION: --run_hamronization not enabled."
    }
    
    // Run Profiling
    if (params.run_profiling && !params.target_species) {
    ch_profiling_input = ch_final_polished_assembly.map { meta, assembly -> 
        [meta, [assembly]]  //  
      
    }
    ch_profiling_input.view { "Profiling input: $it" }
    
    databases_ch = Channel.fromPath(params.databases)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def db_meta = [:]
            db_meta.tool = row.tool
            db_meta.db_name = row.db_name
            db_meta.db_params = row.db_params
            [db_meta, file(row.db_path)]
        }

    PROFILING(
        ch_profiling_input,
        databases_ch
    )
    ch_versions = ch_versions.mix(PROFILING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PROFILING.out.multiqc_files.map { it[1] }.ifEmpty([]))
    

    if (params.perform_assembly && params.run_centrifuge) {
        COMBINE_CONTIGS_AND_SPECIES(
            PROFILING.out.centrifuge_results.join(PROFILING.out.centrifuge_report)
        )
    } else {
        log.info "Skipping COMBINE_CONTIGS_AND_SPECIES: requires both --perform_assembly and --run_centrifuge."
    }

    }
    ch_centrifuge_species = (params.run_profiling && !params.target_species && params.run_centrifuge && params.perform_assembly) \
        ? COMBINE_CONTIGS_AND_SPECIES.out.contigs_species_table \
        : Channel.empty()
    ch_centrifuge_results = ch_centrifuge_species
    

    // Run Target Species
    if (params.target_species && params.perform_assembly) {
        error """
        --target_species is READ-BASED and cannot be used with assembly.

        Please disable:
            --perform_assembly
            --perform_polish_assembly

        Target species mode only works on:
            raw / trimmed / host-removed reads
        """
    }
    ch_reads_for_target =
        params.perform_hostremoval ? ch_hostremoved :
        params.perform_trim       ? ch_processed_reads :
                                    ch_samplesheet
    
    if (params.target_species) {

        TARGET_SPECIES_AMR(
            ch_reads_for_target,
            ch_databases,
            PREPARE_TOOL_DBS.out.resfinder_db,
            params.target_species
        )

        ch_versions = ch_versions.mix(TARGET_SPECIES_AMR.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TARGET_SPECIES_AMR.out.multiqc_files.map { it[1] }.ifEmpty([]))

    }



    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    
    
    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    
    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    if (params.perform_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(READS_HOSTREMOVAL.out.mqc.collect{it[1]}.ifEmpty([]))
    }
    if (params.perform_assembly) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_quast.collect{it[1]}.ifEmpty([]))
    }
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() 
    versions       = ch_versions.ifEmpty(null)               

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/