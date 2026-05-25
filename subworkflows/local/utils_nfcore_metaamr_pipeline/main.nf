//
// Subworkflow with functionality specific to the nf-core/metaamr pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    
    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )
    

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, reads ->
            return [ meta + [ single_end:true ], reads ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report.toList()
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    // Assembly-dependent tools need --perform_assembly
    def assembly_required = (
        params.run_rgi          ||
        params.run_amrfinderplus ||
        params.run_abricate     ||
        params.run_plasmidfinder ||
        params.run_plasclass
    )
    if (assembly_required && !params.perform_assembly) {
        error("Assembly-dependent tools (RGI, AMRFinderPlus, Abricate, PlasmidFinder, PlasClass) require --perform_assembly to be enabled.")
    }

    // target_species is incompatible with assembly
    if (params.target_species && params.perform_assembly) {
        error("--target_species is read-based and cannot be used with --perform_assembly.")
    }

    // target_species needs a database CSV
    if (params.target_species && !params.databases) {
        error("--target_species mode requires a Centrifuge database provided via --databases CSV.")
    }

    // Host removal needs a reference
    if (params.perform_hostremoval && !params.hostremoval_reference) {
        error("--perform_hostremoval requires a host reference genome via --hostremoval_reference.")
    }

    // Profiling tool flags without --run_profiling
    if (params.run_centrifuge && !params.run_profiling) {
        log.warn "--run_centrifuge specified but --run_profiling not enabled. Centrifuge will not run. Add --run_profiling."
    }
    if (params.run_kaiju && !params.run_profiling) {
        log.warn "--run_kaiju specified but --run_profiling not enabled. Kaiju will not run. Add --run_profiling."
    }
    if (params.run_profiling && !params.run_centrifuge && !params.run_kaiju) {
        log.warn "--run_profiling enabled but neither --run_centrifuge nor --run_kaiju specified. No profiling will be performed."
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        params.perform_trim                                       ? "Porechop_ABI (Bonenfant et al. 2023), Filtlong (Wick 2021),"  : "",
        params.perform_hostremoval                                ? "Minimap2 (Li 2018), SAMtools (Danecek et al. 2021)," : "",
        params.perform_assembly                                   ? "Flye (Kolmogorov et al. 2019),"                      : "",
        params.perform_assembly             ? "QUAST (Gurevich et al. 2013),"                       : "",
        params.perform_polish_assembly                            ? "Racon (Vaser et al. 2017),"                          : "",
        params.run_abricate                                       ? "Abricate (Seemann 2020),"                            : "",
        params.run_rgi                                            ? "RGI (Alcock et al. 2023),"                           : "",
        params.run_amrfinderplus                                  ? "AMRFinderPlus (Feldgarden et al. 2021),"             : "",
        params.run_resfinder || params.target_species             ? "ResFinder (Bortolaia et al. 2020),"                  : "",
        params.run_hamronization                                  ? "hAMRonization (Maguire et al. 2023),"                : "",
        params.run_profiling && params.run_centrifuge             ? "Centrifuge (Kim et al. 2016),"                       : "",
        params.run_profiling && params.run_kaiju                  ? "Kaiju (Menzel et al. 2016),"                         : "",
        (params.run_profiling && !params.skip_krona) || params.target_species ? "Krona (Ondov et al. 2011),"             : "",
        params.run_plasmidfinder                                  ? "PlasmidFinder (Carattoli et al. 2014),"              : "",
        params.run_plasclass                                      ? "PlasClass (Pellow et al. 2020),"                     : "",
        "MultiQC (Ewels et al. 2016)."
    ].findAll { it }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>",
        params.perform_trim ? "<li>Bonenfant Q et al., (2023) Porechop_ABI: discovering unknown adapters in Oxford Nanopore Technology sequencing reads for downstream trimming. Bioinformatics Advances, 3(1):vbac085. doi: 10.1093/bioadv/vbac085.</li>" : "",
        params.perform_trim ? "<li>Wick RR, (2021) Filtlong, URL: https://github.com/rrwick/Filtlong.</li>" : "",
        params.perform_hostremoval ? "<li>Li H, (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.</li>" : "",
        params.perform_hostremoval ? "<li>Danecek P et al., (2021) Twelve years of SAMtools and BCFtools. Gigascience, 10(2):giab008. doi: 10.1093/gigascience/giab008.</li>" : "",
        params.perform_assembly ? "<li>Kolmogorov M et al., (2019) Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol, 37:540-546. doi: 10.1038/s41587-019-0072-8.</li>" : "",
        params.perform_assembly ? "<li>Gurevich A et al., (2013) QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29(8):1072-1075. doi: 10.1093/bioinformatics/btt086.</li>" : "",
        params.perform_polish_assembly ? "<li>Vaser R et al., (2017) Fast and accurate de novo genome assembly from long uncorrected reads. Genome Res, 27(5):737-746. doi: 10.1101/gr.214270.116.</li>" : "",
        params.run_abricate ? "<li>Seemann T, (2020) Abricate, URL: https://github.com/tseemann/abricate.</li>" : "",
        params.run_rgi ? "<li>Alcock BP et al., (2023) CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Res, 51(D1):D690-D699. doi: 10.1093/nar/gkac920.</li>" : "",
        params.run_amrfinderplus ? "<li>Feldgarden M et al., (2021) AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep, 11:12728. doi: 10.1038/s41598-021-91456-0.</li>" : "",
        params.run_resfinder || params.target_species ? "<li>Bortolaia V et al., (2020) ResFinder 4.0 for predictions of phenotypes from genotypes. J Antimicrob Chemother, 75(12):3491-3500. doi: 10.1093/jac/dkaa345.</li>" : "",
        params.run_hamronization ? "<li>Maguire M et al., (2023) hAMRonization: Enhancing antimicrobial resistance prediction using a comprehensive tool for standardization of AMR gene detection results. J Antimicrob Chemother. doi: 10.1093/jac/dkad327.</li>" : "",
        params.run_profiling && params.run_centrifuge || params.target_species ? "<li>Kim D et al., (2016) Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome Res, 26(12):1721-1729. doi: 10.1101/gr.210641.116.</li>" : "",
        params.run_profiling && params.run_kaiju ? "<li>Menzel P et al., (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat Commun, 7:11257. doi: 10.1038/ncomms11257.</li>" : "",
        (params.run_profiling && !params.skip_krona) || params.target_species ? "<li>Ondov BD et al., (2011) Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12:385. doi: 10.1186/1471-2105-12-385.</li>" : "",
        params.run_plasmidfinder ? "<li>Carattoli A et al., (2014) In Silico Detection and Typing of Plasmids using PlasmidFinder and Plasmid Multilocus Sequence Typing. Antimicrob Agents Chemother, 58(7):3895-3903. doi: 10.1128/AAC.02412-14.</li>" : "",
        params.run_plasclass ? "<li>Pellow D et al., (2020) PlasClass improves plasmid sequence classification. PLOS Comput Biol, 16(4):e1007781. doi: 10.1371/journal.pcbi.1007781.</li>" : "",
        "<li>Ewels P et al., (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19):3047-3048. doi: 10.1093/bioinformatics/btw354.</li>"
    ].findAll { it }.join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

