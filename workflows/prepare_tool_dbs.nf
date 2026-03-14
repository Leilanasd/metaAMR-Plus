include { DOWNLOAD_DB as RESFINDER_DB_DOWNLOAD } from '../modules/local/download_db'
include { DOWNLOAD_DB as RGI_DB_DOWNLOAD } from '../modules/local/download_db'
include { RESFINDER_INDEX } from '../modules/local/RESFINDER_INDEX'
include { AMRFINDERPLUS_UPDATE } from '../modules/nf-core/amrfinderplus/update/main'
include { DOWNLOAD_DB as PLASMIDFINDER_DB_DOWNLOAD } from '../modules/local/download_db'

/*
 * Read database path for one tool from database.csv
 * Expected CSV header: tool,db_name,db_params,db_path
 */
def get_db_path_from_csv(csv_path, tool_name) {
    if (!csv_path) {
        return null
    }

    def csv_file = file(csv_path)
    if (!csv_file.exists()) {
        log.debug "Database CSV file not found: ${csv_path}"
        return null
    }

    def lines = csv_file.readLines()
    if (lines.size() < 2) {
        log.debug "Database CSV file is empty or has no data rows: ${csv_path}"
        return null
    }

    def header = lines[0].split(',', -1)*.trim()
    def tool_idx = header.indexOf('tool')
    def path_idx = header.indexOf('db_path')

    if (tool_idx == -1 || path_idx == -1) {
        log.warn "Database CSV missing required columns 'tool' or 'db_path': ${csv_path}"
        return null
    }

    for (line in lines.drop(1)) {
        if (!line?.trim()) {
            continue
        }

        def cols = line.split(',', -1)*.trim()
        if (cols.size() <= Math.max(tool_idx, path_idx)) {
            continue
        }

        if (cols[tool_idx] == tool_name && cols[path_idx]) {
            log.debug "Found ${tool_name} database in CSV: ${cols[path_idx]}"
            return cols[path_idx]
        }
    }

    return null
}

workflow PREPARE_TOOL_DBS {

    main:

    /*
     * Resolve DB paths from database.csv
     */
    def csv_resfinder_db     = get_db_path_from_csv(params.databases, 'resfinder')
    def csv_rgi_db           = get_db_path_from_csv(params.databases, 'rgi')
    def csv_amrfinderplus_db = get_db_path_from_csv(params.databases, 'amrfinderplus')
    def csv_plasmidfinder_db = get_db_path_from_csv(params.databases, 'plasmidfinder')

    /*
     * ResFinder
     */
    if (csv_resfinder_db) {
        ch_resfinder_db_final = Channel.value(file(csv_resfinder_db, checkIfExists: true))
    } else if (params.resfinder_db) {
        ch_resfinder_db_final = Channel.value(file(params.resfinder_db, checkIfExists: true))
    } else if (params.download_resfinder_db) {
        RESFINDER_DB_DOWNLOAD(Channel.of('resfinder'))
        ch_resfinder_db_final = RESFINDER_INDEX(RESFINDER_DB_DOWNLOAD.out.db)
            .map { db_files -> file(db_files[0]).parent }
    } else {
        ch_resfinder_db_final = Channel.empty()
    }

    /*
     * RGI
     */
    if (csv_rgi_db) {
        ch_rgi_db_final = Channel.value(file(csv_rgi_db, checkIfExists: true))
    } else if (params.rgi_db) {
        ch_rgi_db_final = Channel.value(file(params.rgi_db, checkIfExists: true))
    } else if (params.download_rgi_db) {
        RGI_DB_DOWNLOAD(Channel.of('rgi'))
        ch_rgi_db_final = RGI_DB_DOWNLOAD.out.db
    } else {
        ch_rgi_db_final = Channel.empty()
    }

    /*
     * AMRFinderPlus
     */
    if (csv_amrfinderplus_db) {
        ch_amrfinderplus_db_final = Channel.value(file(csv_amrfinderplus_db, checkIfExists: true))
    } else if (params.amrfinderplus_db) {
        ch_amrfinderplus_db_final = Channel.value(file(params.amrfinderplus_db, checkIfExists: true))
    } else if (params.download_amrfinderplus_db) {
        AMRFINDERPLUS_UPDATE()
        ch_amrfinderplus_db_final = AMRFINDERPLUS_UPDATE.out.db
    } else {
        ch_amrfinderplus_db_final = Channel.empty()
    }

    /*
     * PlasmidFinder
     */
    if (csv_plasmidfinder_db) {
        ch_plasmidfinder_db_final = Channel.value(file(csv_plasmidfinder_db, checkIfExists: true))
    } else if (params.plasmidfinder_db) {
        ch_plasmidfinder_db_final = Channel.value(file(params.plasmidfinder_db, checkIfExists: true))
    } else if (params.download_plasmidfinder_db) {
        PLASMIDFINDER_DB_DOWNLOAD(Channel.of('plasmidfinder'))
        ch_plasmidfinder_db_final = PLASMIDFINDER_DB_DOWNLOAD.out.db
    } else {
        ch_plasmidfinder_db_final = Channel.empty()
    }

    emit:
    resfinder_db     = ch_resfinder_db_final
    rgi_db           = ch_rgi_db_final
    amrfinderplus_db = ch_amrfinderplus_db_final
    plasmidfinder_db = ch_plasmidfinder_db_final
}