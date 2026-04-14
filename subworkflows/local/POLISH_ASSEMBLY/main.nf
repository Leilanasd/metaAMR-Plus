include { MINIMAP2_ALIGN as MINIMAP2_POLISH } from '../../../modules/nf-core/minimap2/align/main'
include { RACON } from '../../../modules/nf-core/racon/main'

workflow POLISH_ASSEMBLY {
    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    // Ensure input assemblies are .gz compressed
    ch_prepped_assembly = ch_input.map { meta, reads, assembly_file ->
        def final_assembly = file("${workDir}/${meta.id}.assembly.fasta.gz")

        if (!assembly_file.name.endsWith(".gz")) {
            "gzip -c ${assembly_file} > ${final_assembly}".execute().waitFor()
        } else {
            assembly_file.copyTo(final_assembly)
        }

        assert final_assembly.exists() : "Assembly file ${final_assembly.toAbsolutePath()} does not exist"

        return [meta, reads, final_assembly]
    }

    // Minimap2 alignment for polishing
    MINIMAP2_POLISH(
        ch_prepped_assembly.map { meta, reads, assembly -> [meta, reads] },
        ch_prepped_assembly.map { meta, reads, assembly -> [[id: assembly.baseName], assembly] },
        false, "bai", true, false
    )
    // topic channel: MINIMAP2_POLISH versions

    ch_minimap2_output = MINIMAP2_POLISH.out.paf
        .map { it -> sleep(100); it }

    // Prepare input for Racon
    ch_racon_input = ch_prepped_assembly
        .join(ch_minimap2_output)
        .map { meta, reads, assembly, paf ->
            [meta, reads instanceof List ? reads[0] : reads, assembly, paf]
        }

    // Racon polishing
    RACON(ch_racon_input)
    ch_versions = ch_versions.mix(RACON.out.versions)

    ch_racon_gzipped = RACON.out.improved_assembly.map { meta, racon_fasta ->
        def racon_gz = file("${workDir}/${meta.id}.racon.fasta.gz")

        if (!racon_fasta.name.endsWith(".gz")) {
            "gzip -c ${racon_fasta} > ${racon_gz}".execute().waitFor()
        } else {
            racon_fasta.copyTo(racon_gz)
        }

        assert racon_gz.exists() : "Racon output file ${racon_gz.toAbsolutePath()} does not exist"

        return [meta, racon_gz]
    }

    ch_versions = ch_versions.mix(Channel.topic('versions'))

    emit:
    polished_assembly_1 = ch_racon_gzipped
    versions            = ch_versions
}