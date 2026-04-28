process DOWNLOAD_DB {
    tag "$tool"
    label 'process_high'
    label 'error_retry'
    publishDir "${params.outdir}/databases", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.equals('versions.yml')) return null
            if (params.save_databases)           return "${tool}/$filename"
            return null
        }
        
    input:
    val tool

    output:
    path "${tool}_db", emit: db
    path "${tool}_db.tar.gz", emit: rgi_archive, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${tool}_db
    case $tool in
        resfinder)
            git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git ${tool}_db
            TOOL_VERSION=\$(cat ${tool}_db/VERSION 2>/dev/null)
            TOOL_VERSION=\${TOOL_VERSION:-unknown}
            printf '"%s":\n    resfinder: %s\n' "${task.process}" "\$TOOL_VERSION" > versions.yml
            ;;
        
        rgi)
            wget https://card.mcmaster.ca/latest/data -O ${tool}_db.tar.gz
            tar -xvf ${tool}_db.tar.gz -C ${tool}_db
            TOOL_VERSION=\$(python3 -c "import json; d=json.load(open('${tool}_db/card.json')); print(d.get('_version',''))" 2>/dev/null)
            TOOL_VERSION=\${TOOL_VERSION:-unknown}
            printf '"%s":\n    card: %s\n' "${task.process}" "\$TOOL_VERSION" > versions.yml
            ;;
        
        *)
            echo "Unknown tool: $tool"
            exit 1
            ;;
    esac
    """
}