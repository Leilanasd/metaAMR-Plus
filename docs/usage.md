# metaAMR-Plus: Usage

## Introduction

metaAMR-Plus is a Nextflow pipeline for AMR detection, plasmid identification, virulence factor detection, and taxonomic profiling of long-read Nanopore metagenomic data.

Core tools run by default; assembly and assembly-based tools must be explicitly enabled.This allows you to run the exact combination of tools you need without unnecessary computation.

## Samplesheet

Prepare a CSV file with the following columns:

```csv
sample,reads
LH040,/path/to/LH040.fastq.gz
LH085,/path/to/LH085.fastq.gz
LH225,/path/to/LH225.fastq.gz
```

> **Important:** The `sample` name must match the base name of the FASTQ file (without `.fastq.gz`). Mismatched names cause duplicate rows in the MultiQC report.

Pass it to the pipeline with `--input samplesheet.csv`.

## Database 

For taxonomic profiling tools (Centrifuge, Kaiju), existing database is required.
For AMR detection: If you have existing local databases for RGI, AMRFinderPlus, or ResFinder, provide the paths directly in your database.csv file. If not, use the auto-download flags --download_rgi_db, --download_amrfinderplus_db, and/or --download_resfinder_db and the pipeline will download them automatically before running.

```csv

tool,db_name,db_params,db_path
kaiju,kaiju_db,,/path/to/kaiju_db
centrifuge,centrifuge_db,,/path/to/centrifuge_db

```

Pass it with `--databases database.csv`.

## Running the pipeline

### Test run (verify installation)

To verify the pipeline is installed correctly using the built-in test dataset:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile test,singularity \
    --outdir results_test
```

### Default run

The default run executes FastQC + ResFinder + Centrifuge + Krona + MultiQC + HTML report. It requires only your samplesheet and database file:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --databases database.csv \
```


### Full run (all tools enabled)

The following enables all tools that are OFF by default:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --databases database.csv \
    --hostremoval_reference /path/to/host.fa \
    --perform_trim \
    --perform_hostremoval \
    --perform_assembly \
    --perform_polish_assembly \
    --run_abricate \
    --run_plasmidfinder \
    --run_plasclass \
    --run_rgi \
    --run_amrfinderplus \
    --run_kaiju \
    --download_rgi_db \
    --download_amrfinderplus_db
```


### Target species mode

To restrict AMR analysis to reads classified as specific organisms:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --perform_trim \
    --perform_hostremoval \
    --hostremoval_reference /path/to/host_genome.fa \
    --target_species "Klebsiella pneumoniae" \
    --databases database.csv \
    --download_resfinder_db
```

Multiple target species:

```bash
--target_species "Klebsiella pneumoniae,Enterococcus faecium"
```

In target species mode, the pipeline uses Centrifuge to classify reads, extracts reads matching the target species, then runs ResFinder on those reads. Assembly-based tools are not run.

> **Note:** In target species mode, regular ResFinder is automatically skipped regardless of the `--run_resfinder` setting. ResFinder runs automatically on extracted target species reads via the TARGET_SPECIES_AMR subworkflow. No additional flags are needed.

## Parameters

### Input/output

| Parameter | Description | Default |
|---|---|---|
| `--input` | Path to samplesheet CSV | required |
| `--outdir` | Output directory | required |
| `--databases` | Path to database CSV | required for profiling |
| `--run_name` | Optional label shown in the HTML report header | — |

### Pre-processing

| Parameter | Description | Default |
|---|---|---|
| `--perform_trim` | Enable adapter trimming (Porechop_ABI) and quality filtering (Filtlong) | false |
| `--perform_hostremoval` | Enable host read removal (Minimap2) | false |
| `--hostremoval_reference` | Path to host reference genome FASTA | — |
| `--hostremoval_index` | Path to pre-built Minimap2 index (skips indexing step) | — |
| `--filtlong_minlength` | Minimum read length to retain | 1000 |
| `--filtlong_keeppercent` | Percentage of best reads to keep | 80 |
| `--filtlong_targetbases` | Target number of bases to keep | 200000000 |
| `--skip_fastqc` | Skip FastQC quality control | false |
| `--save_filtered_reads` | Save Filtlong-filtered reads to results | false |
| `--save_porechop_reads` | Save Porechop-trimmed reads to results | false |
| `--save_hostremoval_bam` | Save host removal BAM file | false |
| `--save_analysis_ready_fastqs` | Save analysis-ready FASTQs after host removal | false |

### Assembly

| Parameter | Description | Default |
|---|---|---|
| `--perform_assembly` | Enable de novo assembly with Flye | false |
| `--perform_polish_assembly` | Enable assembly polishing with Racon | false |
| `--flye_mode` | Flye assembly mode | `--nano-hq` |

> Use `--nano-hq` for recent ONT data (R10.4+, Q20+). Use `--nano-raw` for older chemistry.

### AMR detection

| Parameter | Description | Default |
|---|---|---|
| `--run_resfinder` | Run ResFinder | true |
| `--resfinder_db` | Path to existing ResFinder database | — |
| `--download_resfinder_db` | Download ResFinder database automatically | false |
| `--resfinder_point_mutations` | Include point mutations in ResFinder analysis | false |
| `--run_rgi` | Run RGI (CARD database) | false |
| `--rgi_db` | Path to existing RGI CARD database | — |
| `--download_rgi_db` | Download RGI database automatically | false |
| `--rgi_min_coverage` | Minimum coverage threshold for RGI | 0.8 |
| `--rgi_min_identity` | Minimum identity threshold for RGI | 0.8 |
| `--run_amrfinderplus` | Run AMRFinderPlus | false |
| `--amrfinderplus_db` | Path to existing AMRFinderPlus database | — |
| `--download_amrfinderplus_db` | Download AMRFinderPlus database automatically | false |
| `--run_abricate` | Run Abricate (virulence factor detection) | false |
| `--abricate_minid` | Minimum identity for Abricate | 80 |
| `--abricate_mincov` | Minimum coverage for Abricate | 80 |
| `--run_hamronization` | Run hAMRonization to integrate AMR results | true |

> **Note:** Abricate runs against the VFDB database for virulence factor detection. Results appear in the Virulence tab of the HTML report, not the AMR tab. hAMRonization integrates results from RGI, AMRFinderPlus, and Abricate.

### Plasmid detection

| Parameter | Description | Default |
|---|---|---|
| `--run_plasmidfinder` | Run PlasmidFinder (replicon typing) | false |
| `--run_plasclass` | Run PlasClass (sequence composition classification) | false |
| `--plasclass_threshold` | PlasClass classification threshold | 0.9 |

### Taxonomic profiling

| Parameter | Description | Default |
|---|---|---|
| `--run_profiling` | Enable taxonomic profiling | true |
| `--run_centrifuge` | Run Centrifuge | true |
| `--run_kaiju` | Run Kaiju | false |
| `--skip_krona` | Skip Krona visualisation | false |
| `--kaiju_taxon_rank` | Taxonomic rank for Kaiju output | species |

### Target species mode

| Parameter | Description | Default |
|---|---|---|
| `--target_species` | Comma-separated list of target species | — |

### Database management

| Parameter | Description | Default |
|---|---|---|
| `--save_databases` | Save downloaded databases to results directory | false |

### HTML report

| Parameter | Description | Default |
|---|---|---|
| `--run_name` | Label shown in the report header | — |

## Resource requirements

Recommended minimum resources per sample:

| Step | CPUs | Memory |
|---|---|---|
| Porechop | 8 | 16 GB |
| Filtlong | 4 | 8 GB |
| Minimap2 (host removal) | 8 | 16 GB |
| Flye (assembly) | 16 | 64 GB |
| Racon (polishing) | 8 | 32 GB |
| RGI | 8 | 32 GB |
| Centrifuge | 8 | 32 GB |
| Kaiju | 8 | 16 GB |

For HPC submission, use `-profile singularity`.

## Resuming runs
Nextflow caches completed processes. To resume an interrupted run:

```bash
nextflow run /path/to/metaAMR-Plus [options] -resume
```

> **Note:** `-resume` reuses cached process results but does **not** automatically regenerate the HTML report. To update the report after a resumed run, regenerate it manually:
>
> ```bash
> python3 bin/generate_report.py \
>     --results_dir <outdir> \
>     --outdir <outdir>/report \
>     --run_name "My run"
> ```



## Running in the background

```bash
nextflow run /path/to/metaAMR-Plus [options] -bg
tail -f .nextflow.log
```
