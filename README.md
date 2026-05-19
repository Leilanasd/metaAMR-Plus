# metaAMR-Plus

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.15682600-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.15682600)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

**metaAMR-Plus** is a Nextflow pipeline for comprehensive analysis of long-read Nanopore metagenomic data. It detects antimicrobial resistance (AMR) genes, identifies plasmids and virulence factors, and performs taxonomic classification using multiple parallel tools and databases. Results are integrated into a single interactive HTML report.

The pipeline supports two main modes:

- **Standard mode** — full analysis with optional assembly, AMR detection, plasmid detection, and taxonomic profiling
- **Target species mode** — extract and analyse reads from specific target organisms (e.g. known pathogens of interest)

## Pipeline overview

```
Input reads (Nanopore FASTQ)
    │
    ├── Quality control (FastQC)
    ├── Adapter trimming (Porechop_ABI)
    ├── Length filtering (Filtlong)
    └── Host removal (Minimap2)
         │
         ├── [WITH ASSEMBLY]
         │    ├── Assembly (Flye) → QC (QUAST) → Polishing (Racon)
         │    ├── AMR detection:        RGI · AMRFinderPlus · ResFinder
         │    ├── Virulence detection:  Abricate (VFDB)
         │    ├── AMR integration:      hAMRonization
         │    ├── Plasmid detection:    PlasmidFinder · PlasClass
         │    └── Taxonomic profiling:  Centrifuge · Kaiju · Krona
         │
         ├── [WITHOUT ASSEMBLY]
         │    ├── ResFinder (reads-based)
         │    └── Taxonomic profiling:  Centrifuge · Kaiju · Krona
         │
         └── [TARGET SPECIES MODE]
              ├── Classification (Centrifuge)
              ├── Read extraction (target species only)
              └── AMR detection (ResFinder)
                   │
                   ▼
         MultiQC + Interactive HTML report
```

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (≥24.04.2)
2. Install [Singularity](https://docs.sylabs.io) or [Docker](https://www.docker.com/)
3. Prepare your samplesheet (see [usage](docs/usage.md))
4. Run the pipeline:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --perform_trim \
    --perform_hostremoval \
    --hostremoval_reference /path/to/host.fa \
    --perform_assembly \
    --perform_polish_assembly \
    --run_rgi --run_amrfinderplus --run_resfinder --run_abricate \
    --run_hamronization \
    --run_profiling --run_centrifuge --run_kaiju \
    --run_plasmidfinder --run_plasclass \
    --databases databases.csv \
    --download_rgi_db \
    --download_resfinder_db \
    --download_amrfinderplus_db
```

> If you have existing local databases, provide them in `databases.csv` instead of using the download flags.

## Features

- **Multi-tool AMR detection** — RGI (CARD), AMRFinderPlus, and ResFinder run in parallel; RGI and AMRFinderPlus results are integrated via hAMRonization
- **Virulence factor detection** — Abricate (VFDB database); results are included in hAMRonization and shown separately in the Virulence tab
- **Plasmid detection** — PlasmidFinder (replicon typing) + PlasClass (sequence composition)
- **Dual taxonomic profiling** — Centrifuge (k-mer) and Kaiju (protein-level) with Krona visualisation
- **Target species mode** — extract reads classified as specific organisms before AMR analysis
- **Interactive HTML report** — integrated summary with AMR heatmap, VF categorisation, plasmid-taxonomy cross-reference, contig search
- **Adaptive reporting** — report automatically adjusts based on which tools were run
- **Fully optional steps** — every tool is opt-in via `--run_*` flags

## Samplesheet format

```csv
sample,reads
LH040,/path/to/LH040.fastq.gz
LH085,/path/to/LH085.fastq.gz
```

> **Important:** The `sample` field must match the base name of the FASTQ file (without `.fastq.gz`) for consistent MultiQC naming.

## Databases

metaAMR-Plus uses several reference databases. These can be downloaded automatically by the pipeline or provided as pre-built local copies.

### AMR databases

| Tool | Database | Auto-download flag | Notes |
|---|---|---|---|
| RGI | CARD | `--download_rgi_db` | Or provide via `databases.csv` |
| AMRFinderPlus | NCBI AMRFinder | `--download_amrfinderplus_db` | Or provide via `databases.csv` |
| ResFinder | ResFinder DB | `--download_resfinder_db` | Or provide via `databases.csv` |

To save downloaded databases for reuse in future runs:

```bash
--save_databases
```

To use existing local databases (skip download), provide them all via the database.csv :

```csv
tool,db_name,db_params,db_path
kaiju,kaiju_db,,/path/to/kaiju_db
centrifuge,centrifuge_db,,/path/to/centrifuge_db
rgi,rgi_db,,/path/to/rgi_db
resfinder,resfinder_db,,/path/to/resfinder_db
amrfinderplus,amrfinderplus_db,,/path/to/amrfinderplus_db
```

Pass it with `--databases databases.csv`. No download flags are needed when paths are provided in the database.csv.

### Taxonomy and profiling databases

Centrifuge and Kaiju databases are always provided via the database.csv (see above). These cannot be auto-downloaded and must be built or obtained separately before running the pipeline.

### Built-in databases (no setup required)

| Tool | Database | Notes |
|---|---|---|
| Abricate | VFDB | Bundled with Abricate container |
| PlasmidFinder | PlasmidFinder DB | Bundled with PlasmidFinder container |
| PlasClass | — | Model-based, no external DB needed |



```csv
tool,db_name,db_params,db_path
kaiju,kaiju_db,,/path/to/kaiju_db
centrifuge,centrifuge_db,,/path/to/centrifuge_db

```

## Target species mode

To run analysis restricted to specific organisms:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results_kleb \
    --target_species "Klebsiella pneumoniae" \
    --databases databases.csv \
    --download_resfinder_db
```

Multiple species can be specified as a comma-separated string:

```bash
--target_species "Klebsiella pneumoniae,Enterococcus faecium"
```

## Output

The main outputs are in `<outdir>/`:

| Directory | Contents |
|---|---|
| `report/` | Interactive HTML report |
| `multiqc/` | MultiQC quality report |
| `assemblies/` | Raw assembly FASTA (Flye) |
| `polished_assemblies/` | Polished assembly FASTA (Racon) |
| `quast/` | Assembly quality statistics |
| `resfinder/` | ResFinder AMR results |
| `rgi/` | RGI (CARD) results |
| `amrfinderplus/` | AMRFinderPlus results |
| `abricate/` | Abricate VF results |
| `hamronization/` | hAMRonization integrated summary |
| `plasmidfinder/` | PlasmidFinder replicon results |
| `plasclass/` | PlasClass contig classification |
| `centrifuge/` | Centrifuge taxonomic profiles |
| `kaiju/` | Kaiju taxonomic profiles |
| `krona/` | Krona interactive visualisations |

See [output documentation](docs/output.md) for full details.

## Credits

metaAMR-Plus was developed at Clinical Genomics Linköping, Sweden. If you use this pipeline, please cite:

> metaAMR-Plus: a Nextflow pipeline for comprehensive AMR analysis of Nanopore metagenomic data.
> DOI: [10.5281/zenodo.15682600](https://doi.org/10.5281/zenodo.15682600)

Please also cite the individual tools used. A full list is available in [CITATIONS.md](CITATIONS.md).

## Contributions and support

If you encounter issues or have suggestions, please open a GitHub issue.
