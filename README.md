# metaAMR-Plus

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.15682600-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.15682600)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

**metaAMR-Plus** is a Nextflow pipeline for comprehensive analysis of long-read Nanopore metagenomic data. It detects antimicrobial resistance (AMR) genes, identifies plasmids and virulence factors, and performs taxonomic classification using multiple parallel tools and databases. Results are integrated into a single interactive HTML report.

The pipeline supports two main modes:

- **Standard mode** — full analysis with assembly, AMR detection, plasmid detection, and taxonomic profiling
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
         │    ├── Assembly (Flye) → [Polishing (Racon) →] QC (QUAST)
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
## Features

- **Multi-tool AMR detection** — [RGI](https://github.com/arpcard/rgi) (CARD), [AMRFinderPlus](https://github.com/ncbi/amr), and [ResFinder](https://github.com/cadms/resfinder) run in parallel; RGI, AMRFinderPlus, and Abricate results are integrated via [hAMRonization](https://github.com/pha4ge/hAMRonization)
- **Virulence factor detection** — [Abricate](https://github.com/tseemann/abricate) (VFDB database); results are included in hAMRonization
- **Plasmid detection** — [PlasmidFinder](https://github.com/cadms/plasmidfinder) (replicon typing) + [PlasClass](https://github.com/Shamir-Lab/PlasClass) (sequence composition)
- **Dual taxonomic profiling** — [Centrifuge](https://github.com/DaehwanKimLab/centrifuge) (k-mer) and [Kaiju](https://github.com/bioinformatics-centre/kaiju) (protein-level) with [Krona](https://github.com/marbl/Krona) visualisation
- **Target species mode** — extract reads classified as specific organisms before AMR analysis
- **Interactive HTML report** — integrated summary with AMR heatmap, VF categorisation, plasmid-taxonomy cross-reference, contig search
- **Adaptive reporting** — report automatically adjusts based on which tools were run
- **Flexible tool selection** — core tools run by default; enable assembly and assembly-based tools (RGI, AMRFinderPlus, PlasmidFinder, PlasClass) explicitly via `--perform_assembly` and `--run_*` flags

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (≥24.04.2)
2. Install [Singularity](https://docs.sylabs.io) or [Docker](https://www.docker.com/)
3. Prepare your samplesheet (see [usage](docs/usage.md))
4. Run the pipeline:

**Minimal run**
```bash
nextflow run Leilanasd/metaAMR-Plus \
    -profile <singularity/docker>
    --outdir results_test
```

**Default run** (FastQC + ResFinder + Centrifuge + Krona + MultiQC + HTML report)
```bash
nextflow run Leilanasd/metaAMR-Plus \
    -profile <singularity/docker> \
    --input samplesheet.csv \
    --outdir results \
    --databases database.csv
```

**Full run** (all tools enabled)
```bash
nextflow run Leilanasd/metaAMR-Plus \
    -profile <singularity/docker> \
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

> If you have existing local databases, provide them in `database.csv` instead of using the download flags.


## Target species mode

To run analysis restricted to specific organisms:

```bash
nextflow run /path/to/metaAMR-Plus \
    -profile <singularity/docker> \
    --input samplesheet.csv \
    --outdir results_kleb \
    --target_species "Klebsiella pneumoniae" \
    --databases database.csv \
    --download_resfinder_db
```

Multiple species can be specified as a comma-separated string:

```bash
--target_species "Klebsiella pneumoniae,Enterococcus faecium"
```


## Samplesheet format

```csv
sample,reads
LH040,/path/to/LH040.fastq.gz
LH024,/path/to/LH024.fastq.gz
LH077,/path/to/LH077.fastq.gz
LH085,/path/to/LH085.fastq.gz
```

> **Important:** The `sample` field must match the base name of the FASTQ file (without `.fastq.gz`) for consistent MultiQC naming.

## Databases

metaAMR-Plus uses several reference databases. These can be downloaded automatically by the pipeline or provided as pre-built local copies.

### AMR databases

| Tool | Database | Auto-download flag | Notes |
|---|---|---|---|
| RGI | CARD | `--download_rgi_db` | Or provide via `database.csv` |
| AMRFinderPlus | NCBI AMRFinder | `--download_amrfinderplus_db` | Or provide via `database.csv` |
| ResFinder | ResFinder DB | `--download_resfinder_db` | Or provide via `database.csv` |

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

Pass it with `--databases database.csv`. No download flags are needed when paths are provided in the database.csv.

### Taxonomy and profiling databases

Centrifuge and Kaiju databases are always provided via the database.csv (see above). These cannot be auto-downloaded and must be built or obtained separately before running the pipeline.

### Built-in databases (no setup required)

| Tool | Database | Notes |
|---|---|---|
| Abricate | VFDB | Bundled with Abricate container |
| PlasmidFinder | PlasmidFinder DB | Bundled with PlasmidFinder container |
| PlasClass | — | Model-based, no external DB needed |



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

metaAMR-Plus was originally written by Leila Nasirzadeh and Jyotirmoy Das.

### Team

- [Leila Nasirzadeh](https://github.com/Leilanasd) — Department of Biomedical and Clinical Sciences, Linköping University
- [Jyotirmoy Das](https://github.com/JD2112) — Department of Biomedical and Clinical Sciences, Linköping University

We thank the following people for their contributions to the development of this pipeline:

- Malgorzata Lysaik
- Jenny Welander
- Haiko Schurz
- Ida Karlsson
- Olivia Andersson
- Samuel Lampa
- Mårten Lindqvist

### Citation

If you use this pipeline, please cite:

> metaAMR-Plus: a Nextflow pipeline for comprehensive AMR analysis of Nanopore metagenomic data.
> DOI: [10.5281/zenodo.15682600](https://doi.org/10.5281/zenodo.15682600)

Please also cite the individual tools used. A full list is available in [CITATIONS.md](CITATIONS.md).

## Acknowledgements
This pipeline was developed at the Bioinformatics Core Facility, Faculty of Medicine and Health Sciences,Linköping University, and Clinical Genomics Linköping, SciLifeLab. 
We are also grateful to LiU-IT (DIGIT) for their support on Fraka (HPC) during development and testing of metaAMR-Plus.


## Contributions and support

If you encounter issues or have suggestions, please open a GitHub issue.
