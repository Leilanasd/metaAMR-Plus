# metaAMR-Plus: Output

## Introduction

This document describes the output files produced by metaAMR-Plus. All outputs are saved in the directory specified by `--outdir`.

## Report

### `report/`

| File | Description |
|---|---|
| `metaamr_report.html` | Interactive HTML report integrating all analysis results |

The HTML report is the primary output of the pipeline. It is a single self-contained file that can be opened in any modern browser without an internet connection. It includes:

- **Summary** — AMR heatmap (drug class × sample) and virulence factor heatmap
- **AMR tab** — All resistance genes grouped by drug class, with identity, coverage, tool, and contig columns
- **Virulence tab** — Virulence factors grouped by functional category (Toxins, Iron Acquisition, Immune Evasion, Adherence, etc.)
- **Plasmids tab** — PlasmidFinder replicon hits cross-referenced with PlasClass classifications, contig sizes, and species of origin
- **Taxonomy tab** — Centrifuge and Kaiju results with contig-level assignments and interactive contig search
- **Assembly tab** — Assembly quality statistics (contigs, N50, N90, total length, GC content)

The report automatically adapts based on which tools were run. Sections and columns that require tools not run in a given analysis are hidden automatically.

To generate or regenerate the report manually:

```bash
python3 bin/generate_report.py \
    --results_dir <outdir> \
    --outdir <outdir>/report \
    --run_name "My run"
```

---

## Pre-processing

### `fastqc/`

FastQC quality reports for raw input reads.

### `porechop/`

Adapter-trimmed reads (only saved if `--save_porechop_reads` is set).

### `filtlong/`

Length- and quality-filtered reads (only saved if `--save_filtered_reads` is set).

### `hostremoval/`

| File | Description |
|---|---|
| `*.bam` | Alignment to host genome (only with `--save_hostremoval_bam`) |
| `*_unmapped.fastq.gz` | Reads with host sequences removed (only with `--save_analysis_ready_fastqs`) |

---

## Assembly

### `assemblies/`

Raw Flye assembly FASTA files (`*_assembly.fasta.gz`).

### `polished_assemblies/`

Racon-polished assembly FASTA files (`*_assembly_consensus.fasta.gz`). These are used as input for AMR and plasmid detection  and taxonomy tools.

### `quast/`

QUAST assembly quality assessment per sample. Key metrics include:

- Number of contigs
- Total assembly length
- N50 and N90
- Largest contig
- GC content

---

## AMR detection

### `resfinder/`

ResFinder results per sample:

| File | Description |
|---|---|
| `ResFinder_results_tab.txt` | Tab-delimited AMR gene hits |
| `ResFinder_results.txt` | Detailed results |
| `pheno_table.txt` | Phenotype predictions |

In target species mode, ResFinder results are found in `target_species/amr_results/`.

### `rgi/`

RGI (Resistance Gene Identifier) results using the CARD database:

| File | Description |
|---|---|
| `*.txt` | Main RGI results table |
| `*.json` | JSON format results |

### `amrfinderplus/`

AMRFinderPlus results:

| File | Description |
|---|---|
| `*_amrfinderplus.tsv` | Tab-delimited AMR gene hits |

### `abricate/`

Abricate results using the VFDB database for virulence factor detection:

| File | Description |
|---|---|
| `*.txt` | Tab-delimited virulence factor hits |

### `hamronization/`

hAMRonization integrated results combining RGI, AMRFinderPlus, and Abricate outputs:

| File | Description |
|---|---|
| `summary/hamronization_combined_report.tsv` | Combined AMR results in standardised format |

---

## Plasmid detection

### `plasmidfinder/`

PlasmidFinder replicon typing results:

| File | Description |
|---|---|
| `*_plasmidfinder.tsv` | Detected replicons with identity and contig information |
| `*_plasmidfinder.json` | JSON format results |

### `plasclass/`

PlasClass sequence composition-based plasmid classification:

| File | Description |
|---|---|
| `*.plasclass_classified.txt` | Per-contig classification (plasmid/chromosome/unknown) with probability scores |

---

## Taxonomic profiling

### `centrifuge/`

Centrifuge k-mer based taxonomic classification:

| File | Description |
|---|---|
| `*_centrifuge_report.txt` | Species-level summary report |
| `*_centrifuge_results.txt` | Per-read classification |
| `*_contigs_species.tsv` | Contig-level species assignments (when assembly was run) |

### `kaiju/`

Kaiju protein-level taxonomic classification:

| File | Description |
|---|---|
| `*.txt` | Species summary report |
| `*.tsv` | Per-contig classification with taxon IDs |

### `krona/`

Krona interactive HTML visualisations of taxonomic composition per sample.

---

## Target species mode

### `target_species/`

Output specific to target species mode:

| Directory | Contents |
|---|---|
| `classification/` | Per-sample species detection summaries (`*.species_summary.txt`) |
| `extracted_reads/` | Reads extracted for each target species |
| `amr_results/` | ResFinder results on target species reads |
| `centrifuge/` | Centrifuge classification used for read extraction |

The `*.species_summary.txt` files report detected target species with read counts and confidence levels (High ≥10 reads, Low 1–9 reads).

---

## Quality reporting

### `multiqc/`

| File | Description |
|---|---|
| `multiqc_report.html` | MultiQC interactive quality report |
| `multiqc_data/` | Raw data used by MultiQC |

### `pipeline_info/`

Nextflow execution reports:

| File | Description |
|---|---|
| `execution_timeline_*.html` | Timeline of process execution |
| `execution_report_*.html` | Resource usage report |
| `execution_trace_*.txt` | Detailed process trace |
| `pipeline_dag_*.html` | Pipeline DAG visualisation |
