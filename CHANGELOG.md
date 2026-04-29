# nf-core/metaamrplus: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [2026-04-29]

Initial release of nf-core/metaamrplus, created with the [nf-core](https://nf-co.re/) template.

- Nanopore metagenomics AMR detection pipeline with standard assembly-based and target-species read-based modes
- Trimming and quality filtering with Porechop_ABI and Filtlong
- Host read removal with Minimap2 and SAMtools
- Metagenome assembly with Flye, quality assessment with QUAST, and polishing with Racon
- AMR detection with Abricate, RGI (CARD), AMRFinderPlus, and ResFinder
- Standardised AMR reporting with hAMRonization
- Taxonomic classification and profiling with Centrifuge and Kaiju, visualised with Krona
- Target-species mode: species-specific read filtering followed by ResFinder analysis
- Plasmid detection with PlasmidFinder and PlasClass
- Flexible database handling via CSV sheet, direct params, or auto-download flags

### `Added`

### `Fixed`

### `Dependencies`

| Tool | Version |
|------|---------|
| Porechop_ABI | 0.5.0 |
| Filtlong | 0.2.1 |
| Minimap2 | 2.28 |
| SAMtools | 1.23.1 |
| Flye | 2.9.5 |
| QUAST | 5.3.0 |
| Racon | 1.4.20 |
| Abricate | 1.0.1 |
| RGI | 6.0.5 |
| AMRFinderPlus | 4.2.7 |
| ResFinder | 4.1.11 |
| hAMRonization | 1.1.9 |
| Centrifuge | 1.0.4 |
| Kaiju | 1.10.1 |
| PlasmidFinder | 2.1.6 |
| PlasClass | 0.1.1 |
| FastQC | 0.12.1 |
| MultiQC | 1.25.1 |

### `Deprecated`
