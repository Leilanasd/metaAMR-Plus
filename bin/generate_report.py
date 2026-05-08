#!/usr/bin/env python3
"""
generate_report.py — metaAMR-Plus HTML Report Generator
========================================================
Generates a single self-contained, offline-capable HTML report for
clinical use (Karolinska). No external dependencies beyond Python stdlib.

Usage:
    python generate_report.py \\
        --results_dir /path/to/results \\
        --outdir      /path/to/report_output \\
        --run_name    "Run-2026-05-07"
"""

import os
import sys
import argparse
import json
import csv
import re
import glob
from datetime import datetime
from collections import defaultdict

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

VERSION = "1.0.0"

# Drug classes flagged as clinically critical (case-insensitive substring match)
CRITICAL_DRUG_CLASSES = {
    "carbapenem", "colistin", "polymyxin", "vancomycin",
    "linezolid", "tigecycline", "methicillin", "fosfomycin",
    "cefiderocol", "ceftazidime", "meropenem", "imipenem",
}

# Canonical display name for each tool
TOOL_LABELS = {
    "rgi":            "RGI",
    "rgi_main":       "RGI",
    "amrfinderplus":  "AMRFinder",
    "amrfinder":      "AMRFinder",
    "abricate":       "Abricate",
    "resfinder":      "ResFinder",
}

STATUS_CRITICAL = "critical"
STATUS_WARNING  = "warning"
STATUS_CLEAN    = "clean"
STATUS_FAILED   = "failed"

# Extensions to strip from input_file_name to recover sample name
STRIP_EXTS = (
    ".fasta", ".fa", ".fna", ".ffn",
    ".fastq.gz", ".fastq", ".fq.gz", ".fq",
    "_polished", "_assembly", "_contigs",
)


# ─────────────────────────────────────────────────────────────────────────────
# Data parsing
# ─────────────────────────────────────────────────────────────────────────────

def _strip_sample_name(raw):
    """Strip path and common suffixes to recover a bare sample name."""
    name = os.path.basename(raw)
    for ext in STRIP_EXTS:
        if name.endswith(ext):
            name = name[: -len(ext)]
    return name


def parse_hamronization(results_dir):
    """
    Parse the hAMRonization combined report.
    Returns: dict[sample_name -> list[amr_entry_dict]]
    """
    filepath = os.path.join(
        results_dir, "hamronization", "summary", "hamronization_combined_report.tsv"
    )
    amr_by_sample = defaultdict(list)

    if not os.path.exists(filepath):
        return amr_by_sample

    with open(filepath, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw    = row.get("input_file_name", "")
            sample = _strip_sample_name(raw)
            if not sample:
                continue

            # drug class: prefer dedicated column, fall back to antimicrobial_agent
            drug_class = (
                row.get("drug_class", "").strip()
                or row.get("antimicrobial_agent", "").strip()
                or "Unknown"
            )

            tool_raw = (row.get("analysis_software_name", "") or "").lower().strip()
            tool     = TOOL_LABELS.get(tool_raw, row.get("analysis_software_name", tool_raw))

            # Skip Abricate VFDB entries — they belong in the VF tab not AMR
            if tool_raw == "abricate":
                continue

            def _pct(key):
                try:
                    return round(float(row.get(key) or 0), 1)
                except (ValueError, TypeError):
                    return 0.0

            amr_by_sample[sample].append({
                "gene":       (row.get("gene_symbol", "") or "").strip(),
                "gene_name":  (row.get("gene_name", "") or "").strip(),
                "drug_class": drug_class,
                "mechanism":  (row.get("resistance_mechanism", "") or "").strip()
                              or _infer_mechanism(row.get("gene_symbol", "") or ""),
                "identity":   _pct("identity_percentage"),
                "coverage":   _pct("coverage_percentage"),
                "tool":       tool,
                "db":         (row.get("reference_database_id", "") or "").strip(),
                "db_version": (row.get("reference_database_version", "") or "").strip(),
                "contig":     (row.get("input_sequence_id", "") or "").strip(),
            })

    return amr_by_sample


def _match_sample(name, known_samples):
    """
    Try to match a derived name to the set of known samples.
    Exact match first; then check if any known sample is a substring.
    """
    if name in known_samples:
        return name
    for s in known_samples:
        if s in name or name in s:
            return s
    return name


def parse_centrifuge(results_dir, sample):
    """
    Parse Centrifuge species report.
    Returns list of top-25 species dicts, sorted by read count desc.
    """
    filepath = os.path.join(
        results_dir, "centrifuge", sample,
        f"{sample}_centrifuge_report.txt"
    )
    if not os.path.exists(filepath):
        filepath = os.path.join(
            results_dir, "target_species", "centrifuge", sample,
            f"{sample}_centrifuge_report.txt"
        )
    taxa = []

    if not os.path.exists(filepath):
        return taxa

    try:
        with open(filepath, encoding="utf-8") as fh:
            lines = [l for l in fh if not l.startswith("#")]
        import io
        with io.StringIO("".join(lines)) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    reads = int(row.get("numReads", 0) or 0)
                except ValueError:
                    reads = 0
                if reads == 0:
                    continue
                taxa.append({
                    "name":      (row.get("name", "") or "").strip(),
                    "rank":      (row.get("taxRank", "") or "").strip(),
                    "reads":     reads,
                    "abundance": 0.0,
                })
        total = sum(t["reads"] for t in taxa)
        for t in taxa:
            t["abundance"] = round(t["reads"] / total * 100, 2) if total > 0 else 0.0
        taxa.sort(key=lambda x: x["reads"], reverse=True)
        return taxa[:25]
    except Exception as e:
        print(f"[generate_report] WARNING: parse_centrifuge failed for sample: {e}")
        return taxa



def parse_kaiju(results_dir, sample):
    """
    Parse Kaiju summary report ({sample}.txt).
    Returns list of top-25 species dicts, sorted by reads desc.
    """
    filepath = os.path.join(results_dir, "kaiju", sample, f"{sample}.txt")
    taxa = []

    if not os.path.exists(filepath):
        return taxa

    try:
        with open(filepath, encoding="utf-8") as fh:
            lines = [l for l in fh if not l.startswith("#")]
        import io
        with io.StringIO("".join(lines)) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = (row.get("taxon_name", "") or "").strip()
                if not name or "unclassified" in name.lower() or "cannot be assigned" in name.lower() or name.strip() == "-":
                    continue
                try:
                    reads = int(float(row.get("reads", 0) or 0))
                except (ValueError, TypeError):
                    reads = 0
                if reads == 0:
                    continue
                try:
                    abundance = round(float(row.get("percent", 0) or 0), 2)
                except (ValueError, TypeError):
                    abundance = 0.0
                taxa.append({
                    "name":      name,
                    "rank":      "species",
                    "reads":     reads,
                    "abundance": abundance,
                    "contigs":   [],
                })
        taxa.sort(key=lambda x: x["reads"], reverse=True)
        return taxa[:25]
    except Exception as e:
        print(f"[generate_report] WARNING: parse_kaiju failed for {sample}: {e}")
        return taxa



def parse_kaiju_contigs(results_dir, sample):
    """
    Build contig→species mapping from Kaiju output.
    Uses {sample}.tsv (contig→taxid) + {sample}.txt (taxid→name).
    Returns dict: {contig_id -> species_name}
    """
    summary_path = os.path.join(results_dir, "kaiju", sample, f"{sample}.txt")
    percontig_path = os.path.join(results_dir, "kaiju", sample, f"{sample}.tsv")
    mapping = {}

    if not os.path.exists(summary_path) or not os.path.exists(percontig_path):
        return mapping

    try:
        # Build taxon_id → name lookup from summary
        taxid_name = {}
        with open(summary_path, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                tid  = (row.get("taxon_id", "") or "").strip()
                name = (row.get("taxon_name", "") or "").strip()
                if tid and name:
                    taxid_name[tid] = name

        # Build contig → species from per-contig file
        with open(percontig_path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                status   = parts[0].strip()
                contig   = parts[1].strip()
                taxid    = parts[2].strip()
                if status == "C" and contig and taxid in taxid_name:
                    species = taxid_name[taxid]
                    if "unclassified" not in species.lower() and                        "cannot be assigned" not in species.lower():
                        mapping[contig] = species
    except Exception as e:
        print(f"[generate_report] WARNING: parse_kaiju_contigs failed for {sample}: {e}")

    return mapping


def parse_quast(results_dir, sample):
    """
    Parse QUAST report.tsv.
    Returns dict of assembly statistics; empty dict if file missing.
    """
    filepath = os.path.join(results_dir, "quast", sample, "report.tsv")
    raw = {}

    if not os.path.exists(filepath):
        return {}

    try:
        with open(filepath, encoding="utf-8") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    raw[parts[0].strip()] = parts[-1].strip()

        def _int(key):
            try:
                return int(raw[key].replace(",", "").replace(" ", ""))
            except (KeyError, ValueError, AttributeError):
                return None

        def _float(key):
            try:
                return float(raw[key])
            except (KeyError, ValueError, AttributeError):
                return None

        return {
            "contigs":      _int("# contigs"),
            "largest":      _int("Largest contig"),
            "total_length": _int("Total length"),
            "n50":          _int("N50"),
            "n90":          _int("N90"),
            "gc":           _float("GC (%)"),
        }
    except Exception as e:
        print(f"[generate_report] WARNING: parse_quast failed for sample: {e}")
        return {}


def parse_plasmidfinder(results_dir, sample):
    """Parse PlasmidFinder results_tab.tsv. Returns list of hit dicts."""
    filepath = os.path.join(
        results_dir, "plasmidfinder", f"{sample}_plasmidfinder.tsv"
    )
    results = []

    if not os.path.exists(filepath):
        return results

    try:
        seen = set()
        with open(filepath, encoding="utf-8") as fh:
            lines = [l for l in fh if not l.startswith("#")]
        import io
        with io.StringIO("".join(lines)) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                plasmid = (row.get("Plasmid", "") or "").strip()
                if not plasmid or "no hit" in plasmid.lower():
                    continue
                contig = (row.get("Contig", "") or "").strip()
                if (plasmid.lower(), contig.lower()) in seen:
                    continue
                seen.add((plasmid.lower(), contig.lower()))
                results.append({
                    "plasmid":   plasmid,
                    "identity":  (row.get("Identity", "") or "").strip(),
                    "database":  (row.get("Database", "") or "").strip(),
                    "contig":    contig,
                    "accession": (row.get("Accession number", "") or "").strip(),
                })
        return results
    except Exception as e:
        print(f"[generate_report] WARNING: parse_plasmidfinder failed for sample: {e}")
        return results


def parse_plasclass(results_dir, sample):
    """
    Parse PlasClass classified output.
    Returns dict with plasmid/chromosome/unknown counts AND per-contig classification.
    """
    filepath = os.path.join(
        results_dir, "plasclass",
        f"{sample}.plasclass_classified.txt"
    )
    result = {"plasmid": 0, "chromosome": 0, "unknown": 0, "contigs": {}}

    if not os.path.exists(filepath):
        return result

    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            contig_id = parts[0].strip()
            cls       = parts[1].strip().lower()
            if contig_id == "Contig_ID" or cls == "classification":
                continue  # skip header
            if cls in ("plasmid", "chromosome"):
                result[cls] += 1
            else:
                result["unknown"] += 1
            result["contigs"][contig_id] = cls

    return result


# Gene prefix → resistance mechanism (for ResFinder which lacks mechanism column)
RESFINDER_MECHANISMS = {
    "bla":   "Antibiotic inactivation (beta-lactamase)",
    "mec":   "Antibiotic target alteration",
    "tet":   "Antibiotic efflux / target protection",
    "msr":   "Antibiotic efflux",
    "mef":   "Antibiotic efflux",
    "erm":   "Antibiotic target methylation",
    "van":   "Antibiotic target alteration",
    "aac":   "Antibiotic inactivation (aminoglycoside)",
    "aph":   "Antibiotic inactivation (aminoglycoside)",
    "ant":   "Antibiotic inactivation (aminoglycoside)",
    "aad":   "Antibiotic inactivation (aminoglycoside)",
    "qnr":   "Antibiotic target protection",
    "cat":   "Antibiotic inactivation (chloramphenicol)",
    "cml":   "Antibiotic efflux",
    "sul":   "Antibiotic target replacement",
    "dfr":   "Antibiotic target replacement",
    "mph":   "Antibiotic inactivation",
    "lnu":   "Antibiotic inactivation",
    "mcr":   "Antibiotic target alteration (colistin)",
    "fos":   "Antibiotic inactivation",
    "cfr":   "Antibiotic target alteration",
    "mdt":   "Efflux pump",
    "oqx":   "Efflux pump",
    "optra": "Antibiotic target protection",
    "poxt":  "Antibiotic target protection",
    "lmr":   "Efflux pump",
}

def _infer_mechanism(gene):
    """Infer resistance mechanism from gene name prefix."""
    g = gene.lower().lstrip("(")
    for prefix, mech in RESFINDER_MECHANISMS.items():
        if g.startswith(prefix):
            return mech
    return "Resistance gene"


def parse_resfinder(results_dir, sample):
    """
    Parse ResFinder_results_tab.txt.
    Returns findings in same format as hAMRonization entries.
    Deduplicates by (gene, drug_class) keeping highest identity hit.
    """
    filepath = os.path.join(
        results_dir, "resfinder", sample, "ResFinder_results_tab.txt"
    )
    results = []
    seen = set()

    if not os.path.exists(filepath):
        # Try target species flat file format
        ts_path = os.path.join(
            results_dir, "target_species", "amr_results",
            f"{sample}.resfinder_results.txt"
        )
        if os.path.exists(ts_path):
            filepath = ts_path
    if not os.path.exists(filepath):
        return results
    json_phenotypes = {}
    json_path = os.path.join(results_dir, "resfinder", sample,
                             "std_format_under_development.json")
    if os.path.exists(json_path):
        try:
            import json as _json
            jdata = _json.load(open(json_path, encoding="utf-8"))
            for ginfo in jdata.get("genes", {}).values():
                name   = ginfo.get("name", "")
                phenos = ginfo.get("phenotypes", [])
                if name and phenos:
                    json_phenotypes.setdefault(name, [p.capitalize() for p in phenos])
        except Exception:
            pass

    # Read file, skipping # comment lines (target species mode adds them)
    with open(filepath, encoding="utf-8") as fh:
        clean_lines = [l for l in fh if not l.startswith("#")]
    import io
    with io.StringIO("".join(clean_lines)) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("Resistance gene", "") or "").strip()
            phenotype = (row.get("Phenotype", "") or "").strip()

            if not gene:
                continue

            # Extract drug class from phenotype string
            # e.g. "Tetracycline resistance" → "Tetracycline"
            # "Macrolide, Lincosamide and Streptogramin B resistance" → keep as-is
            drug_class = phenotype.replace(" resistance", "").strip() or "Unknown"

            # Deduplicate — same gene can appear with multiple accessions
            key = (gene.lower(), drug_class.lower())
            if key in seen:
                continue
            seen.add(key)

            try:
                identity = round(float((row.get("Identity", "") or "0").strip()), 1)
            except ValueError:
                identity = 0.0

            try:
                coverage = round(float((row.get("Coverage", "") or "0").strip()), 1)
            except ValueError:
                coverage = 0.0

            results.append({
                "gene":       gene,
                "gene_name":  ", ".join(json_phenotypes.get(gene, [])[:4]) if json_phenotypes.get(gene) else gene,
                "drug_class": drug_class,
                "mechanism":  _infer_mechanism(gene),
                "identity":   identity,
                "coverage":   coverage,
                "tool":       "ResFinder",
                "db":         "ResFinder",
                "db_version": "",
                "contig":     (row.get("Contig", "") or "").split()[0].strip(),
            })

    return results
# ─────────────────────────────────────────────────────────────────────────────
# Sample discovery & aggregation
# ─────────────────────────────────────────────────────────────────────────────

def discover_samples(results_dir, amr_by_sample):
    """
    Discover all sample names from subdirectory listings + hAMRonization data.
    Returns a sorted list of unique sample names.
    """
    samples = set(amr_by_sample.keys())

    for subdir in [
        "quast", "centrifuge", "plasmidfinder", "plasclass",
        "abricate", "rgi", "amrfinderplus", "resfinder",
    ]:
        path = os.path.join(results_dir, subdir)
        if os.path.isdir(path):
            for entry in os.listdir(path):
                if os.path.isdir(os.path.join(path, entry)):
                    samples.add(entry)
    # Target species mode — flat files
    ts_amr = os.path.join(results_dir, "target_species", "amr_results")
    if os.path.isdir(ts_amr):
        for fname in os.listdir(ts_amr):
            if fname.endswith(".resfinder_results.txt"):
                samples.add(fname.replace(".resfinder_results.txt", ""))
    ts_cent = os.path.join(results_dir, "target_species", "centrifuge")
    if os.path.isdir(ts_cent):
        for entry in os.listdir(ts_cent):
            if os.path.isdir(os.path.join(ts_cent, entry)):
                samples.add(entry)

    # Filter obvious non-sample names
    samples = {s for s in samples if s and not s.startswith(".")}
    return sorted(samples)



def detect_run_config(results_dir):
    """
    Detect which pipeline tools were run based on output directory/file existence.
    Returns dict of booleans used by the HTML report to show/hide sections.
    """
    def has_dir(subdir):
        return os.path.isdir(os.path.join(results_dir, subdir))

    def has_file(*parts):
        return os.path.exists(os.path.join(results_dir, *parts))

    return {
        "assembly":      has_dir("quast"),
        "hamronization": has_file("hamronization", "summary",
                                  "hamronization_combined_report.tsv"),
        "centrifuge":    has_dir("centrifuge"),
        "kaiju":         has_dir("kaiju"),
        "plasmidfinder": has_dir("plasmidfinder"),
        "plasclass":     has_dir("plasclass"),
        "resfinder":     has_dir("resfinder") or has_dir(
                             os.path.join("target_species", "amr_results")),
        "rgi":           has_dir("rgi"),
        "amrfinderplus": has_dir("amrfinderplus"),
        "abricate":      has_dir("abricate"),
        "target_species": has_dir("target_species"),
    }


def get_sample_status(amr_list, assembly, assembly_skipped=False):
    """
    Determine traffic-light status for a sample:
      critical → any critical drug class detected
      warning  → AMR genes found but none critical
      clean    → no AMR genes found
      failed   → no AMR and no assembly data at all
    """
    if not amr_list:
        if assembly_skipped:
            return STATUS_CLEAN
        return STATUS_FAILED if not assembly else STATUS_CLEAN

    for entry in amr_list:
        dc = (entry.get("drug_class", "") or "").lower()
        if any(c in dc for c in CRITICAL_DRUG_CLASSES):
            return STATUS_CRITICAL

    return STATUS_WARNING


def build_summary_matrix(samples, all_amr):
    """
    Build a drug-class × sample count matrix.
    Returns: (sorted_drug_class_list, {sample: {dc: count}})
    Critical classes are sorted to the front.
    """
    all_dcs = set()
    for amr_list in all_amr.values():
        for e in amr_list:
            dc = (e.get("drug_class", "") or "Unknown").strip()
            all_dcs.add(dc)

    def dc_sort_key(dc):
        is_crit = any(c in dc.lower() for c in CRITICAL_DRUG_CLASSES)
        return (0 if is_crit else 1, dc.lower())

    drug_classes = sorted(all_dcs, key=dc_sort_key)

    matrix = {}
    for sample in samples:
        matrix[sample] = {}
        for dc in drug_classes:
            matrix[sample][dc] = sum(
                1 for e in all_amr.get(sample, [])
                if (e.get("drug_class", "") or "").strip() == dc
            )

    return drug_classes, matrix


# ─────────────────────────────────────────────────────────────────────────────
# HTML generation
# ─────────────────────────────────────────────────────────────────────────────

# HTML template — uses __PLACEHOLDER__ markers to avoid f-string brace escaping.
# All JavaScript uses normal {  } syntax.
HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>metaAMR-Plus — __RUN_DATE__</title>
<style>
/* ── Reset ─────────────────────────────────────────────────────────── */
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html,body{height:100%}
body{
  font-family:'Segoe UI',system-ui,-apple-system,sans-serif;
  font-size:13px;line-height:1.5;color:#1a1f2e;
  background:#eef0f4;display:flex;flex-direction:column;overflow:hidden;height:100vh
}

/* ── Header ─────────────────────────────────────────────────────────── */
.hdr{
  background:#0f2744;color:#fff;
  display:flex;align-items:center;gap:24px;
  padding:0 20px;height:52px;flex-shrink:0;
  border-bottom:3px solid #b91c1c;
}
.hdr-logo{display:flex;align-items:center;gap:10px}
.hdr-logo-icon{
  width:30px;height:30px;background:#b91c1c;border-radius:4px;
  display:flex;align-items:center;justify-content:center;font-weight:800;font-size:13px
}
.hdr-title{font-size:15px;font-weight:700;letter-spacing:.2px}
.hdr-subtitle{font-size:10px;opacity:.55;letter-spacing:.3px;margin-top:1px}
.hdr-stats{display:flex;gap:1px;margin-left:auto;background:rgba(255,255,255,.08);border-radius:6px;overflow:hidden}
.hdr-stat{
  padding:6px 16px;text-align:center;
  border-right:1px solid rgba(255,255,255,.08)
}
.hdr-stat:last-child{border-right:none}
.hdr-stat-val{font-size:18px;font-weight:700;line-height:1}
.hdr-stat-val.alarm{color:#fca5a5}
.hdr-stat-lbl{font-size:9px;text-transform:uppercase;letter-spacing:.6px;opacity:.55;margin-top:1px}

/* ── Layout ─────────────────────────────────────────────────────────── */
.layout{display:flex;flex:1;overflow:hidden}

/* ── Sidebar ────────────────────────────────────────────────────────── */
.sidebar{
  width:256px;background:#fff;
  border-right:1px solid #dce1e9;
  display:flex;flex-direction:column;
  flex-shrink:0;overflow:hidden
}

/* section headers inside sidebar */
.sb-hdr{
  padding:7px 12px;font-size:10px;font-weight:700;
  text-transform:uppercase;letter-spacing:.6px;color:#6b7280;
  background:#f9fafb;border-bottom:1px solid #e5e7eb;
  display:flex;justify-content:space-between;align-items:center;
  cursor:pointer;user-select:none;flex-shrink:0
}
.sb-hdr:hover{background:#f3f4f6}
.sb-caret{transition:transform .18s;font-size:9px}
.sb-hdr.collapsed .sb-caret{transform:rotate(-90deg)}

/* summary heatmap */
.sb-summary-wrap{overflow:auto;max-height:220px;flex-shrink:0}
.summary-tbl{border-collapse:collapse;font-size:10.5px;width:max-content;min-width:100%}
.summary-tbl th{
  padding:3px 6px;background:#f9fafb;border:1px solid #e5e7eb;
  font-weight:600;white-space:nowrap;position:sticky;top:0;z-index:1;font-size:10px
}
.summary-tbl th.sc{position:sticky;left:0;z-index:2;background:#f9fafb}
.summary-tbl td{padding:3px 5px;border:1px solid #f0f0f0;text-align:center;cursor:pointer}
.summary-tbl td.sn{
  text-align:left;font-weight:500;position:sticky;left:0;
  background:#fff;z-index:1;padding-left:10px;white-space:nowrap
}
.summary-tbl td.sn:hover{color:#0f2744;text-decoration:underline}
/* cell heat levels */
.c0{background:#fff}
.c1{background:#fef9c3}
.c2{background:#fde68a}
.c3{background:#fbbf24;color:#1a1f2e}
.ch{background:#ef4444;color:#fff;font-weight:700}
.cc{background:#7f1d1d;color:#fecaca;font-weight:700}   /* critical class + hits */

/* sample list */
.sample-list{flex:1;overflow-y:auto}
.si{
  padding:7px 12px 7px 8px;cursor:pointer;
  display:flex;align-items:center;gap:8px;
  border-left:3px solid transparent;
  border-bottom:1px solid #f5f5f5;
  transition:background .1s
}
.si:hover{background:#f8f9fb}
.si.active{background:#eff4ff;border-left-color:#0f2744;font-weight:600}
.dot{width:8px;height:8px;border-radius:50%;flex-shrink:0}
.dot-critical{background:#dc2626}
.dot-warning{background:#f97316}
.dot-clean{background:#16a34a}
.dot-failed{background:#9ca3af}
.si-count{margin-left:auto;font-size:10px;color:#9ca3af}
.si-badge{
  margin-left:2px;font-size:9px;padding:1px 5px;border-radius:8px;
  background:#fee2e2;color:#991b1b;font-weight:700
}

/* ── Content ────────────────────────────────────────────────────────── */
.content{flex:1;display:flex;flex-direction:column;overflow:hidden;min-width:0}

.content-hdr{
  background:#fff;border-bottom:1px solid #e5e7eb;
  padding:10px 20px 0;flex-shrink:0
}
.content-sample-name{font-size:15px;font-weight:700;color:#0f2744;margin-bottom:7px}
.content-sample-name span{font-weight:400;font-size:12px;color:#6b7280;margin-left:8px}

/* tabs */
.tabs{display:flex;gap:0}
.tab{
  padding:6px 16px;cursor:pointer;border-bottom:2px solid transparent;
  margin-bottom:-1px;font-size:12.5px;color:#6b7280;font-weight:500;
  transition:color .12s;white-space:nowrap
}
.tab:hover{color:#0f2744}
.tab.active{color:#0f2744;border-bottom-color:#0f2744;font-weight:700}
.tb{
  display:inline-block;margin-left:4px;padding:1px 6px;border-radius:9px;
  font-size:10px;background:#e5e7eb;color:#374151;font-weight:600
}
.tb.alarm{background:#fca5a5;color:#7f1d1d}

/* content body */
.cbody{flex:1;overflow-y:auto;padding:16px 20px}

/* ── Toolbar ────────────────────────────────────────────────────────── */
.toolbar{display:flex;gap:8px;align-items:center;margin-bottom:14px;flex-wrap:wrap}
.toolbar select,.toolbar input[type=text]{
  padding:5px 9px;border:1px solid #d1d5db;border-radius:5px;
  font-size:12px;background:#fff;color:#1a1f2e;outline:none
}
.toolbar select:focus,.toolbar input:focus{border-color:#0f2744}
.toolbar input[type=text]{min-width:200px}
.tlbl{font-size:11px;font-weight:600;color:#6b7280}

/* ── AMR groups ─────────────────────────────────────────────────────── */
.dc-group{margin-bottom:14px;border:1px solid #e5e7eb;border-radius:6px;overflow:hidden}
.dc-hdr{
  display:flex;align-items:center;gap:8px;padding:7px 12px;
  background:#f3f6fb;font-weight:600;font-size:12px;cursor:pointer;
  user-select:none;border-bottom:1px solid #e5e7eb
}
.dc-hdr.crit{background:#fff0f0;border-bottom-color:#fecaca}
.dc-hdr-icon{font-size:13px}
.dc-cnt{
  margin-left:auto;font-size:10px;padding:1px 7px;
  border-radius:9px;background:#0f2744;color:#fff;font-weight:700
}
.dc-cnt.crit{background:#dc2626}
.dc-body{display:block}

/* AMR table */
.amr-tbl{width:100%;border-collapse:collapse;font-size:12px}
.amr-tbl th{
  padding:5px 10px;background:#f9fafb;border-bottom:1px solid #e5e7eb;
  text-align:left;font-weight:600;font-size:10.5px;white-space:nowrap;color:#4b5563
}
.amr-tbl td{padding:6px 10px;border-bottom:1px solid #f3f4f6;vertical-align:top}
.amr-tbl tr:last-child td{border-bottom:none}
.amr-tbl tr:hover td{background:#f8faff}
.gene-name{font-weight:700;font-size:12.5px}
.gene-full{font-size:10.5px;color:#6b7280;margin-top:1px}

/* identity coloring */
.id-hi{color:#15803d;font-weight:700}
.id-md{color:#b45309;font-weight:600}
.id-lo{color:#9ca3af}

/* tool badges */
.tbadge{
  display:inline-block;padding:2px 7px;border-radius:4px;
  font-size:10px;font-weight:700;letter-spacing:.2px;white-space:nowrap
}
.tbadge-rgi     {background:#d1fae5;color:#065f46}
.tbadge-amr     {background:#dbeafe;color:#1e3a8a}
.tbadge-abr     {background:#fce7f3;color:#9d174d}
.tbadge-res     {background:#fef3c7;color:#92400e}
.tbadge-other   {background:#f3f4f6;color:#374151}

/* ── Taxonomy ────────────────────────────────────────────────────────── */
.tax-tbl{width:100%;border-collapse:collapse;font-size:12px}
.tax-tbl th{padding:5px 10px;background:#f9fafb;border-bottom:1px solid #e5e7eb;
  font-weight:600;font-size:10.5px;color:#4b5563}
.tax-tbl td{padding:6px 10px;border-bottom:1px solid #f3f4f6}
.tax-tbl tr:hover td{background:#f8faff}
.tax-rank{font-size:10px;color:#9ca3af}
.bar-cell{min-width:120px}
.bar-wrap{background:#e5e7eb;border-radius:3px;height:8px;overflow:hidden}
.bar-fill{background:#0f2744;height:8px;transition:width .3s}

/* ── Assembly stats ──────────────────────────────────────────────────── */
.stat-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(155px,1fr));gap:10px;margin-bottom:18px}
.stat-card{
  background:#f9fafb;border:1px solid #e5e7eb;border-radius:6px;
  padding:13px 15px;
}
.stat-lbl{font-size:10px;text-transform:uppercase;letter-spacing:.5px;color:#6b7280;font-weight:700}
.stat-val{font-size:22px;font-weight:800;color:#0f2744;margin-top:3px;line-height:1}
.stat-unit{font-size:11px;color:#9ca3af;font-weight:400}

/* ── Plasmids ────────────────────────────────────────────────────────── */
.pc-grid{display:flex;gap:10px;margin-bottom:14px;flex-wrap:wrap}
.pc-card{
  background:#f9fafb;border:1px solid #e5e7eb;border-radius:6px;
  padding:10px 18px;text-align:center;min-width:110px
}
.section-title{
  font-size:10.5px;font-weight:700;text-transform:uppercase;letter-spacing:.5px;
  color:#6b7280;margin-bottom:8px;padding-bottom:6px;border-bottom:1px solid #e5e7eb
}
.pf-tbl{width:100%;border-collapse:collapse;font-size:12px}
.pf-tbl th{padding:5px 10px;background:#f9fafb;border-bottom:1px solid #e5e7eb;
  font-weight:600;font-size:10.5px;color:#4b5563}
.pf-tbl td{padding:6px 10px;border-bottom:1px solid #f3f4f6}

/* ── Empty state ─────────────────────────────────────────────────────── */
.empty{text-align:center;padding:48px 20px;color:#9ca3af}
.empty svg{margin-bottom:10px;opacity:.35}
.empty-msg{font-size:13px;line-height:1.6}

/* ── Welcome ─────────────────────────────────────────────────────────── */
.welcome{
  display:flex;flex-direction:column;align-items:center;justify-content:center;
  height:100%;color:#9ca3af;text-align:center;gap:10px
}
.welcome h2{font-size:17px;color:#6b7280;font-weight:600}
.welcome p{font-size:12.5px;max-width:340px;line-height:1.6}

/* ── Print ───────────────────────────────────────────────────────────── */
@media print{
  .sidebar,.toolbar{display:none!important}
  .layout{display:block}
  body,.content{overflow:visible;height:auto}
  .cbody{overflow:visible;padding:0}
  .tab-content{display:block!important}
  .tabs,.hdr{display:none}
  .dc-body{display:block!important}
  .amr-tbl tr:hover td{background:none}
}
</style>
</head>
<body>

<!-- ─── Header ──────────────────────────────────────────────────────── -->
<header class="hdr">
  <div class="hdr-logo">
    <div class="hdr-logo-icon">AMR</div>
    <div>
      <div class="hdr-title">metaAMR-Plus</div>
      <div class="hdr-subtitle">CLINICAL REPORT · __RUN_DATE__</div>
    </div>
  </div>
  <div id="hdr-tools" style="display:flex;gap:6px;font-size:10px;align-items:center;
    margin-left:16px;flex-wrap:wrap"></div>
  <div class="hdr-stats">
    <div class="hdr-stat">
      <div class="hdr-stat-val" id="hdr-samples">__N_SAMPLES__</div>
      <div class="hdr-stat-lbl">Samples</div>
    </div>
    <div class="hdr-stat">
      <div class="hdr-stat-val" id="hdr-amr">__N_AMR__</div>
      <div class="hdr-stat-lbl">AMR Genes</div>
    </div>
    <div class="hdr-stat">
      <div class="hdr-stat-val __CRITICAL_CLASS__" id="hdr-critical">__N_CRITICAL__</div>
      <div class="hdr-stat-lbl">Critical</div>
    </div>
  </div>
</header>

<div class="layout">

  <!-- ─── Sidebar ─────────────────────────────────────────────────── -->
  <aside class="sidebar">

    <!-- Summary heatmap -->
    <div id="summary-section">
      <div class="sb-hdr" id="summary-hdr" onclick="toggleSummary()">
        Summary <span class="sb-caret" id="summary-caret">▾</span>
      </div>
      <div class="sb-summary-wrap" id="summary-body">
        <table class="summary-tbl" id="summary-table"></table>
      </div>
    </div>

    <!-- Sample list -->
    <div class="sb-hdr" style="cursor:default">Samples</div>
    <div class="sample-list" id="sample-list"></div>

  </aside>

  <!-- ─── Content ─────────────────────────────────────────────────── -->
  <main class="content">

    <div class="content-hdr" id="content-hdr" style="display:none">
      <div class="content-sample-name" id="content-title"></div>
      <div class="tabs">
        <div class="tab active" id="tab-btn-amr"      onclick="switchTab('amr')">AMR <span class="tb" id="tb-amr">0</span></div>
        <div class="tab"        id="tab-btn-vf"        onclick="switchTab('vf')">Virulence <span class="tb" id="tb-vf">0</span></div>
        <div class="tab"        id="tab-btn-plasmids"  onclick="switchTab('plasmids')">Plasmids <span class="tb" id="tb-plasmids">0</span></div>
        <div class="tab"        id="tab-btn-taxonomy"  onclick="switchTab('taxonomy')">Taxonomy <span class="tb" id="tb-taxonomy">0</span></div>
        <div class="tab"        id="tab-btn-assembly"  onclick="switchTab('assembly')">Assembly <span class="tb" id="tb-assembly">—</span></div>
      </div>
    </div>

    <div class="cbody" id="cbody">
      <div class="welcome" id="welcome-panel">
        <h2>Select a sample</h2>
        <p>Click any sample in the sidebar to view its AMR findings, taxonomy, plasmid detection, and assembly statistics.</p>
      </div>
      <div id="tab-amr"      class="tab-content" style="display:none"></div>
      <div id="tab-vf"       class="tab-content" style="display:none"></div>
      <div id="tab-plasmids" class="tab-content" style="display:none"></div>
      <div id="tab-taxonomy" class="tab-content" style="display:none"></div>
      <div id="tab-assembly" class="tab-content" style="display:none"></div>
    </div>

  </main>
</div>

<script>
// ── Embedded pipeline data ───────────────────────────────────────────
const DATA = __DATA_JSON__;
const CRITICAL_DRUG_CLASSES = __CRITICAL_CLASSES_JSON__;
const CONFIG = DATA.config || {};

// ── App state ────────────────────────────────────────────────────────
let currentSample = null;
let currentTab    = 'amr';
let filterDC      = '';
let filterQuery   = '';

// ── Utility functions ────────────────────────────────────────────────
function isCrit(dc) {
  const s = (dc || '').toLowerCase();
  return CRITICAL_DRUG_CLASSES.some(c => s.includes(c));
}

function esc(str) {
  return String(str ?? '')
    .replace(/&/g,'&amp;')
    .replace(/</g,'&lt;')
    .replace(/>/g,'&gt;')
    .replace(/"/g,'&quot;');
}

function fmtN(n) {
  if (n == null) return '—';
  return Number(n).toLocaleString();
}

function idClass(pct) {
  if (pct >= 95) return 'id-hi';
  if (pct >= 80) return 'id-md';
  return 'id-lo';
}

function toolBadgeClass(tool) {
  const t = (tool || '').toLowerCase();
  if (t.includes('rgi'))         return 'tbadge-rgi';
  if (t.includes('amrfinder'))   return 'tbadge-amr';
  if (t.includes('abricate'))    return 'tbadge-abr';
  if (t.includes('resfinder'))   return 'tbadge-res';
  return 'tbadge-other';
}

function emptyState(msg) {
  return `<div class="empty">
    <svg width="44" height="44" viewBox="0 0 24 24" fill="none"
         stroke="currentColor" stroke-width="1.5" stroke-linecap="round">
      <circle cx="12" cy="12" r="10"/>
      <line x1="12" y1="8" x2="12" y2="12"/>
      <line x1="12" y1="16" x2="12.01" y2="16"/>
    </svg>
    <div class="empty-msg">${msg}</div>
  </div>`;
}

// ── Summary table ────────────────────────────────────────────────────
function buildSummaryTable() {
  const tbl = document.getElementById('summary-table');
  const dcs = DATA.drug_classes;
  if (!dcs || dcs.length === 0) {
    tbl.innerHTML = '<tr><td style="padding:8px;color:#9ca3af;font-size:11px">No AMR data.</td></tr>';
    return;
  }

  // Header row with rotated drug class labels
  let html = '<thead><tr><th class="sc">Sample</th>';
  dcs.forEach(dc => {
    const crit = isCrit(dc);
    const label = dc.length > 16 ? dc.slice(0,14) + '…' : dc;
    html += `<th title="${esc(dc)}" style="padding:3px 4px;` +
      (crit ? 'color:#dc2626;' : '') + '">' +
      `<div style="writing-mode:vertical-lr;transform:rotate(180deg);` +
      `white-space:nowrap;max-height:80px;font-size:9.5px">${esc(label)}</div></th>`;
  });
  html += '</tr></thead><tbody>';

  DATA.samples.forEach(s => {
    html += `<tr><td class="sn" onclick="selectSample('${esc(s)}')">${esc(s)}</td>`;
    dcs.forEach(dc => {
      const cnt  = (DATA.summary_matrix[s] || {})[dc] || 0;
      const crit = isCrit(dc) && cnt > 0;
      let cls = 'c0';
      if (crit)       cls = 'cc';
      else if (cnt >= 6) cls = 'ch';
      else if (cnt === 3) cls = 'c3';
      else if (cnt === 2) cls = 'c2';
      else if (cnt === 1) cls = 'c1';
      html += `<td class="${cls}"
        onclick="selectSampleDC('${esc(s)}','${esc(dc)}')"
        title="${esc(s)} — ${esc(dc)}: ${cnt} gene(s)">${cnt || ''}</td>`;
    });
    html += '</tr>';
  });

  html += '</tbody>';
  tbl.innerHTML = html;
}

// ── Sample list ──────────────────────────────────────────────────────
function buildSampleList() {
  const list = document.getElementById('sample-list');
  let html = '';
  DATA.samples.forEach(s => {
    const status   = DATA.status[s] || 'clean';
    const amrCount = (DATA.amr[s] || []).length;
    const hasCrit  = (DATA.amr[s] || []).some(e => isCrit(e.drug_class));
    html += `
    <div class="si" id="si-${esc(s)}" onclick="selectSample('${esc(s)}')">
      <div class="dot dot-${esc(status)}"></div>
      <span>${esc(s)}</span>
      ${amrCount > 0 ? `<span class="si-count">${amrCount}</span>` : ''}
      ${hasCrit ? `<span class="si-badge">!</span>` : ''}
    </div>`;
  });
  list.innerHTML = html;
}

// ── Select sample ────────────────────────────────────────────────────
function selectSample(s) {
  if (currentSample) {
    const prev = document.getElementById('si-' + currentSample);
    if (prev) prev.classList.remove('active');
  }
  currentSample = s;
  filterDC    = '';
  filterQuery = '';

  const el = document.getElementById('si-' + s);
  if (el) { el.classList.add('active'); el.scrollIntoView({block:'nearest'}); }

  document.getElementById('welcome-panel').style.display = 'none';
  document.getElementById('content-hdr').style.display   = '';

  const status = DATA.status[s] || 'clean';
  const statusLabels = {critical:'⚠ Critical resistance', warning:'Resistance detected',
                        clean:'No resistance detected', failed:'Assembly/data missing'};
  document.getElementById('content-title').innerHTML =
    `${esc(s)} <span>${esc(statusLabels[status] || '')}</span>`;

  updateBadges();
  renderTab(currentTab);
}

function selectSampleDC(s, dc) {
  selectSample(s);
  switchTab('amr');
  filterDC = dc;
  filterQuery = '';
  renderAmr();
}

// ── Tabs ─────────────────────────────────────────────────────────────
function switchTab(tab) {
  currentTab = tab;
  ['amr','vf','plasmids','taxonomy','assembly'].forEach(t => {
    document.getElementById('tab-' + t).style.display        = t === tab ? '' : 'none';
    document.getElementById('tab-btn-' + t).classList.toggle('active', t === tab);
  });
  renderTab(tab);
}

function renderTab(tab) {
  if (!currentSample) return;
  switch(tab) {
    case 'amr':      renderAmr();      break;
    case 'vf':       renderVF();       break;
    case 'plasmids': renderPlasmids(); break;
    case 'taxonomy': renderTaxonomy(); break;
    case 'assembly': renderAssembly(); break;
  }
}

function updateBadges() {
  if (!currentSample) return;
  const amr   = DATA.amr[currentSample]      || [];
  const vfall = DATA.vf[currentSample]       || [];
  const taxa  = DATA.taxonomy[currentSample] || [];
  const pf    = (DATA.plasmids[currentSample] || {}).plasmidfinder || [];
  const asm   = DATA.assembly[currentSample]  || {};

  const tbAMR = document.getElementById('tb-amr');
  tbAMR.textContent = amr.length;
  tbAMR.className   = 'tb' + (amr.some(e => isCrit(e.drug_class)) ? ' alarm' : '');

  document.getElementById('tb-vf').textContent       = vfall.length || '0';
  document.getElementById('tb-plasmids').textContent = pf.length || '0';
  document.getElementById('tb-taxonomy').textContent = taxa.length || '0';
  document.getElementById('tb-assembly').textContent =
    asm.n50 ? fmtN(asm.n50) + ' N50' : (Object.keys(asm).length ? 'OK' : '—');
}

// ── AMR tab ──────────────────────────────────────────────────────────
function renderAmr() {
  const el     = document.getElementById('tab-amr');
  const allAmr = DATA.amr[currentSample] || [];

  // Build context banners
  let banners = '';
  if (!CONFIG.hamronization) {
    const tsMsg = CONFIG.target_species
      ? 'Target species mode — showing ResFinder results for targeted species only. Assembly-based tools (RGI, AMRFinderPlus, Abricate) were not run.'
      : CONFIG.resfinder
        ? 'hAMRonization was not run. Showing ResFinder results only. Enable <code>--run_hamronization true</code> for full AMR detection.'
        : 'No AMR data available. Run with <code>--run_hamronization true</code> and/or <code>--run_resfinder true</code>.';
    banners += `<div style="background:#fef3c7;border-left:3px solid #f59e0b;
      padding:9px 14px;margin-bottom:12px;border-radius:3px;font-size:12px">
      ⓘ ${tsMsg}
    </div>`;
  }
  if (!CONFIG.assembly && CONFIG.resfinder && !CONFIG.target_species) {
    banners += `<div style="background:#eff6ff;border-left:3px solid #3b82f6;
      padding:9px 14px;margin-bottom:12px;border-radius:3px;font-size:12px">
      ⓘ Assembly was skipped — showing ResFinder results only.
      Assembly-based tools (RGI, AMRFinderPlus, Abricate) were not run.
    </div>`;
  }

  if (allAmr.length === 0) {
    el.innerHTML = banners + emptyState('No AMR findings for this sample.');
    return;
  }

  // Drug class options for select
  const uniqDC = [...new Set(allAmr.map(e => e.drug_class || 'Unknown'))].sort((a,b) => {
    const ca = isCrit(a), cb = isCrit(b);
    return ca && !cb ? -1 : !ca && cb ? 1 : a.localeCompare(b);
  });

  let dcOptions = '<option value="">All drug classes</option>';
  uniqDC.forEach(dc => {
    dcOptions += `<option value="${esc(dc)}"${filterDC===dc?' selected':''}>${esc(dc)}</option>`;
  });

  let html = `<div class="toolbar">
    <span class="tlbl">Class:</span>
    <select id="dc-sel" onchange="filterDC=this.value;filterQuery='';renderAmr()">
      ${dcOptions}
    </select>
    <span class="tlbl">Gene:</span>
    <input type="text" id="gene-search" placeholder="Search gene symbol…"
      value="${esc(filterQuery)}"
      oninput="filterQuery=this.value;renderAmr()">
  </div>`;

  // Apply filters
  const q = (filterQuery || '').toLowerCase();
  const filtered = allAmr.filter(e => {
    if (filterDC && e.drug_class !== filterDC) return false;
    if (q && !(
      (e.gene || '').toLowerCase().includes(q) ||
      (e.gene_name || '').toLowerCase().includes(q)
    )) return false;
    return true;
  });

  if (filtered.length === 0) {
    el.innerHTML = html + emptyState('No findings match the current filter.');
    return;
  }

  // Group by drug class, sort critical first
  const groups = {};
  filtered.forEach(e => {
    const dc = e.drug_class || 'Unknown';
    (groups[dc] = groups[dc] || []).push(e);
  });
  const sortedDCs = Object.keys(groups).sort((a,b) => {
    const ca = isCrit(a), cb = isCrit(b);
    return ca && !cb ? -1 : !ca && cb ? 1 : a.localeCompare(b);
  });

  sortedDCs.forEach(dc => {
    const entries = groups[dc].sort((a,b) => b.identity - a.identity);
    const crit    = isCrit(dc);
    html += `<div class="dc-group">
      <div class="dc-hdr${crit?' crit':''}"
           onclick="this.nextElementSibling.style.display=this.nextElementSibling.style.display==='none'?'':'none'">
        <span class="dc-hdr-icon">${crit ? '⚠' : '▸'}</span>
        <span>${esc(dc)}</span>
        <span class="dc-cnt${crit?' crit':''}">${entries.length}</span>
      </div>
      <div class="dc-body">
        <table class="amr-tbl">
          <thead><tr>
            <th>Gene</th>
            <th>Resistance mechanism</th>
            <th style="width:80px">Identity</th>
            <th style="width:80px">Coverage</th>
            <th style="width:90px">Tool</th>
            ${CONFIG.centrifuge ? '<th style="width:110px">Contig</th>' : ''}
          </tr></thead>
          <tbody>`;

    entries.forEach(e => {
      html += `<tr>
        <td>
          <div class="gene-name">${esc(e.gene || '—')}</div>
          ${e.gene_name ? `<div class="gene-full">${esc(e.gene_name)}</div>` : ''}
        </td>
        <td style="max-width:230px;font-size:11.5px">${esc(e.mechanism || '—')}</td>
        <td class="${idClass(e.identity)}">${e.identity ? e.identity.toFixed(1)+'%' : '—'}</td>
        <td style="color:#6b7280">${e.coverage ? e.coverage.toFixed(1)+'%' : '—'}</td>
        <td><span class="tbadge ${toolBadgeClass(e.tool)}">${esc(e.tool || '—')}</span></td>
        ${CONFIG.centrifuge ? `<td style="font-size:10.5px;font-family:monospace;color:#6b7280">${esc(e.contig || '—')}</td>` : ''}
      </tr>`;
    });

    html += '</tbody></table></div></div>';
  });

  el.innerHTML = html;

  // Restore cursor in search field
  if (q) {
    const inp = document.getElementById('gene-search');
    if (inp) { inp.focus(); inp.setSelectionRange(inp.value.length, inp.value.length); }
  }
}

// ── Taxonomy tab ─────────────────────────────────────────────────────
let taxSource = 'centrifuge';  // 'centrifuge' or 'kaiju'

function switchTaxSource(src) {
  taxSource = src;
  renderTaxonomy();
}

function updateToggleButtons() {
  ['centrifuge','kaiju'].forEach(s => {
    const btn = document.getElementById('tax-btn-' + s);
    if (!btn) return;
    const active = s === taxSource;
    btn.style.background = active ? '#0f2744' : '#fff';
    btn.style.color      = active ? '#fff'    : '#374151';
    btn.style.fontWeight = active ? '700'     : '500';
    btn.style.border     = '2px solid ' + (active ? '#0f2744' : '#d1d5db');
  });
  const note = document.getElementById('tax-src-note');
  if (note) note.textContent = CONFIG.assembly
    ? 'ⓘ Contig-level classification — contig search available'
    : 'ⓘ Read-level classification — no assembly was run';
}

function renderTaxonomy() {
  const el = document.getElementById('tab-taxonomy');

  // Determine which source to show
  const useCentrifuge = taxSource === 'centrifuge' || !CONFIG.kaiju;
  const taxa = useCentrifuge
    ? (DATA.taxonomy[currentSample] || [])
    : (DATA.kaiju[currentSample]    || []);

  // Source toggle
  let toggleHtml = '';
  if (CONFIG.centrifuge && CONFIG.kaiju) {
    toggleHtml = `<div style="background:#f8f9fa;border:1px solid #e5e7eb;
      border-radius:6px;padding:10px 14px;margin-bottom:14px">
      <div style="font-size:11px;font-weight:700;color:#6b7280;
        text-transform:uppercase;letter-spacing:.5px;margin-bottom:8px">
        Classification source
      </div>
      <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap">
        <button id="tax-btn-centrifuge"
          onclick="switchTaxSource('centrifuge')"
          style="padding:5px 14px;border-radius:5px;cursor:pointer;font-size:12px;
                 font-weight:${taxSource==='centrifuge'?'700':'500'};
                 border:2px solid ${taxSource==='centrifuge'?'#0f2744':'#d1d5db'};
                 background:${taxSource==='centrifuge'?'#0f2744':'#fff'};
                 color:${taxSource==='centrifuge'?'#fff':'#374151'}">
          ● Centrifuge
        </button>
        <button id="tax-btn-kaiju"
          onclick="switchTaxSource('kaiju')"
          style="padding:5px 14px;border-radius:5px;cursor:pointer;font-size:12px;
                 font-weight:${taxSource==='kaiju'?'700':'500'};
                 border:2px solid ${taxSource==='kaiju'?'#0f2744':'#d1d5db'};
                 background:${taxSource==='kaiju'?'#0f2744':'#fff'};
                 color:${taxSource==='kaiju'?'#fff':'#374151'}">
          ● Kaiju
        </button>
        <span style="font-size:11px;color:#6b7280;margin-left:4px">
          ${CONFIG.assembly ? (useCentrifuge ? 'ⓘ Contig-level — contig search available' : 'ⓘ Contig-level') : 'ⓘ Read-level — no assembly'}


        </span>
      </div>
    </div>`;
  }

  // Target species classification — show prominently if available
  const targetClass = (DATA.target_classification || {})[currentSample] || [];
  let targetHtml = '';
  if (CONFIG.target_species && targetClass.length > 0) {
    targetHtml = '<div style="margin-bottom:16px;border:1px solid #e5e7eb;border-radius:6px;overflow:hidden">'
      + '<div style="background:#0f2744;color:#fff;padding:8px 14px;font-size:12px;font-weight:700">'
      + 'Target Species Classification</div>'
      + '<table style="width:100%;border-collapse:collapse;font-size:12px">'
      + '<thead><tr>'
      + '<th style="padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e5e7eb;text-align:left">Species</th>'
      + '<th style="padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e5e7eb;text-align:right">Reads</th>'
      + '<th style="padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e5e7eb">Status</th>'
      + '<th style="padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e5e7eb">Confidence</th>'
      + '</tr></thead><tbody>';
    targetClass.forEach(t => {
      const highConf = t.confidence === 'High';
      targetHtml += '<tr>'
        + '<td style="padding:7px 10px;border-bottom:1px solid #f3f4f6"><em>' + esc(t.species) + '</em></td>'
        + '<td style="padding:7px 10px;border-bottom:1px solid #f3f4f6;text-align:right">' + fmtN(t.count) + '</td>'
        + '<td style="padding:7px 10px;border-bottom:1px solid #f3f4f6">'
        + '<span style="background:#d1fae5;color:#065f46;padding:2px 8px;border-radius:4px;font-size:10px;font-weight:700">'
        + esc(t.status) + '</span></td>'
        + '<td style="padding:7px 10px;border-bottom:1px solid #f3f4f6">'
        + '<span style="background:' + (highConf ? '#d1fae5;color:#065f46' : '#fef3c7;color:#92400e')
        + ';padding:2px 8px;border-radius:4px;font-size:10px;font-weight:700">'
        + esc(t.confidence) + '</span></td>'
        + '</tr>';
    });
    targetHtml += '</tbody></table>'
      + '<div style="padding:6px 14px;font-size:10.5px;color:#9ca3af;background:#f9fafb">'
      + 'ⓘ High confidence: ≥10 reads. Low confidence: 1–9 reads (possible misclassification).'
      + '</div></div>';
  }

  if (taxa.length === 0 && targetClass.length === 0) {
    el.innerHTML = toggleHtml + emptyState('No taxonomy data available.<br>Profiling may have been skipped.');
    return;
  }

  // Build full contig→species lookup — Centrifuge or Kaiju
  const fullContigSp = useCentrifuge
    ? (DATA.contig_species || {})[currentSample] || {}
    : (DATA.kaiju_contig_species || {})[currentSample] || {};

  const maxReads = Math.max(...taxa.map(t => t.reads), 1);
  let html = `<table class="tax-tbl">
    <thead><tr>
      <th style="width:28px">#</th>
      <th>Species</th>
      <th style="width:60px">Rank</th>
      <th style="width:70px">Reads</th>
      <th style="width:70px">Abundance</th>
      <th class="bar-cell">Relative abundance</th>
      ${(CONFIG.centrifuge || CONFIG.kaiju) && CONFIG.assembly ? '<th>Contigs</th>' : ''}
    </tr></thead><tbody>`;

  taxa.forEach((t, i) => {
    const barPct = Math.round(t.reads / maxReads * 100);
    const contigs = t.contigs || [];
    const contigDisplay = contigs.length > 0
      ? (contigs.slice(0,3).map(c => esc(c)).join(', ')
         + (contigs.length > 3 ? ` <span style="color:#9ca3af">+${contigs.length-3} more</span>` : ''))
      : '—';
    html += `<tr class="tax-row">
      <td style="color:#9ca3af;text-align:right">${i+1}</td>
      <td><em>${esc(t.name)}</em></td>
      <td class="tax-rank">${esc(t.rank)}</td>
      <td style="text-align:right">${fmtN(t.reads)}</td>
      <td style="text-align:right">${t.abundance != null ? t.abundance.toFixed(2)+'%' : '—'}</td>
      <td class="bar-cell">
        <div class="bar-wrap">
          <div class="bar-fill" style="width:${barPct}%"></div>
        </div>
      </td>
      ${(CONFIG.centrifuge || CONFIG.kaiju) && CONFIG.assembly
        ? `<td style="font-size:10px;font-family:monospace;color:#6b7280">${contigDisplay}</td>`
        : ''}
    </tr>`;
  });

  html += '</tbody></table>';

  // Cross-reference note
  if (CONFIG.centrifuge && CONFIG.assembly) {
    html += `<div style="margin-top:8px;font-size:11px;color:#9ca3af;
      border-top:1px solid #f3f4f6;padding-top:8px">
      ⓘ Contigs shown are Centrifuge-assigned — indicative only.
      Use contig IDs to manually cross-reference with AMR and Plasmids tabs.
    </div>`;
  }

  // Contig search box
  let searchHtml = toggleHtml;
  if ((CONFIG.centrifuge || CONFIG.kaiju) && CONFIG.assembly) {
    searchHtml += `<div style="margin-bottom:12px;display:flex;gap:8px;align-items:center">
      <span style="font-size:11px;font-weight:600;color:#6b7280">Search contig:</span>
      <input type="text" id="contig-search-input"
        placeholder="e.g. contig_572"
        style="padding:5px 9px;border:1px solid #d1d5db;border-radius:5px;
               font-size:12px;font-family:monospace;min-width:180px"
        oninput="searchContig(this.value)">
      <span id="contig-search-result" style="font-size:11px;color:#6b7280"></span>
    </div>`;
  }

  el.innerHTML = targetHtml + searchHtml + html;
  updateToggleButtons();
}

function searchContig(query) {
  const q = (query || '').trim().toLowerCase();
  const resultEl = document.getElementById('contig-search-result');
  const rows = document.querySelectorAll('#tab-taxonomy .tax-row');

  // Reset all highlights
  rows.forEach(r => {
    r.style.background = '';
    r.style.fontWeight = '';
    r.style.boxShadow = '';
  });

  if (!q) { if (resultEl) resultEl.textContent = ''; return; }

  const fullContigSp = taxSource === "kaiju"
    ? (DATA.kaiju_contig_species || {})[currentSample] || {}
    : (DATA.contig_species || {})[currentSample] || {};
  const assignedSp   = fullContigSp[query] || fullContigSp[query.toLowerCase()] || '';

  // Check if contig is in any top-25 species
  const taxa = taxSource === "kaiju"
    ? (DATA.kaiju[currentSample] || [])
    : (DATA.taxonomy[currentSample] || []);
  let foundInTop = false;

  rows.forEach((row, i) => {
    const taxon  = taxa[i];
    if (!taxon) return;
    const contigs = (taxon.contigs || []).map(c => c.toLowerCase());
    if (contigs.includes(q)) {
      row.style.background = '#bfdbfe';
      row.style.fontWeight = '700';
      row.style.boxShadow = 'inset 3px 0 0 #2563eb';
      row.scrollIntoView({ block: 'nearest' });
      foundInTop = true;
    }
  });

  if (foundInTop) {
    if (resultEl) resultEl.innerHTML =
      `<span style="color:#065f46">✓ Found in table below</span>`;
  } else if (assignedSp) {
    if (resultEl) resultEl.innerHTML =
      `<span style="color:#92400e">⚠ Assigned to <em>${esc(assignedSp)}</em>
       — not in top 25 species</span>`;
  } else if (query) {
    if (resultEl) resultEl.innerHTML =
      `<span style="color:#9ca3af">— Not found in Centrifuge classification</span>`;
  }
}

// ── Assembly tab ─────────────────────────────────────────────────────
function renderAssembly() {
  const el  = document.getElementById('tab-assembly');
  const asm = DATA.assembly[currentSample] || {};

  if (Object.keys(asm).length === 0 || Object.values(asm).every(v => v == null)) {
    if (CONFIG.assembly === false) {
      el.innerHTML = emptyState(
        'Assembly was skipped.<br>Enable <code>--perform_assembly</code> to get assembly statistics, ' +
        'plasmid detection and contig-level analysis.'
      );
    } else {
      el.innerHTML = emptyState('No assembly data.<br>Assembly may have failed for this sample.');
    }
    return;
  }

  function sc(label, val, unit) {
    const v = val != null ? fmtN(val) : '—';
    return `<div class="stat-card">
      <div class="stat-lbl">${label}</div>
      <div class="stat-val">${v}<span class="stat-unit">${unit ? ' '+unit : ''}</span></div>
    </div>`;
  }

  let html = '<div class="stat-grid">';
  html += sc('Contigs',        asm.contigs,      '');
  html += sc('Total Length',   asm.total_length, 'bp');
  html += sc('N50',            asm.n50,          'bp');
  html += sc('N90',            asm.n90,          'bp');
  html += sc('Largest Contig', asm.largest,      'bp');
  if (asm.gc != null) html += sc('GC Content', asm.gc, '%');
  html += '</div>';

  el.innerHTML = html;
}

// ── Virulence tab ────────────────────────────────────────────────────
function renderVF() {
  const el    = document.getElementById('tab-vf');
  const allVF = DATA.vf[currentSample] || [];

  if (!CONFIG.abricate) {
    el.innerHTML = emptyState(
      'Virulence detection was not run.<br>' +
      'Enable <code>--run_abricate</code> to detect virulence factors.'
    );
    return;
  }

  if (allVF.length === 0) {
    el.innerHTML = emptyState('No virulence factors detected for this sample.');
    return;
  }

  // Group by category, sort by category name
  const groups = {};
  allVF.forEach(e => {
    const cat = e.category || 'Other Virulence';
    (groups[cat] = groups[cat] || []).push(e);
  });

  const sortedCats = Object.keys(groups).sort();

  let html = '';

  sortedCats.forEach(cat => {
    const entries = groups[cat].sort((a,b) => b.identity - a.identity);
    html += `<div class="dc-group">
      <div class="dc-hdr" onclick="this.nextElementSibling.style.display=
        this.nextElementSibling.style.display==='none'?'':'none'">
        <span>▸</span>
        <span>${esc(cat)}</span>
        <span class="dc-cnt">${entries.length}</span>
      </div>
      <div class="dc-body">
        <table class="amr-tbl">
          <thead><tr>
            <th>Gene</th>
            <th>VF Name</th>
            <th style="width:80px">Identity</th>
            <th style="width:80px">Coverage</th>
            ${CONFIG.centrifuge ? '<th style="width:110px">Contig</th>' : ''}
          </tr></thead>
          <tbody>`;

    entries.forEach(e => {
      html += `<tr>
        <td>
          <div class="gene-name">${esc(e.gene || '—')}</div>
          ${e.product ? `<div class="gene-full" title="${esc(e.product)}">${
            esc(e.product.length > 60 ? e.product.slice(0,58)+'…' : e.product)
          }</div>` : ''}
        </td>
        <td style="font-weight:600;font-size:12px">${esc(e.vf_name || '—')}</td>
        <td class="${idClass(e.identity)}">${e.identity ? e.identity.toFixed(1)+'%' : '—'}</td>
        <td style="color:#6b7280">${e.coverage ? e.coverage.toFixed(1)+'%' : '—'}</td>
        ${CONFIG.centrifuge ? `<td style="font-size:10.5px;font-family:monospace;
          color:#6b7280">${esc(e.contig || '—')}</td>` : ''}
      </tr>`;
    });

    html += '</tbody></table></div></div>';
  });

  // Cross-reference note
  if (CONFIG.centrifuge) {
    html += `<div style="margin-top:10px;font-size:11px;color:#9ca3af;
      border-top:1px solid #f3f4f6;padding-top:8px">
      ⓘ Use the <strong>Contig</strong> column to cross-reference with
      AMR and Plasmids tabs. A virulence gene on a plasmid contig may
      indicate mobile virulence.
    </div>`;
  }

  el.innerHTML = html;
}

// ── Plasmids tab ─────────────────────────────────────────────────────
function renderPlasmids() {
  const el        = document.getElementById('tab-plasmids');
  const pData     = DATA.plasmids[currentSample] || {};
  const pf        = pData.plasmidfinder || [];
  const pc        = pData.plasclass     || {};
  const pcContigs = pc.contigs          || {};

  let html = '';

  // Assembly required check
  if (CONFIG.assembly === false) {
    el.innerHTML = emptyState(
      'Plasmid detection requires assembly.<br>' +
      'Enable <code>--perform_assembly</code> to run PlasClass and PlasmidFinder.'
    );
    return;
  }

  // Tool explanation banner
  html += `<div style="background:#f0f4fa;border-left:3px solid #0f2744;padding:10px 14px;
    margin-bottom:14px;border-radius:3px;font-size:12px;line-height:1.7">
    <strong>PlasClass</strong> — classifies every contig as plasmid-like or chromosome-like
    based on sequence composition (k-mer signatures).<br>
    <strong>PlasmidFinder</strong> — searches for known plasmid replicons from a curated database.<br>
    <span style="color:#6b7280;font-size:11px">These tools are complementary.
    Discordance may indicate novel/unknown plasmids or misclassification — see status column below.</span>
  </div>`;

  // PlasClass summary cards
  html += `<div class="pc-grid">
    <div class="pc-card">
      <div class="stat-lbl">Plasmid-like contigs</div>
      <div class="stat-val">${fmtN(pc.plasmid || 0)}</div>
      <div style="font-size:10px;color:#9ca3af;margin-top:2px">PlasClass</div>
    </div>
    <div class="pc-card">
      <div class="stat-lbl">Chromosome-like contigs</div>
      <div class="stat-val">${fmtN(pc.chromosome || 0)}</div>
      <div style="font-size:10px;color:#9ca3af;margin-top:2px">PlasClass</div>
    </div>
    ${pc.unknown ? `<div class="pc-card">
      <div class="stat-lbl">Unclassified</div>
      <div class="stat-val">${fmtN(pc.unknown)}</div>
      <div style="font-size:10px;color:#9ca3af;margin-top:2px">PlasClass</div>
    </div>` : ''}
  </div>`;

  // Cross-reference PlasmidFinder hits with PlasClass
  const pfContigIds = new Set(pf.map(p => (p.contig || '').split(' ')[0]));
  const putative    = Object.entries(pcContigs)
    .filter(([cid, cls]) => cls === 'plasmid' && !pfContigIds.has(cid))
    .map(([cid]) => cid);

  const nConfirmed = pf.filter(p => pcContigs[(p.contig||'').split(' ')[0]] === 'plasmid').length;
  const nConflict  = pf.filter(p => pcContigs[(p.contig||'').split(' ')[0]] === 'chromosome').length;

  // Interpretation line
  const interp = [];
  if (nConfirmed > 0)   interp.push(`<span style="color:#065f46">✓ ${nConfirmed} replicon(s) confirmed by both tools</span>`);
  if (nConflict > 0)    interp.push(`<span style="color:#92400e">⚠ ${nConflict} replicon(s) on chromosome-classified contigs</span>`);
  if (putative.length > 0) interp.push(`<span style="color:#1e3a8a">◉ ${putative.length} putative plasmid contig(s) — no known replicon</span>`);
  if (interp.length === 0 && (pc.plasmid||0) === 0)
    interp.push('<span style="color:#9ca3af">No plasmid-like contigs detected</span>');
  if (interp.length === 0 && pf.length === 0 && (pc.plasmid||0) > 0)
    interp.push(`<span style="color:#1e3a8a">◉ ${pc.plasmid} plasmid-like contigs — no known replicons identified</span>`);

  if (interp.length > 0) {
    html += `<div style="margin-bottom:14px;padding:8px 12px;background:#f9fafb;
      border:1px solid #e5e7eb;border-radius:4px;font-size:12px;line-height:1.8">
      <strong>Summary: </strong>${interp.join(' &nbsp;·&nbsp; ')}
    </div>`;
  }

  // PlasmidFinder table with PlasClass cross-reference
  html += `<div class="section-title">PlasmidFinder — ${pf.length} replicon hit${pf.length !== 1 ? 's' : ''}</div>`;
  if (pf.length === 0) {
    html += '<p style="color:#9ca3af;font-size:12px;margin-bottom:14px">No plasmid replicons detected by PlasmidFinder.</p>';
  } else {
    html += `<table class="pf-tbl" style="margin-bottom:16px">
      <thead><tr>
        <th>Plasmid replicon</th><th>Identity</th><th>Database</th>
        <th>Contig</th><th>PlasClass</th><th>Status</th>
      </tr></thead><tbody>`;

    pf.forEach(p => {
      const cid     = (p.contig || '').split(' ')[0];
      const plasCls = pcContigs[cid];
      let statusBadge, classBadge;

      if (plasCls === 'plasmid') {
        statusBadge = '<span style="background:#d1fae5;color:#065f46;padding:2px 8px;border-radius:4px;font-size:10px;font-weight:700">✓ Confirmed</span>';
        classBadge  = '<span style="background:#d1fae5;color:#065f46;padding:1px 6px;border-radius:3px;font-size:10px">plasmid</span>';
      } else if (plasCls === 'chromosome') {
        statusBadge = '<span style="background:#fef3c7;color:#92400e;padding:2px 8px;border-radius:4px;font-size:10px;font-weight:700">⚠ Conflicting</span>';
        classBadge  = '<span style="background:#fef3c7;color:#92400e;padding:1px 6px;border-radius:3px;font-size:10px">chromosome</span>';
      } else {
        statusBadge = '<span style="background:#e5e7eb;color:#6b7280;padding:2px 8px;border-radius:4px;font-size:10px">Not in PlasClass</span>';
        classBadge  = '<span style="color:#9ca3af;font-size:10px">—</span>';
      }

      html += `<tr>
        <td><strong>${esc(p.plasmid)}</strong></td>
        <td>${esc(p.identity)}</td>
        <td>${esc(p.database)}</td>
        <td style="font-size:10.5px;color:#6b7280">${esc(cid)}</td>
        <td>${classBadge}</td>
        <td>${statusBadge}</td>
      </tr>`;
    });
    html += '</tbody></table>';
  }

  // Putative plasmids — PlasClass=plasmid but no PlasmidFinder hit
  if (putative.length > 0) {
    html += `<div class="section-title">Putative plasmids — ${putative.length} contig${putative.length !== 1 ? 's' : ''} (PlasClass only)</div>`;
    html += `<p style="font-size:11.5px;color:#6b7280;margin-bottom:8px;line-height:1.6">
      These contigs were classified as plasmid-like by PlasClass but have no matching
      replicon in PlasmidFinder. They may carry novel or uncharacterised plasmids,
      or represent PlasClass misclassification.
    </p>`;
    html += `<div style="font-size:11px;font-family:monospace;background:#f9fafb;
      border:1px solid #e5e7eb;border-radius:4px;padding:10px;
      max-height:140px;overflow-y:auto;color:#374151">`;
    putative.slice(0, 100).forEach(cid => { html += esc(cid) + '<br>'; });
    if (putative.length > 100)
      html += `<span style="color:#9ca3af">… and ${putative.length - 100} more</span>`;
    html += '</div>';
  }

  el.innerHTML = html;
}

// ── Summary section toggle ───────────────────────────────────────────
let summaryOpen = true;
function toggleSummary() {
  summaryOpen = !summaryOpen;
  document.getElementById('summary-body').style.display = summaryOpen ? '' : 'none';
  const hdr = document.getElementById('summary-hdr');
  hdr.classList.toggle('collapsed', !summaryOpen);
  document.getElementById('summary-caret').textContent = summaryOpen ? '▾' : '▸';
}

// ── Tool status badges in header ────────────────────────────────────
function buildToolBadges() {
  const tools = [
    {key:'hamronization', label:'hAMRonize'},
    {key:'rgi',           label:'RGI'},
    {key:'amrfinderplus', label:'AMRFinder'},
    {key:'resfinder',     label:'ResFinder'},
    {key:'abricate',      label:'Abricate'},
    {key:'plasmidfinder', label:'PlasmidFinder'},
    {key:'plasclass',     label:'PlasClass'},
    {key:'centrifuge',    label:'Centrifuge'},
    {key:'kaiju',         label:'Kaiju'},
    {key:'assembly',      label:'Assembly'},
  ];
  const el = document.getElementById('hdr-tools');
  if (!el) return;
  el.innerHTML = tools
    .filter(t => CONFIG[t.key])
    .map(t => `<span style="background:rgba(255,255,255,0.15);padding:2px 7px;
      border-radius:3px;white-space:nowrap">${t.label}</span>`)
    .join('');
}

// ── Init ─────────────────────────────────────────────────────────────
buildSummaryTable();
buildSampleList();
buildToolBadges();

// Auto-select first sample
if (DATA.samples.length > 0) {
  selectSample(DATA.samples[0]);
}


</script>
</body>
</html>"""


def generate_html(data):
    """Fill the HTML template with pipeline data and return the complete HTML string."""
    n_samples  = len(data["samples"])
    n_amr      = sum(len(v) for v in data["amr"].values())
    n_critical = sum(1 for s in data["samples"] if data["status"].get(s) == STATUS_CRITICAL)

    html = HTML_TEMPLATE
    html = html.replace("__RUN_DATE__",          data.get("run_date", ""))
    html = html.replace("__N_SAMPLES__",         str(n_samples))
    html = html.replace("__N_AMR__",             str(n_amr))
    html = html.replace("__N_CRITICAL__",        str(n_critical))
    html = html.replace("__CRITICAL_CLASS__",    "alarm" if n_critical > 0 else "")
    html = html.replace("__DATA_JSON__",         json.dumps(data, ensure_ascii=False))
    html = html.replace("__CRITICAL_CLASSES_JSON__",
                        json.dumps(sorted(CRITICAL_DRUG_CLASSES)))
    return html


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────




def parse_target_classification(results_dir, sample):
    """
    Parse target species classification summary.
    Returns list of detected species with confidence.
    """
    filepath = os.path.join(
        results_dir, "target_species", "classification",
        f"{sample}.species_summary.txt"
    )
    results = []

    if not os.path.exists(filepath):
        return results

    try:
        with open(filepath, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5 or parts[0] == "TaxID":
                    continue
                try:
                    count = int(parts[2])
                except ValueError:
                    count = 0
                results.append({
                    "taxid":      parts[0].strip(),
                    "species":    parts[1].strip(),
                    "count":      count,
                    "status":     parts[3].strip(),
                    "confidence": parts[4].strip(),
                })
        results.sort(key=lambda x: x["count"], reverse=True)
    except Exception as e:
        print(f"[generate_report] WARNING: parse_target_classification failed for {sample}: {e}")

    return results


def parse_contig_species(results_dir, sample):
    """
    Parse Centrifuge contig-level species assignment.
    Returns dict: {contig_id -> species_name}
    """
    filepath = os.path.join(
        results_dir, "centrifuge", sample,
        f"{sample}_contigs_species.tsv"
    )
    mapping = {}

    if not os.path.exists(filepath):
        return mapping

    try:
        with open(filepath, encoding="utf-8") as fh:
            lines = [l for l in fh if not l.startswith("#")]
        import io
        with io.StringIO("".join(lines)) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig  = (row.get("Contig_ID", "") or "").strip()
                species = (row.get("Species",   "") or "").strip()
                if contig and species and species.lower() not in ("unknown", ""):
                    mapping[contig] = species
    except Exception as e:
        print(f"[generate_report] WARNING: parse_contig_species failed for {sample}: {e}")

    return mapping



# VF category keyword mapping
VF_CATEGORIES = [
    ("Toxins & Enzymes",    ["toxin", "hemolysin", "leukotoxin", "enterotoxin",
                              "cytolysin", "hyaluronidase", "protease", "lipase",
                              "phospholipase", "collagenase", "nuclease", "rtx"]),
    ("Adherence & Biofilm", ["pili", "pilus", "fimbri", "adhesin", "efaa", "ebp",
                              "biofilm", "attachment", "agglutinin", "hemagglutinin"]),
    ("Immune Evasion",      ["capsule", "serum", "complement", "immune", "evasion",
                              "opsonin", "antiphagocytic", "m protein"]),
    ("Iron Acquisition",    ["iron", "siderophore", "ferric", "heme", "hemoglobin",
                              "transferrin", "isd", "lactoferrin", "aerobactin"]),
    ("Invasion",            ["invasin", "invasion", "internalin", "penetrat"]),
    ("Secretion System",    ["secretion", "effector", "type iii", "type iv",
                              "type vi", "t3ss", "t4ss", "t6ss"]),
    ("Regulation",          ["regulator", "transcriptional", "two-component",
                              "sensor kinase", "response regulator"]),
]


def _extract_vf_name(product):
    """Extract VF group name from VFDB product string.
    e.g. '(efaA) desc [EfaA (VF0354)] [Organism]' → 'EfaA'
    """
    m = re.search(r"\[([^\[]+?)\s*\(VF\d+\)\]", product)
    return m.group(1).strip() if m else ""


def _categorize_vf(vf_name, product):
    """Map VF name/product to clinical functional category."""
    text = (vf_name + " " + product).lower()
    for category, keywords in VF_CATEGORIES:
        if any(k in text for k in keywords):
            return category
    return "Other Virulence"


def parse_abricate_vf(results_dir, sample):
    """
    Parse Abricate VFDB output for virulence factor detection.
    Returns list of VF entry dicts.
    """
    filepath = os.path.join(results_dir, "abricate", sample, f"{sample}.txt")
    results  = []

    if not os.path.exists(filepath):
        return results

    try:
        with open(filepath, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 14:
                    continue

                gene      = parts[5].strip()
                contig    = parts[1].strip()
                product   = parts[13].strip()
                database  = parts[11].strip().lower()
                accession = parts[12].strip()

                if database != "vfdb":
                    continue

                try:
                    identity = round(float(parts[10].strip()), 1)
                except (ValueError, IndexError):
                    identity = 0.0
                try:
                    coverage = round(float(parts[9].strip()), 1)
                except (ValueError, IndexError):
                    coverage = 0.0

                vf_name  = _extract_vf_name(product)
                category = _categorize_vf(vf_name, product)

                results.append({
                    "gene":      gene,
                    "vf_name":   vf_name or gene,
                    "category":  category,
                    "product":   product,
                    "identity":  identity,
                    "coverage":  coverage,
                    "contig":    contig,
                    "accession": accession,
                })
    except Exception as e:
        print(f"[generate_report] WARNING: parse_abricate_vf failed for {sample}: {e}")

    return results


def _build_identity_index(results_dir, sample):
    """
    Read raw tool outputs to get identity % missing from hAMRonization.
    Returns dict: {gene_name_lower: identity_float}
    """
    index = {}

    # RGI — Best_Identities column
    for fpath in glob.glob(os.path.join(results_dir, "rgi", sample, "*.txt")):
        try:
            with open(fpath, encoding="utf-8") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    gene = (row.get("Best_Hit_ARO", "") or "").strip().lower()
                    val  = row.get("Best_Identities", "") or ""
                    if gene and val:
                        index[gene] = round(float(val), 1)
        except Exception:
            pass

    # AMRFinderPlus — % Identity to reference column
    for fpath in glob.glob(os.path.join(results_dir, "amrfinderplus", sample, "*.tsv")):
        try:
            with open(fpath, encoding="utf-8") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    gene = (row.get("Element symbol", "") or "").strip().lower()
                    val  = row.get("% Identity to reference", "") or ""
                    if gene and val:
                        index[gene] = round(float(val), 1)
        except Exception:
            pass

    # Abricate — %IDENTITY column
    for fpath in glob.glob(os.path.join(results_dir, "abricate", sample, "*.txt")):
        try:
            with open(fpath, encoding="utf-8") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    gene = (row.get("GENE", "") or "").strip().lower()
                    val  = row.get("%IDENTITY", "") or ""
                    if gene and val:
                        index[gene] = round(float(val), 1)
        except Exception:
            pass

    return index


def main():
    parser = argparse.ArgumentParser(
        description="Generate metaAMR-Plus clinical HTML report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results_dir", required=True,
        help="Path to pipeline results directory (e.g. results/)"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for the HTML report"
    )
    parser.add_argument(
        "--run_name", default="",
        help="Optional run label shown in the report"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.results_dir):
        sys.exit(f"ERROR: results_dir not found: {args.results_dir}")

    print(f"[generate_report] Results dir : {args.results_dir}")
    print(f"[generate_report] Output dir  : {args.outdir}")

    # ── Parse data ──────────────────────────────────────────────────────────
    print("[generate_report] Parsing hAMRonization…")
    all_amr = parse_hamronization(args.results_dir)

    print("[generate_report] Discovering samples…")
    samples = discover_samples(args.results_dir, all_amr)
    print(f"[generate_report]   Found {len(samples)} sample(s): {', '.join(samples)}")

    run_config = detect_run_config(args.results_dir)
    print("[generate_report] Parsing per-sample data…")
    taxonomy                = {}
    kaiju_taxonomy          = {}
    kaiju_contig_species_map = {}
    assembly                = {}
    plasmids           = {}
    vf                 = {}
    status             = {}
    contig_species_map      = {}
    target_classification  = {}

    for s in samples:
        taxonomy[s]              = parse_centrifuge(args.results_dir, s)
        target_classification[s] = parse_target_classification(args.results_dir, s)
        kaiju_taxonomy[s]           = parse_kaiju(args.results_dir, s)
        kaiju_contig_species_map[s] = parse_kaiju_contigs(args.results_dir, s)
        # Add contig lists to Kaiju taxonomy entries
        kj_sp_contigs = {}
        for cid, sp in kaiju_contig_species_map[s].items():
            kj_sp_contigs.setdefault(sp, []).append(cid)
            parts = sp.split()
            if len(parts) > 2:
                kj_sp_contigs.setdefault(' '.join(parts[:2]), []).append(cid)
        for taxon in kaiju_taxonomy[s]:
            taxon['contigs'] = kj_sp_contigs.get(taxon.get('name', ''), [])

        # Add contig lists to taxonomy entries (contig → species inverted)
        contig_sp = parse_contig_species(args.results_dir, s)
        contig_species_map[s] = contig_sp
        sp_contigs = {}
        for cid, sp in contig_sp.items():
            sp_contigs.setdefault(sp, []).append(cid)
            # Also map strain-level names to parent species
            parts = sp.split()
            if len(parts) > 2:
                sp_contigs.setdefault(" ".join(parts[:2]), []).append(cid)
        for taxon in taxonomy[s]:
            taxon["contigs"] = sp_contigs.get(taxon.get("name", ""), [])

        assembly[s] = parse_quast(args.results_dir, s)
        vf[s]       = parse_abricate_vf(args.results_dir, s)
        plasmids[s] = {
            "plasmidfinder": parse_plasmidfinder(args.results_dir, s),
            "plasclass":     parse_plasclass(args.results_dir, s),
        }
        rf = parse_resfinder(args.results_dir, s)
        if rf:
            all_amr[s] = all_amr.get(s, []) + rf
        status[s]   = get_sample_status(
            all_amr.get(s, []), assembly[s],
            assembly_skipped=not run_config["assembly"]
        )

    # Fill in missing identity values from raw tool outputs
    print("[generate_report] Filling missing identity values...")
    for s in samples:
        idx = _build_identity_index(args.results_dir, s)
        for entry in all_amr.get(s, []):
            if entry["identity"] == 0.0 and entry["gene"]:
                key = entry["gene"].lower()
                if key in idx:
                    entry["identity"] = idx[key]

    drug_classes, summary_matrix = build_summary_matrix(samples, all_amr)

    # ── Assemble data object ─────────────────────────────────────────────────
    report_data = {
        "run_date":       datetime.now().strftime("%Y-%m-%d"),
        "run_name":       args.run_name,
        "pipeline_ver":   VERSION,
        "config":         run_config,
        "samples":        samples,
        "drug_classes":   drug_classes,
        "summary_matrix": summary_matrix,
        "amr":            {s: all_amr.get(s, []) for s in samples},
        "vf":             vf,
        "kaiju":          kaiju_taxonomy,
        "target_classification": target_classification,
        "kaiju_contig_species": kaiju_contig_species_map,
        "contig_species": contig_species_map,
        "taxonomy":       taxonomy,
        "assembly":       assembly,
        "plasmids":       plasmids,
        "status":         status,
    }

    # ── Generate HTML ────────────────────────────────────────────────────────
    print("[generate_report] Generating HTML…")
    html = generate_html(report_data)

    os.makedirs(args.outdir, exist_ok=True)
    outfile = os.path.join(args.outdir, "metaamr_report.html")
    with open(outfile, "w", encoding="utf-8") as fh:
        fh.write(html)

    size_kb = os.path.getsize(outfile) / 1024
    print(f"[generate_report] Done! Report written to: {outfile}  ({size_kb:.1f} KB)")
    print(f"[generate_report]   Samples : {len(samples)}")
    print(f"[generate_report]   AMR genes: {sum(len(v) for v in all_amr.values())}")
    print(f"[generate_report]   Critical : {sum(1 for s in samples if status[s] == STATUS_CRITICAL)}")


if __name__ == "__main__":
    main()
