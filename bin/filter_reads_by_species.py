#!/usr/bin/env python3

import sys
import re
from collections import defaultdict


def clean_species_name(name):
    """Remove trailing score like ' (1.00)'."""
    return re.sub(r"\s*\([^)]*\)$", "", name).strip()


def normalize_species(name):
    """Make matching tolerant to spaces, underscores, and hyphens."""
    name = clean_species_name(name)
    name = name.lower().replace("_", " ").replace("-", " ")
    name = re.sub(r"\s+", " ", name).strip()
    return name


def get_species_taxids(report_file, species_list):
    taxid_to_name = {}
    entries = []

    with open(report_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            raw_name = parts[0].strip()
            clean_name = clean_species_name(raw_name)
            normalized_name = normalize_species(raw_name)
            taxid = parts[1].strip()
            rank = parts[2].strip().lower()

            # Skip header if present
            if normalized_name == "name" and taxid.lower() == "taxid":
                continue

            # Keep species-level and more specific labels that may represent strains/subspecies
            if rank in {"species", "strain", "subspecies"}:
                entries.append((normalized_name, taxid, clean_name))

    selected_taxids = set()

    for s in species_list:
        target = normalize_species(s)

        matched = False
        for name_norm, taxid, original_name in entries:
            if name_norm == target or name_norm.startswith(target + " "):
                selected_taxids.add(taxid)
                taxid_to_name[taxid] = original_name
                matched = True

        if not matched:
            print(f"WARNING: Species '{s}' not found in report.", file=sys.stderr)

    return selected_taxids, taxid_to_name

def extract_reads(results_file, taxid_set, taxid_to_name, output_file, summary_file):
    species_read_counts = defaultdict(int)

    if not taxid_set:
        print("No matching taxIDs found — creating empty output.", file=sys.stderr)
        with open(output_file, 'w') as out:
            pass
        with open(summary_file, 'w') as summary:
            summary.write("# Confidence: High = >=10 reads, Low = 1-9 reads\n")
            summary.write("# Low confidence detections may represent misclassification or sequencing noise\n")
            summary.write("#\n")
            summary.write("TaxID\tSpecies\tCount\tStatus\tConfidence\n")
            summary.write("NA\tNA\t0\tAbsent\tNA\n")
        return

    with open(results_file) as infile, open(output_file, 'w') as out:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            if parts[0].lower() in {"readid", "read_id", "name"}:
                continue
            read_id = parts[0]
            taxid = parts[2]
            if taxid in taxid_set:
                out.write(f"{read_id}\t{taxid}\n")
                species_read_counts[taxid] += 1

    with open(summary_file, 'w') as summary:
        summary.write("# Confidence: High = >=10 reads, Low = 1-9 reads\n")
        summary.write("# Low confidence detections may represent misclassification or sequencing noise\n")
        summary.write("#\n")
        summary.write("TaxID\tSpecies\tCount\tStatus\tConfidence\n")
        for taxid in sorted(taxid_set):
            name = taxid_to_name.get(taxid, "Unknown")
            count = species_read_counts.get(taxid, 0)
            status = "Present" if count > 0 else "Absent"
            confidence = "High" if count >= 10 else "Low"
            summary.write(f"{taxid}\t{name}\t{count}\t{status}\t{confidence}\n")
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage:", file=sys.stderr)
        print(
            "  python filter_reads_by_species.py <report.txt> <results.txt> <output_reads.txt> <summary.txt> 'Species1,Species2,...'",
            file=sys.stderr,
        )
        sys.exit(1)

    report_file = sys.argv[1]
    results_file = sys.argv[2]
    output_file = sys.argv[3]
    summary_file = sys.argv[4]
    species_input = sys.argv[5]
    species_list = [s.strip() for s in species_input.split(',') if s.strip()]

    selected_taxids, taxid_to_name = get_species_taxids(report_file, species_list)
    extract_reads(results_file, selected_taxids, taxid_to_name, output_file, summary_file)
