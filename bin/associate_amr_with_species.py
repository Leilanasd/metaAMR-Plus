#!/usr/bin/env python3

import sys


def load_species_info(summary_file):
    species_counts = {}
    total = 0
    with open(summary_file) as f:
        for line in f:
            if line.lower().startswith("taxid"):  # Skip header
                continue
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            taxid, name, count, status = parts
            count = int(count)
            if status == "Present":
                species_counts[name] = count
                total += count
    return species_counts, total


def associate_amr(amr_file, summary_file):
    species_counts, total_reads = load_species_info(summary_file)

    with open(amr_file) as f:
        for line in f:
            line = line.rstrip('\n')

            if not line.strip():
                continue

            # Header line
            if line.startswith('#') or line.lower().startswith('resistance gene'):
                print(line + '\tPredicted Species')
            else:
                if total_reads > 0 and species_counts:
                    entries = list(species_counts.keys())
                    print(line + '\t' + ', '.join(entries))
                else:
                    print(line + '\tNA')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: associate_amr_with_species.py <amr_file> <species_summary_file>")
        sys.exit(1)

    associate_amr(sys.argv[1], sys.argv[2])