#!/usr/bin/env python3
import sys
import os

def load_species_info(summary_file):
    species_counts = {}
    total = 0
    with open(summary_file) as f:
        for line in f:
            if line.startswith('#'):        # Skip comment lines
                continue
            if line.lower().startswith("taxid"):  # Skip header
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:             # Accept 4 or 5 columns
                continue
            taxid, name, count, status = parts[0], parts[1], parts[2], parts[3]
            count = int(count)
            if status == "Present":
                species_counts[name] = count
                total += count
    return species_counts, total
    
def associate_amr(amr_file, summary_file, target_species=None, sample_id=None):
    species_counts, total_reads = load_species_info(summary_file)
    detected_species = list(species_counts.keys())

    # Filter detected species by target if provided
    if target_species:
        targets = [t.lower().strip() for t in target_species.split(',')]
        relevant_species = [sp for sp in detected_species
                           if any(t in sp.lower().strip() for t in targets)]
    else:
        relevant_species = detected_species

    # Write AMR results file
    species_str = ', '.join(relevant_species) if relevant_species else 'None detected'
    
    with open(amr_file) as f:
        lines = f.readlines()

    # Print header note
    print(f"# AMR genes detected in reads classified as: {species_str}")
    print(f"# Note: gene-to-species attribution requires assembly-based analysis")
    print(f"# Sample: {sample_id or os.path.basename(amr_file)}")
    print(f"# Total classified reads: {total_reads}")
    print("#")

    for line in lines:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        # Print header and data rows without species column
        if line.startswith('#') or line.lower().startswith('resistance gene'):
            print(line)
        else:
            print(line)

if __name__ == "__main__":
    if len(sys.argv) not in (3, 4, 5):
        print("Usage: associate_amr_with_species.py <amr_file> <species_summary_file> [target_species] [sample_id]")
        sys.exit(1)

    target = sys.argv[3] if len(sys.argv) >= 4 else None
    sample = sys.argv[4] if len(sys.argv) == 5 else None
    associate_amr(sys.argv[1], sys.argv[2], target, sample)