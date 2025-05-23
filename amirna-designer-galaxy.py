#!/usr/bin/env python3
"""
Artificial miRNA Designer for Galaxy
------------------------------------
A Galaxy tool to design artificial microRNAs (amiRNAs) for targeted gene silencing,
similar to the functionality of http://wmd3.weigelworld.org/cgi-bin/webapp.cgi

This script:
1. Takes a target sequence as input
2. Identifies potential amiRNA sites based on empirical design rules
3. Scores candidate sequences based on hybridization properties
4. Returns optimized amiRNA sequences for cloning

Modified version that doesn't require BioPython
"""

import sys
import re
import random
import os

# Constants for miRNA design - Default values
DEFAULT_MIR319A_BACKBONE = """CTTGCATATGTAGGCGTCTAATTAAGTGCTATGGTGAAGAAGATGAGCATGATTTGATTCATTTCATTGACTATGGAGAGAAAGGTTGTTTTCT
CAATACAGGTGAATCTTTCATTGAAATTGCTTACAGTGTTATGCATAATAAAAGAGAATGAAGTAGTTAATGAAAAGAAAGACATTTTACTTCG
TTTCTATGTGATCTTTAAGAGTGGGACTTTGGCCTTAGTTTTGTTGATTCATTAGATGGTTTTTGCTTTGTATCTGTGGGAAAAGAATTGAAGG
GTGTTATTTTAGTACAAAGTAGAGTGTTTCAGAAACAATTGTTGTGTTGTAAAACTGTTAAGAAATCTTTGATTTCTTTTTTACC"""

DEFAULT_AMIRNA_5P_SEQ = "GTTGTTGTCCATAGTCGACTTAAGGTA"
DEFAULT_AMIRNA_3P_SEQ = "TCTAGACAAATTCGTAGTTTAAGGTA"

# Fallback calculation function when ViennaRNA is not available
def calculate_hybridization_energy_simple(mirna, target):
    """
    A simplified algorithm to estimate hybridization energy when ViennaRNA is not available.
    This is a rough approximation based on base pairing rules.
    """
    mirna = mirna.upper()
    target = target.upper()

    # Score different types of base pairings
    energy_contribution = {
        'GC': -3.0,  # G-C base pair (strongest)
        'CG': -3.0,  # C-G base pair
        'AU': -2.0,  # A-U base pair
        'UA': -2.0,  # U-A base pair
        'GU': -1.0,  # G-U wobble pair
        'UG': -1.0,  # U-G wobble pair
        'MM': 0.5    # Mismatch (unfavorable)
    }

    total_energy = 0

    # Calculate energy based on simple base pairing
    for i in range(min(len(mirna), len(target))):
        pair = f"{mirna[i]}{target[i]}"
        if pair in energy_contribution:
            total_energy += energy_contribution[pair]
        else:
            total_energy += energy_contribution['MM']

    # Add a penalty for length mismatch
    if len(mirna) != len(target):
        total_energy += abs(len(mirna) - len(target)) * 1.0

    return total_energy

# Custom functions to replace BioPython functionality
def parse_fasta(fasta_file):
    """Custom FASTA parser function."""
    records = []
    current_id = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    records.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]  # Get sequence ID
                current_seq = []
            elif line and current_id:  # Sequence line
                current_seq.append(line)
    
    # Add the last record
    if current_id:
        records.append((current_id, ''.join(current_seq)))
        
    return records

def get_complement(sequence):
    """Get the complement of a DNA/RNA sequence."""
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a'
    }
    return ''.join(complement_map.get(base, base) for base in sequence)

# Try to import RNA, but provide a fallback if not available
try:
    import RNA  # ViennaRNA package
    calculate_hybridization_energy = lambda mirna, target: RNA.duplexfold(mirna, target).energy
    has_vienna = True
except ImportError:
    calculate_hybridization_energy = calculate_hybridization_energy_simple
    has_vienna = False
    sys.stderr.write("Warning: ViennaRNA package not found. Using simplified scoring method.\n")

def check_target_sequence(sequence):
    """Validates if the input sequence contains only valid nucleotides."""
    if not re.match(r'^[ATGCNatgcn]+$', sequence):
        return False
    return True

def find_potential_targets(sequence, length=21):
    """Find all potential target sites of specified length in the sequence."""
    sequence = sequence.upper()
    potential_targets = []

    for i in range(len(sequence) - length + 1):
        target = sequence[i:i+length]
        # Skip sequences with Ns
        if 'N' in target:
            continue
        potential_targets.append({
            'target_seq': target,
            'position': i,
            'score': 0
        })

    return potential_targets

def design_amirna(target_seq):
    """Design complementary amiRNA based on target sequence."""
    # Convert target to RNA and take complement
    target_rna = target_seq.replace('T', 'U')

    # Create complementary amiRNA sequence with mismatches at positions 1 and 21
    # (following empirical design rules)
    complementary = get_complement(target_rna)

    # Add mismatches at specific positions based on empirical design rules
    amirna_seq = list(complementary)

    # Mismatch at position 1
    amirna_seq[0] = get_mismatch(target_rna[0])

    # Mismatch at position 21 (last position)
    amirna_seq[-1] = get_mismatch(target_rna[-1])

    # G-U wobble at position 20
    if target_rna[-2] == 'A':
        amirna_seq[-2] = 'G'
    elif target_rna[-2] == 'C':
        amirna_seq[-2] = 'A'

    return ''.join(amirna_seq)

def get_mismatch(nucleotide):
    """Return a nucleotide that mismatches with the input nucleotide."""
    nucleotide = nucleotide.upper()
    if nucleotide == 'A':
        return random.choice(['C', 'G'])
    elif nucleotide == 'U':
        return random.choice(['A', 'G'])
    elif nucleotide == 'G':
        return random.choice(['A', 'C'])
    elif nucleotide == 'C':
        return random.choice(['G', 'U'])
    return 'A'  # Default

def score_targets(targets, simple_scoring=False):
    """Score target sites based on design rules."""
    scored_targets = []

    for target in targets:
        target_seq = target['target_seq']
        target_rna = target_seq.replace('T', 'U')

        # Design amiRNA
        amirna = design_amirna(target_seq)

        # Calculate hybridization energy
        energy = calculate_hybridization_energy(amirna, target_rna)

        # Calculate GC content
        gc_content = (target_seq.count('G') + target_seq.count('C')) / len(target_seq)

        # Apply scoring rules (lower is better)
        score = abs(energy)  # Strong binding is good

        # Penalty for extreme GC content
        if gc_content < 0.3 or gc_content > 0.7:
            score -= 5

        # Penalty for targets with runs of 4+ identical nucleotides
        if re.search(r'([ATGC])\1{3,}', target_seq):
            score -= 3

        # Final target information
        scored_target = {
            'target_seq': target_seq,
            'position': target['position'],
            'amirna_seq': amirna,
            'energy': energy,
            'gc_content': gc_content,
            'score': score
        }

        scored_targets.append(scored_target)

    # Sort by score (higher is better)
    return sorted(scored_targets, key=lambda x: x['score'], reverse=True)

def generate_cloning_sequence(amirna_seq, mir_backbone, amirna_5p, amirna_3p):
    """Generate the complete sequence for cloning in miRNA backbone."""
    # Replace the miR319a sequence with our amiRNA
    mir_star_seq = get_complement(amirna_seq).replace('U', 'T')

    # Create the complete insert sequence
    cloning_seq = f"{amirna_5p}{amirna_seq.replace('U', 'T')}{amirna_3p}{mir_star_seq}"

    return cloning_seq

def main():
    # Galaxy passes arguments directly to the script
    try:
        input_fasta = sys.argv[1]  # Input FASTA file
        output_file = sys.argv[2]  # Output file
        num_candidates = int(sys.argv[3])  # Number of amiRNA candidates to return

        # Check if we're using custom or default backbone
        use_custom_backbone = sys.argv[4].lower() == "true"

        if use_custom_backbone:
            mir_backbone = sys.argv[5]
            amirna_5p = sys.argv[6]
            amirna_3p = sys.argv[7]
        else:
            mir_backbone = DEFAULT_MIR319A_BACKBONE
            amirna_5p = DEFAULT_AMIRNA_5P_SEQ
            amirna_3p = DEFAULT_AMIRNA_3P_SEQ

    except IndexError:
        sys.stderr.write("Error: Not enough arguments provided to the script.\n")
        sys.stderr.write("Usage: python amirna-designer-galaxy.py input.fasta output.tabular num_candidates use_custom_backbone [mir_backbone amirna_5p amirna_3p]\n")
        sys.exit(1)

    target_sequences = []

    # Process input FASTA file using our custom parser
    try:
        target_sequences = parse_fasta(input_fasta)
        # Filter out invalid sequences
        target_sequences = [(seq_id, seq) for seq_id, seq in target_sequences if check_target_sequence(seq)]
    except Exception as e:
        sys.stderr.write(f"Error reading FASTA file: {e}\n")
        sys.exit(1)

    if not target_sequences:
        sys.stderr.write("No valid sequences found in input file\n")
        sys.exit(1)

    # Process each sequence
    results = []

    for seq_id, sequence in target_sequences:
        sys.stderr.write(f"Processing sequence: {seq_id}\n")

        # Find potential targets
        potential_targets = find_potential_targets(sequence)

        if not potential_targets:
            sys.stderr.write(f"No valid target sites found in {seq_id}\n")
            continue

        # Score targets
        scored_targets = score_targets(potential_targets, not has_vienna)

        # Get the top N candidates
        top_candidates = scored_targets[:num_candidates]

        for i, candidate in enumerate(top_candidates):
            # Generate cloning sequence
            cloning_seq = generate_cloning_sequence(
                candidate['amirna_seq'],
                mir_backbone,
                amirna_5p,
                amirna_3p
            )

            results.append({
                'sequence_id': seq_id,
                'candidate_num': i+1,
                'target_position': f"{candidate['position']+1}-{candidate['position']+21}",
                'target_sequence': candidate['target_seq'],
                'amirna_sequence': candidate['amirna_seq'],
                'binding_energy': candidate['energy'],
                'gc_content': candidate['gc_content'],
                'cloning_sequence': cloning_seq
            })

    # Write output file
    try:
        with open(output_file, 'w') as out_file:
            # Write header
            out_file.write("Sequence_ID\tCandidate\tTarget_Position\tTarget_Sequence\tamiRNA_Sequence\tBinding_Energy\tGC_Content\tCloning_Sequence\n")

            # Write results
            for result in results:
                out_file.write(f"{result['sequence_id']}\t{result['candidate_num']}\t{result['target_position']}\t{result['target_sequence']}\t")
                out_file.write(f"{result['amirna_sequence']}\t{result['binding_energy']:.2f}\t{result['gc_content']:.2f}\t{result['cloning_sequence']}\n")
    except Exception as e:
        sys.stderr.write(f"Error writing output file: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
