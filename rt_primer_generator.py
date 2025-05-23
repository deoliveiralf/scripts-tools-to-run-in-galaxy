#!/usr/bin/env python

import argparse
import sys

def generate_stem_loop_primer(mirna_sequence):
    """
    Generate a Universal stem-loop RT primer based on Varkonyi-Gasic et al. 2007.
    
    Args:
        mirna_sequence (str): The mature miRNA sequence
    
    Returns:
        str: The complete Universal stem-loop RT primer
    """
    # The backbone sequence
    backbone = "5'-gTCgTATCCAgTgCAgTTAgCCAggTCgACgTgATCCAgTACgAC-3'"
    
    # Clean the backbone sequence (remove 5'- and -3' notation)
    backbone_clean = backbone.replace("5'-", "").replace("-3'", "")
    
    # Extract the last 6 nucleotides from the miRNA sequence
    if len(mirna_sequence) < 6:
        raise ValueError("miRNA sequence must be at least 6 nucleotides long")
    
    last_six_nt = mirna_sequence[-6:]
    
    # Generate the complementary sequence to the last 6 nucleotides
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    
    # Generate the reverse complementary sequence
    complementary_seq = ''.join(complement_map.get(nt, nt) for nt in last_six_nt[::-1])
    
    # Add the complementary sequence to the 3' end of the backbone
    full_primer = f"5'-{backbone_clean}{complementary_seq}-3'"
    
    return full_primer

def main():
    # Set up argument parser for Galaxy
    parser = argparse.ArgumentParser(description="Generate Universal stem-loop RT primer for miRNA")
    parser.add_argument('--input_source', required=True, choices=['file', 'direct'], 
                        help='Source of miRNA sequence (file or direct input)')
    parser.add_argument('--input_file', help='Input file containing miRNA sequence')
    parser.add_argument('--input_text', help='Direct input of miRNA sequence')
    parser.add_argument('--output', required=True, help='Output file for RT primer')
    parser.add_argument('--details', required=True, help='Output file for detailed information')
    
    args = parser.parse_args()
    
    # Get miRNA sequence from either file or direct input
    if args.input_source == 'file':
        if not args.input_file:
            sys.stderr.write("Error: Input file is required when input source is 'file'\n")
            sys.exit(1)
        with open(args.input_file, 'r') as input_file:
            mirna_seq = input_file.read().strip()
    else:  # direct input
        if not args.input_text:
            sys.stderr.write("Error: Input text is required when input source is 'direct'\n")
            sys.exit(1)
        mirna_seq = args.input_text
    
    # Clean input (remove whitespace and any notations)
    mirna_seq = mirna_seq.strip().replace("5'-", "").replace("-3'", "")
    
    try:
        primer = generate_stem_loop_primer(mirna_seq)
        
        # Write primer to output file
        with open(args.output, 'w') as output_file:
            output_file.write(f"{primer}\n")
        
        # Get the last 6 nt and complementary sequence for detailed output
        last_six = mirna_seq[-6:]
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                          'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        complementary_seq = ''.join(complement_map.get(nt, nt) for nt in last_six[::-1])
        
        # Write detailed information to details file
        with open(args.details, 'w') as details_file:
            details_file.write("Universal stem-loop RT primer generator (Varkonyi-Gasic et al. 2007)\n")
            details_file.write("--------------------------------------------------------------\n")
            details_file.write(f"Input miRNA sequence: {mirna_seq}\n\n")
            details_file.write("Generated Universal stem-loop RT primer:\n")
            details_file.write(f"{primer}\n\n")
            details_file.write(f"Last 6 nucleotides of miRNA: {last_six}\n")
            details_file.write(f"Complementary sequence added: {complementary_seq}\n")
            details_file.write(f"Position: Added to 3' end of the backbone\n")
        
    except ValueError as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
