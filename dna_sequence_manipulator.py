#!/usr/bin/env python

import argparse
import sys

def reverse_sequence(sequence):
    """
    Takes a DNA sequence and returns it reversed.
    """
    return sequence[::-1]

def complement_sequence(sequence):
    """
    Takes a DNA sequence and returns its complement.
    """
    # Define the complement mapping
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n',
                  'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
                  'r': 'y', 'y': 'r', 'm': 'k', 'k': 'm',
                  'S': 'S', 'W': 'W', 'H': 'D', 'D': 'H',
                  's': 's', 'w': 'w', 'h': 'd', 'd': 'h',
                  'B': 'V', 'V': 'B', 's': 's', 'w': 'w',
                  'b': 'v', 'v': 'b'}
    
    # Generate the complement
    comp = ''.join(complement.get(base, base) for base in sequence)
    
    return comp

def reverse_complement(sequence):
    """
    Takes a DNA sequence and returns its reverse complement.
    """
    # First complement, then reverse
    return reverse_sequence(complement_sequence(sequence))

def parse_fasta(input_file):
    """
    Parse a FASTA file and return sequences with their identifiers.
    """
    sequences = []
    current_id = None
    current_seq = []
    
    for line in input_file:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save the previous sequence if it exists
            if current_id is not None:
                sequences.append((current_id, ''.join(current_seq)))
            
            # Start a new sequence
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    
    # Don't forget to add the last sequence
    if current_id is not None:
        sequences.append((current_id, ''.join(current_seq)))
    
    return sequences

def main():
    parser = argparse.ArgumentParser(description='Manipulate DNA sequences: reverse, complement, or reverse-complement.')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='Input file in FASTA format (default: stdin)')
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output file (default: stdout)')
    parser.add_argument('--operation', choices=['reverse', 'complement', 'reverse-complement'], 
                        default='reverse-complement', help='Operation to perform (default: reverse-complement)')
    
    args = parser.parse_args()
    
    sequences = parse_fasta(args.input)
    
    for seq_id, sequence in sequences:
        if args.operation == 'reverse':
            result = reverse_sequence(sequence)
            operation_name = "Reverse"
        elif args.operation == 'complement':
            result = complement_sequence(sequence)
            operation_name = "Complement"
        else:  # reverse-complement
            result = reverse_complement(sequence)
            operation_name = "Reverse Complement"
            
        args.output.write(f'>{seq_id} [{operation_name}]\n')
        
        # Write sequence in chunks of 60 bases per line
        for i in range(0, len(result), 60):
            args.output.write(f'{result[i:i+60]}\n')
    
    args.input.close()
    args.output.close()

if __name__ == '__main__':
    main()
