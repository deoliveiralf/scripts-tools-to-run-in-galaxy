#!/usr/bin/env python

import argparse
import sys
import os
from collections import defaultdict

# Check if Biopython is available and import it
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: Biopython not found, using built-in FASTA parser")

# IUPAC degenerate nucleotide codes
IUPAC_DEGENERATE = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'U': ['U', 'T'],  # Treat U (RNA) as T (DNA)
    'R': ['A', 'G'],  # Purine
    'Y': ['C', 'T'],  # Pyrimidine
    'S': ['G', 'C'],  # Strong
    'W': ['A', 'T'],  # Weak
    'K': ['G', 'T'],  # Keto
    'M': ['A', 'C'],  # Amino
    'B': ['C', 'G', 'T'],  # Not A
    'D': ['A', 'G', 'T'],  # Not C
    'H': ['A', 'C', 'T'],  # Not G
    'V': ['A', 'C', 'G'],  # Not T
    'N': ['A', 'C', 'G', 'T'],  # Any base
    '-': ['-']  # Gap
}

# Complementary bases for generating antisense (reverse complement) sequences
COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'U': 'A',
    'R': 'Y',  # Purine (A,G) -> Pyrimidine (C,T)
    'Y': 'R',  # Pyrimidine (C,T) -> Purine (A,G)
    'S': 'S',  # Strong (G,C) -> Strong (G,C)
    'W': 'W',  # Weak (A,T) -> Weak (A,T)
    'K': 'M',  # Keto (G,T) -> Amino (A,C)
    'M': 'K',  # Amino (A,C) -> Keto (G,T)
    'B': 'V',  # Not A (C,G,T) -> Not T (A,C,G)
    'D': 'H',  # Not C (A,G,T) -> Not G (A,C,T)
    'H': 'D',  # Not G (A,C,T) -> Not C (A,G,T)
    'V': 'B',  # Not T (A,C,G) -> Not A (C,G,T)
    'N': 'N',  # Any base -> Any base
}

def parse_args():
    parser = argparse.ArgumentParser(description='Exact Search: Find matches of query sequences in a database, supporting degenerate nucleotides')
    parser.add_argument('--query', required=True, help='FASTA file containing query sequences')
    parser.add_argument('--database', required=True, help='FASTA file containing database sequences')
    parser.add_argument('--mismatch', type=int, default=0, help='Number of mismatches allowed (default: 0)')
    parser.add_argument('--gap', type=int, default=0, help='Number of gaps allowed (default: 0)')
    parser.add_argument('--antisense', action='store_true', help='Include antisense (reverse complement) search')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    return parser.parse_args()

def reverse_complement(sequence):
    """
    Generate the reverse complement of a DNA sequence, supporting IUPAC codes
    """
    # Convert to uppercase and reverse the sequence
    sequence = sequence.upper()[::-1]
    
    # Replace each base with its complement
    complement = ''.join(COMPLEMENT.get(base, base) for base in sequence)
    
    return complement

def parse_fasta(file_path):
    """Parse a FASTA file without requiring Biopython"""
    sequences = {}
    current_id = None
    current_seq = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # Save the previous sequence if there was one
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start a new sequence
                    current_id = line[1:].split()[0]  # Get ID (everything after '>' until first space)
                    current_seq = []
                else:
                    # Add sequence line
                    current_seq.append(line)
            
            # Save the last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
    
    except Exception as e:
        print(f"Error parsing FASTA file {file_path}: {e}")
        return {}
    
    return sequences

def print_fasta_info(file_path, debug=False):
    """Print information about sequences in a FASTA file"""
    if debug:
        print(f"Checking file: {file_path}")
        
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} does not exist")
        return False

    try:
        if BIOPYTHON_AVAILABLE:
            sequences = list(SeqIO.parse(file_path, "fasta"))
            if debug:
                print(f"  Number of sequences: {len(sequences)}")
                for i, seq in enumerate(sequences[:3]):  # Print info for first 3 sequences
                    print(f"  Sequence {i+1}: ID={seq.id}, Length={len(seq.seq)}, Sequence: {seq.seq}")
                if len(sequences) > 3:
                    print(f"  ... and {len(sequences)-3} more sequences")
            return len(sequences) > 0
        else:
            sequences = parse_fasta(file_path)
            if debug:
                print(f"  Number of sequences: {len(sequences)}")
                for i, (seq_id, seq) in enumerate(list(sequences.items())[:3]):
                    print(f"  Sequence {i+1}: ID={seq_id}, Length={len(seq)}, Sequence: {seq}")
                if len(sequences) > 3:
                    print(f"  ... and {len(sequences)-3} more sequences")
            return len(sequences) > 0
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return False

def is_base_match(query_base, target_base):
    """
    Check if two bases match considering degenerate IUPAC codes
    """
    query_base = query_base.upper()
    target_base = target_base.upper()
    
    # If either base is not in our mapping, assume they don't match
    if query_base not in IUPAC_DEGENERATE or target_base not in IUPAC_DEGENERATE:
        return query_base == target_base
    
    # Get the possible bases for each position
    query_options = IUPAC_DEGENERATE[query_base]
    target_options = IUPAC_DEGENERATE[target_base]
    
    # Check for overlap in possible bases
    return any(base in target_options for base in query_options)

def calculate_degenerate_mismatches(query_seq, target_seq):
    """
    Calculate mismatches between two sequences, considering degenerate bases
    """
    if len(query_seq) != len(target_seq):
        raise ValueError("Sequences must be of equal length for comparison")
    
    mismatches = 0
    for q_base, t_base in zip(query_seq, target_seq):
        if not is_base_match(q_base, t_base):
            mismatches += 1
    
    return mismatches

def align_sequences(query, target, max_mismatch, max_gap, debug=False):
    """
    Align query to target sequence allowing for mismatches and gaps
    Returns a list of (position, mismatch_count, gap_count) tuples
    Supporting degenerate nucleotide codes
    """
    results = []
    query_len = len(query)
    target_len = len(target)
    
    if debug:
        print(f"  Aligning query (len={query_len}) against target (len={target_len})")
        print(f"  Query sequence: {query}")
        print(f"  Target sequence (first 50): {target[:50]}")
    
    # Check for empty sequences
    if query_len == 0 or target_len == 0:
        if debug:
            print("  Empty sequence detected, skipping")
        return results
    
    # Target sequence too short
    if target_len < query_len:
        if debug:
            print("  Target too short for query, skipping")
        return results
    
    # Simple case: Exact match or with mismatches only (no gaps)
    if max_gap == 0:
        for i in range(target_len - query_len + 1):
            target_segment = target[i:i+query_len]
            mismatches = calculate_degenerate_mismatches(query, target_segment)
            
            if mismatches <= max_mismatch:
                results.append((i, mismatches, 0))
                if debug and len(results) <= 3:
                    print(f"  Match at position {i+1} with {mismatches} mismatches")
    else:
        # Basic approach for matches with gaps - simplified for demonstration
        for i in range(target_len - query_len + 1):
            # First try without gaps
            target_segment = target[i:i+query_len]
            mismatches = calculate_degenerate_mismatches(query, target_segment)
            if mismatches <= max_mismatch:
                results.append((i, mismatches, 0))
                if debug and len(results) <= 3:
                    print(f"  Match at position {i+1} with {mismatches} mismatches and 0 gaps")
            
            # Try with gaps, simplified approach
            if i + query_len + max_gap <= target_len:
                extended_segment = target[i:i+query_len+max_gap]
                
                # For each possible gap size
                for gap_size in range(1, max_gap + 1):
                    # For each possible gap starting position
                    for gap_start in range(query_len - 1):
                        # Split the query
                        query_part1 = query[:gap_start]
                        query_part2 = query[gap_start:]
                        
                        # Try different positions for the gap in the target
                        for gap_pos in range(len(query_part1), len(extended_segment) - len(query_part2) + 1):
                            if gap_pos - len(query_part1) > gap_size:
                                continue
                                
                            target_part1 = extended_segment[:gap_pos]
                            target_part2 = extended_segment[gap_pos:]
                            
                            if len(target_part1) >= len(query_part1) and len(target_part2) >= len(query_part2):
                                # Check each part
                                mismatches1 = calculate_degenerate_mismatches(
                                    query_part1, 
                                    target_part1[-len(query_part1):] if len(target_part1) > len(query_part1) else target_part1
                                )
                                
                                mismatches2 = calculate_degenerate_mismatches(
                                    query_part2, 
                                    target_part2[:len(query_part2)]
                                )
                                
                                total_mismatches = mismatches1 + mismatches2
                                
                                if total_mismatches <= max_mismatch:
                                    effective_gap = gap_pos - len(query_part1)
                                    if (i, total_mismatches, effective_gap) not in results:
                                        results.append((i, total_mismatches, effective_gap))
                                        if debug and len(results) <= 3:
                                            print(f"  Match at position {i+1} with {total_mismatches} mismatches and {effective_gap} gaps")
    
    if debug:
        print(f"  Total matches found: {len(results)}")
    
    return results

def search_database(query_file, database_file, max_mismatch, max_gap, include_antisense=False, debug=False):
    """
    Search all query sequences against all database sequences
    Return a dictionary mapping query IDs to lists of matches
    """
    results = defaultdict(list)
    
    # Read query sequences
    if BIOPYTHON_AVAILABLE:
        queries = {}
        for record in SeqIO.parse(query_file, "fasta"):
            queries[record.id] = str(record.seq).upper()
    else:
        queries = {k: v.upper() for k, v in parse_fasta(query_file).items()}
    
    if debug:
        print(f"Loaded {len(queries)} query sequences")
    
    if len(queries) == 0:
        print(f"WARNING: No query sequences found in {query_file}")
        return results
    
    # Read database sequences
    if BIOPYTHON_AVAILABLE:
        database = {}
        for record in SeqIO.parse(database_file, "fasta"):
            database[record.id] = str(record.seq).upper()
    else:
        database = {k: v.upper() for k, v in parse_fasta(database_file).items()}
    
    if debug:
        print(f"Loaded {len(database)} database sequences")
    
    if len(database) == 0:
        print(f"WARNING: No database sequences found in {database_file}")
        return results
    
    # Search each query against each database sequence
    match_count = 0
    for query_id, query_seq in queries.items():
        if debug:
            print(f"Processing query: {query_id}, length: {len(query_seq)}, sequence: {query_seq}")
        
        query_match_count = 0
        for db_id, db_seq in database.items():
            if debug:
                print(f"Comparing against database sequence: {db_id}, length: {len(db_seq)}")
            
            # Search on sense strand
            alignments = align_sequences(query_seq, db_seq, max_mismatch, max_gap, debug)
            
            for position, mismatches, gaps in alignments:
                match_count += 1
                query_match_count += 1
                match_end = min(position+len(query_seq)+gaps, len(db_seq))
                match_seq = db_seq[position:match_end]
                
                results[query_id].append({
                    'database_id': db_id,
                    'position': position + 1,  # 1-based position
                    'mismatches': mismatches,
                    'gaps': gaps,
                    'query_length': len(query_seq),
                    'match_sequence': match_seq,
                    'strand': 'sense'
                })
            
            # If antisense option is enabled, also search on antisense strand
            if include_antisense:
                if debug:
                    print(f"Searching antisense strand for database sequence: {db_id}")
                
                # Generate antisense (reverse complement) of the database sequence
                db_seq_antisense = reverse_complement(db_seq)
                
                if debug:
                    print(f"  Antisense sequence (first 50): {db_seq_antisense[:50]}")
                
                antisense_alignments = align_sequences(query_seq, db_seq_antisense, max_mismatch, max_gap, debug)
                
                for position, mismatches, gaps in antisense_alignments:
                    match_count += 1
                    query_match_count += 1
                    match_end = min(position+len(query_seq)+gaps, len(db_seq_antisense))
                    match_seq = db_seq_antisense[position:match_end]
                    
                    # For antisense matches, position is from the 3' end of the original sequence
                    antisense_position = len(db_seq) - (position + len(match_seq)) + 1
                    
                    results[query_id].append({
                        'database_id': db_id,
                        'position': antisense_position,  # 1-based position from original sequence
                        'mismatches': mismatches,
                        'gaps': gaps,
                        'query_length': len(query_seq),
                        'match_sequence': match_seq,
                        'strand': 'antisense'
                    })
        
        if debug:
            print(f"Found {query_match_count} matches for query {query_id}")
    
    if debug:
        print(f"Total queries with matches: {len(results)}")
        print(f"Total matches found: {match_count}")
    
    return results

def write_results(results, output_file, include_antisense=False):
    """Write search results to output file"""
    with open(output_file, 'w') as out:
        # Write header
        if include_antisense:
            out.write("Query\tDatabase\tPosition\tMismatches\tGaps\tQueryLength\tMatchSequence\tStrand\n")
        else:
            out.write("Query\tDatabase\tPosition\tMismatches\tGaps\tQueryLength\tMatchSequence\n")
        
        # Write results
        match_count = 0
        for query_id, matches in results.items():
            for match in matches:
                match_count += 1
                if include_antisense:
                    out.write(f"{query_id}\t{match['database_id']}\t{match['position']}\t"
                             f"{match['mismatches']}\t{match['gaps']}\t{match['query_length']}\t"
                             f"{match['match_sequence']}\t{match['strand']}\n")
                else:
                    out.write(f"{query_id}\t{match['database_id']}\t{match['position']}\t"
                             f"{match['mismatches']}\t{match['gaps']}\t{match['query_length']}\t"
                             f"{match['match_sequence']}\n")
    
    return match_count

def main():
    args = parse_args()
    debug = args.debug
    include_antisense = args.antisense
    
    try:
        if debug:
            print("Starting exactSearch with parameters:")
            print(f"  Query file: {args.query}")
            print(f"  Database file: {args.database}")
            print(f"  Mismatches allowed: {args.mismatch}")
            print(f"  Gaps allowed: {args.gap}")
            print(f"  Include antisense search: {include_antisense}")
            print(f"  Output file: {args.output}")
            print(f"  Biopython available: {BIOPYTHON_AVAILABLE}")
        
        # Check input files
        if not print_fasta_info(args.query, debug):
            print("Error with query file. Please check format.")
            return 1
        
        if not print_fasta_info(args.database, debug):
            print("Error with database file. Please check format.")
            return 1
        
        # Perform the search
        if debug:
            print("Starting database search with degenerate nucleotide support...")
            if include_antisense:
                print("Antisense (reverse complement) search is enabled")
        
        results = search_database(args.query, args.database, args.mismatch, args.gap, include_antisense, debug)
        
        # Write results to output file
        match_count = write_results(results, args.output, include_antisense)
        
        print(f"Search completed. Found {match_count} matches across {len(results)} queries.")
        if include_antisense:
            print("Results include both sense and antisense matches.")
        print(f"Results written to {args.output}")
        
        if match_count == 0:
            print("WARNING: No matches found. Try increasing the mismatch or gap parameters.")
            if not include_antisense:
                print("You may also try enabling antisense search to look for matches on the complementary strand.")
        
        return 0
        
    except Exception as e:
        sys.stderr.write(f"Error: {str(e)}\n")
        if debug:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
