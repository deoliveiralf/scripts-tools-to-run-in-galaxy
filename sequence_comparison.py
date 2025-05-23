#!/usr/bin/env python3

import argparse
import csv
import sys

def compare_sequences(input_file1, input_file2, output_file, 
                      col1_id_total, col1_length_total, 
                      col2_id_aligned, col2_length_aligned,
                      header1, header2,
                      additional_column_index,
                      only_full_alignment):
    """
    Compare sequence lengths between two input files and calculate alignment percentage.
    
    Args:
    input_file1 (str): Path to the first input file with sequence information
    input_file2 (str): Path to the second input file with sequence information
    output_file (str): Path to the output file
    col1_id_total (int): Column index for sequence ID in first file (0-based)
    col1_length_total (int): Column index for total length in first file (0-based)
    col2_id_aligned (int): Column index for sequence ID in second file (0-based)
    col2_length_aligned (int): Column index for aligned length in second file (0-based)
    header1 (bool): Whether the first file has a header line
    header2 (bool): Whether the second file has a header line
    additional_column_index (int): Column index for additional information from second file
    only_full_alignment (bool): Whether to output only 100% aligned sequences
    """
    # Dictionaries to store sequence information
    total_sequences = {}
    second_input_data = []

    # Read first input file (total sequences)
    try:
        with open(input_file1, 'r') as f1:
            reader1 = csv.reader(f1, delimiter='\t')
            
            # Skip header if present
            if header1:
                next(reader1)
            
            for row in reader1:
                if len(row) > max(col1_id_total, col1_length_total):
                    seq_id = row[col1_id_total]
                    try:
                        total_sequences[seq_id] = int(row[col1_length_total])
                    except ValueError:
                        print(f"Warning: Invalid length value for ID {seq_id} in first file. Skipping.", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: First input file {input_file1} not found.", file=sys.stderr)
        sys.exit(1)

    # Read second input file
    try:
        with open(input_file2, 'r') as f2:
            reader2 = csv.reader(f2, delimiter='\t')
            
            # Skip header if present
            if header2:
                next(reader2)
            
            for row in reader2:
                if len(row) > max(col2_id_aligned, col2_length_aligned, additional_column_index):
                    seq_id = row[col2_id_aligned]
                    
                    # Check if the sequence ID exists in total_sequences
                    if seq_id in total_sequences:
                        try:
                            aligned_length = int(row[col2_length_aligned])
                            additional_info = row[additional_column_index]
                            
                            # Calculate alignment percentage
                            total_length = total_sequences[seq_id]
                            alignment_percentage = (aligned_length / total_length) * 100
                            
                            # Apply full alignment filter if requested
                            if not only_full_alignment or alignment_percentage == 100:
                                second_input_data.append({
                                    'id': seq_id,
                                    'total_length': total_length,
                                    'aligned_length': aligned_length,
                                    'alignment_percentage': alignment_percentage,
                                    'additional_info': additional_info
                                })
                        except (ValueError, IndexError) as e:
                            print(f"Warning: Error processing row {row}: {e}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Second input file {input_file2} not found.", file=sys.stderr)
        sys.exit(1)

    # Write output
    try:
        with open(output_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            
            # Write header
            header = ['ID', 'Total_Nucleotides', 'Aligned_Nucleotides', 'Alignment_Percentage', 'Additional_Column']
            writer.writerow(header)

            # Write rows
            for entry in second_input_data:
                writer.writerow([
                    entry['id'], 
                    entry['total_length'], 
                    entry['aligned_length'], 
                    f"{entry['alignment_percentage']:.2f}", 
                    entry['additional_info']
                ])

        print(f"Comparison complete. Results written to {output_file}", file=sys.stderr)

    except IOError:
        print(f"Error: Unable to write to output file {output_file}", file=sys.stderr)
        sys.exit(1)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Compare sequence lengths and calculate alignment percentage.')
    parser.add_argument('input1', help='First input file with sequence IDs and total lengths')
    parser.add_argument('input2', help='Second input file with sequence IDs and additional information')
    parser.add_argument('output', help='Output file for comparison results')
    parser.add_argument('--col1_id_total', type=int, default=0, help='Column index for sequence ID in first file (0-based)')
    parser.add_argument('--col1_length_total', type=int, default=1, help='Column index for total length in first file (0-based)')
    parser.add_argument('--col2_id_aligned', type=int, default=0, help='Column index for sequence ID in second file (0-based)')
    parser.add_argument('--col2_length_aligned', type=int, default=1, help='Column index for aligned length in second file (0-based)')
    parser.add_argument('--header1', action='store_true', help='First file contains a header line')
    parser.add_argument('--header2', action='store_true', help='Second file contains a header line')
    parser.add_argument('--additional_column_index', type=int, default=1, help='Column index for additional information from second file')
    parser.add_argument('--only_full_alignment', action='store_true', help='Output only 100% aligned sequences')

    # Parse arguments
    args = parser.parse_args()

    # Run comparison
    compare_sequences(
        args.input1, 
        args.input2, 
        args.output, 
        args.col1_id_total,
        args.col1_length_total,
        args.col2_id_aligned,
        args.col2_length_aligned,
        args.header1,
        args.header2,
        args.additional_column_index,
        args.only_full_alignment
    )

if __name__ == '__main__':
    main()
