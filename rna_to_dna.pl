#!/usr/bin/perl -w
# RNA to DNA Converter
# Converts RNA sequences to DNA (U→T) with optional case conversion
# usage: perl rna_to_dna.pl <FASTA input file> <output file> <case conversion: 0=none, 1=upper, 2=lower>

use strict;

# Get command line arguments
my $input_file = $ARGV[0];
my $output_file = $ARGV[1];
my $case_conversion = $ARGV[2] || 0; # Default: no case conversion

# Open input and output files
open (IN, "<$input_file") or die "Cannot open input file $input_file: $!";
open (OUT, ">$output_file") or die "Cannot open output file $output_file: $!";

my $sequence = "";
my $header = "";

while (my $line = <IN>) {
    chomp $line;
    
    # Check if line is a FASTA header
    if ($line =~ /^>/) {
        # If we have a stored sequence from previous header, process it
        if ($sequence ne "" && $header ne "") {
            process_and_write_sequence($sequence, $header);
            $sequence = "";
        }
        
        # Store the new header
        $header = $line;
    } else {
        # Append to the current sequence
        $sequence .= $line;
    }
}

# Process the last sequence
if ($sequence ne "" && $header ne "") {
    process_and_write_sequence($sequence, $header);
}

close(IN);
close(OUT);

sub process_and_write_sequence {
    my ($seq, $head) = @_;
    
    # Convert RNA (U) to DNA (T)
    $seq =~ s/[uU]/T/g;
    
    # Apply case conversion if requested
    if ($case_conversion == 1) {
        # Convert to uppercase
        $seq = uc($seq);
    } elsif ($case_conversion == 2) {
        # Convert to lowercase
        $seq = lc($seq);
    }
    
    # Write to output file
    print OUT "$head\n";
    
    # Print sequence with line wrapping at 60 characters
    while (length($seq) > 60) {
        print OUT substr($seq, 0, 60) . "\n";
        $seq = substr($seq, 60);
    }
    print OUT "$seq\n" if length($seq) > 0;
}
