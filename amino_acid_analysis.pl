#!/usr/bin/perl -w
# Galaxy wrapper for amino acid analysis
# This script calculates percentages and counts of all amino acids for each protein sequence in a FASTA file

use strict;
use warnings;

# Simple Galaxy input/output handling
open (IN, "<$ARGV[0]") or die "Cannot open input file: $!";
open (OUT, ">$ARGV[1]") or die "Cannot open output file: $!";

# Define the standard amino acids
my @amino_acids = qw(A R N D C Q E G H I L K M F P S T W Y V);

# Print header
print OUT "Sequence_ID\tTotal_Length";
foreach my $aa (@amino_acids) {
    print OUT "\t${aa}%";
}
foreach my $aa (@amino_acids) {
    print OUT "\t${aa}_Count";
}
print OUT "\n";

my $seq_id = "";
my %aa_counts = map { $_ => 0 } @amino_acids;
my $length = 0;
my $sequence = "";

while (<IN>) {
    chomp;
    if (m/^>/) {
        # Process previous sequence if not the first header
        if ($sequence) {
            # Process the complete sequence
            $length = length($sequence);
            count_amino_acids($sequence, \%aa_counts);
            process_and_output_sequence($seq_id, $length, \%aa_counts);
        }
        
        # Get new sequence ID (remove '>' and get first word)
        s/^>//;
        $seq_id = (split)[0];
        
        # Reset for new sequence
        %aa_counts = map { $_ => 0 } @amino_acids;
        $sequence = "";
    } else {
        # Add this line to the sequence
        $sequence .= uc($_);
    }
}

# Process the last sequence
if ($sequence) {
    $length = length($sequence);
    count_amino_acids($sequence, \%aa_counts);
    process_and_output_sequence($seq_id, $length, \%aa_counts);
}

close(IN);
close(OUT);

# Subroutine to count amino acids in a sequence
sub count_amino_acids {
    my ($seq, $counts_ref) = @_;
    foreach my $aa (@amino_acids) {
        # Count occurrences of each amino acid
        my $count = () = $seq =~ /$aa/g;
        $counts_ref->{$aa} = $count;
    }
}

# Subroutine to calculate percentages and output results
sub process_and_output_sequence {
    my ($id, $len, $counts_ref) = @_;
    
    print OUT "$id\t$len";
    
    # Print percentages
    foreach my $aa (@amino_acids) {
        my $percent = ($len > 0) ? ($counts_ref->{$aa} / $len) * 100 : 0;
        print OUT "\t" . sprintf("%.3f", $percent);
    }
    
    # Print counts
    foreach my $aa (@amino_acids) {
        print OUT "\t" . $counts_ref->{$aa};
    }
    
    print OUT "\n";
}
