#!/usr/bin/perl -w

# usage : perl nucleotide_analysis.pl <FASTA file> <output file>
# This script calculates GC%, G%, C%, AT%, A%, T% and counts of all nucleotides for each sequence in a FASTA file

open (IN, "<$ARGV[0]") or die "Cannot open input file: $!";
open (OUT, ">$ARGV[1]") or die "Cannot open output file: $!";

# Print header
print OUT "Sequence_ID\tTotal_Count\tGC%\tAT%\tG%\tC%\tA%\tT%\tG_Count\tC_Count\tA_Count\tT_Count\n";

my $seq_id = "";
my $g_count = 0;     # G count
my $c_count = 0;     # C count
my $a_count = 0;     # A count
my $t_count = 0;     # T count (also counts U for RNA)
my $length = 0;      # Sequence length

while (<IN>) {
    chomp;
    if (m/^>/) {
        # Process previous sequence if not the first header
        if ($length > 0) {
            my $gc_count = $g_count + $c_count;
            my $at_count = $a_count + $t_count;
            
            my $gc_percent = ($gc_count / $length) * 100;
            my $at_percent = ($at_count / $length) * 100;
            my $g_percent = ($g_count / $length) * 100;
            my $c_percent = ($c_count / $length) * 100;
            my $a_percent = ($a_count / $length) * 100;
            my $t_percent = ($t_count / $length) * 100;
            
            print OUT "$seq_id\t";
            print OUT "$length\t";
            print OUT sprintf("%.3f", $gc_percent) . "\t";
            print OUT sprintf("%.3f", $at_percent) . "\t";
            print OUT sprintf("%.3f", $g_percent) . "\t";
            print OUT sprintf("%.3f", $c_percent) . "\t";
            print OUT sprintf("%.3f", $a_percent) . "\t";
            print OUT sprintf("%.3f", $t_percent) . "\t";
            print OUT "$g_count\t$c_count\t$a_count\t$t_count\n";
        }
        
        # Get new sequence ID (remove '>' and get first word)
        s/^>//;
        $seq_id = (split)[0];
        
        # Reset counters for new sequence
        $g_count = 0;
        $c_count = 0;
        $a_count = 0;
        $t_count = 0;
        $length = 0;
    } else {
        # Convert to uppercase for case-insensitive matching
        my $uppercase = uc($_);
        
        # Count occurrences of each nucleotide
        my $g_line_count = ($uppercase =~ tr/G//);
        my $c_line_count = ($uppercase =~ tr/C//);
        my $a_line_count = ($uppercase =~ tr/A//);
        my $t_line_count = ($uppercase =~ tr/TU//); # Count both T and U (for RNA)
        
        # Update counters
        $g_count += $g_line_count;
        $c_count += $c_line_count;
        $a_count += $a_line_count;
        $t_count += $t_line_count;
        
        # Update sequence length
        $length += length $_;
    }
}

# Process the last sequence
if ($length > 0) {
    my $gc_count = $g_count + $c_count;
    my $at_count = $a_count + $t_count;
    
    my $gc_percent = ($gc_count / $length) * 100;
    my $at_percent = ($at_count / $length) * 100;
    my $g_percent = ($g_count / $length) * 100;
    my $c_percent = ($c_count / $length) * 100;
    my $a_percent = ($a_count / $length) * 100;
    my $t_percent = ($t_count / $length) * 100;
    
    print OUT "$seq_id\t";
    print OUT "$length\t";
    print OUT sprintf("%.3f", $gc_percent) . "\t";
    print OUT sprintf("%.3f", $at_percent) . "\t";
    print OUT sprintf("%.3f", $g_percent) . "\t";
    print OUT sprintf("%.3f", $c_percent) . "\t";
    print OUT sprintf("%.3f", $a_percent) . "\t";
    print OUT sprintf("%.3f", $t_percent) . "\t";
    print OUT "$g_count\t$c_count\t$a_count\t$t_count\n";
}

close(IN);
close(OUT);
