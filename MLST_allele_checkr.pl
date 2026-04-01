#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($MLST_input,$bam_input,$MLST_ref) = @ARGV;
open (MY_MLST_INPUT, "$MLST_input");
my @MLST_results;
my $input_counter = 1;
while (<MY_MLST_INPUT>) {
    my $line = $_;
    chomp($line);
    if ($input_counter == 1) {
	$input_counter++;
    } elsif ($input_counter == 2) {
	@MLST_results = split('\t',$line);
	last;
    }
}

my $MLST_depth = $MLST_results[-2];
my $MLST_mismatch = $MLST_results[-4];
print "MLST_depth is: $MLST_depth\n";
print "MLST_mismatches: $MLST_mismatch\n\n";



###Subroutines###
sub extractFastaByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
    #print "lookup: $lookup\n";
    local $/ = "\n>";  # read by FASTA record

    my $output;
    while (my $seq = <$fh>) {
        chomp $seq;
        #print "while seq:\n$seq\n";
        my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        if ($id eq $lookup) {
            $seq =~ s/^>*.+\n//;  # remove FASTA header
            #$seq =~ s/\n//g;  # remove endlines
            #print ">$id\n";
            #print "$seq\n";
            $output = ">$id\n$seq\n";
            last;
        }
    }
    return $output;
}



###Start Doing Stuff###
my @mismatch_array = split(';',$MLST_mismatch);
if (! $MLST_mismatch == 0) {
    open ( my $exFile_out, ">>", 'Check_Target_Sequence.txt' ) or die "Could not open file 'Check_Target_Sequence.txt': $!";
    print $exFile_out '#' x 65;
    print $exFile_out "--NEW MLST SEQUENCE--";
    print $exFile_out '#' x 65;
    print $exFile_out "\n\n";
    if ($MLST_depth >= 30) {
	print "Average MLST depth is above threshold.\n";
        system("samtools index $bam_input");
	foreach (@mismatch_array) {
	    print $exFile_out "For MLST Allele/SNP: $_\n";
	    $_ =~ /(.*)\/.*/;
	    my $extract_allele = $1;
	    my $pileup_allele = `samtools mpileup -f $MLST_ref $bam_input -r $extract_allele`;

	    my $MLST_bam = $extract_allele."_".$bam_input;
	    (my $MLST_bai = $MLST_bam) =~ s/\.bam/\.bai/g;
	    (my $MLST_vcf = $MLST_bam) =~ s/\.bam/\.vcf/g;
	    (my $MLST_fna = $MLST_bam) =~ s/\.bam/\.fna/g;
	    system("samtools view -b $bam_input $extract_allele > $MLST_bam");
	    #system("samtools view -h $MLST_bam");
	    system("samtools index $MLST_bam $MLST_bai");
	    open ( my $MLST_out, ">", $MLST_fna ) or die "Could not open file $MLST_fna: $!";
	    my $allele_ref = extractFastaByID($extract_allele,$MLST_ref);
            print $MLST_out "$allele_ref";

            ###Create the variant-called consensus fasta###
	    system("freebayes -q 20 -p 1 -f $MLST_fna $MLST_bam -v $MLST_vcf");
	    system("bgzip $MLST_vcf");
	    system("tabix -p vcf $MLST_vcf.gz");
	    my $MLST_consensus = `cat $MLST_fna | vcf-consensus $MLST_vcf.gz`;
	    #print "MLST consensus:\n$MLST_consensus\n\n";
	    #print "MLST pileup:\n$pileup_allele\n";
	    print $exFile_out "New MLST Allele Consensus:\n$MLST_consensus\n";
	    print $exFile_out "New MLST Allele Pileup:\n$pileup_allele\n\n";
	}
    } else {
	print $exFile_out "The new allele should not be trusted because the average depth was not above the minimum threshold\n";
    }
    print $exFile_out '#' x 151;
    print $exFile_out "\n\n";
    close $exFile_out;
} else {
    print "No new MLST alleles\n";
}

#unlink("$MLST_bam");
