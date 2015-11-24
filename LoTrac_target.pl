#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;
use File::Basename;

###MODULE LOAD###
#module load perl/5.12.3
#module load ncbi-blast+/2.2.29
#module load BEDTools/2.17.0
#module load Python/2.7
#module load prodigal/2.60

sub checkOptions {
    my %opts;
    getopts('h1:2:q:o:n:L:I:S:f', \%opts);
    my ($help, $fastq1, $fastq2, $query, $outDir, $outName, $length, $identity, $gSize, $frag);

    if($opts{h}) {
        $help = $opts{h};
        help();
    }

    if($opts{1}) {
        $fastq1 = $opts{1};
        if( -e $fastq1) {
            print "Paired-end Read 1 is: $fastq1\n";
        } else {
            print "The forward paired-end file name is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "No paired end 1 fastq file path argument given.\n";
        help();
    }

    if($opts{2}) {
        $fastq2 = $opts{2};
        if( -e $fastq2) {
            print "Paired-end Read 2 is: $fastq2\n";
        } else {
            print "The reverse paired-end file name is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "No paired end 2 fastq file path argument given.\n";
        help();
    }

    if($opts{q}) {
        $query = $opts{q};
        if (-e $query) {
            print "File containing the query reference sequence: $query\n";
        } else {
            print "The location given for the query reference sequence is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/query_file).\n";
            help();
        }
    } else {
        print "The location of the query reference sequence (including full path) has not been given.\n";
        help();
    }

    $outDir = "./";
    if($opts{o}) {
        if (-d $opts{o}) {
            $outDir = $opts{o};
            print "The output directory is: $outDir\n";
        } else {
            $outDir = $opts{o};
            mkdir $outDir;
            print "The output directory has been created: $outDir\n";
        }
    } else {
        print "The files will be output into the current directory.\n";
    }

    if($opts{n}) {
        $outName = $opts{n};
        print "The output file name prefix: $outName\n";
    } else {
        $outName = `echo "$fastq1" | awk -F"/" '{print \$(NF)}' | sed 's/_S[0-9]\\+_L[0-9]\\+_R[0-9]\\+.*//g'`;
        print "The default output file name prefix is: $outName";
    }

    if($opts{S}) {
        $gSize = $opts{S};
	if ($gSize =~ /[0-9]+[k|M|G]$/ || $gSize =~ /[0-9]+$/) {
	    print "The given genome size is: $gSize\n";
	} else {
	    print "The genome size is not in the proper format.\n";
	    print "Please provide the genome size in bp (can use k/M/G suffix).\n";
	    help();
        }
    } else {
        print "The genome size argument has not been given\n";
        help();
    }    

    $length = 0.5;
    if($opts{L}) {
	$length = $opts{L};
	print "The alignment length threshold: $length\n";
    } else {
	print "The default length threshold of 0.5 (50%) will be used\n";
    }

    $identity = 0.5;
    if($opts{I}) {
        $length = $opts{I};
        print "The alignment identity threshold: $identity\n";
    } else {
        print "The default identity threshold of 0.5 (50%) will be used\n";
    }
    
    if($opts{f}) {
	$frag = "yes";
	print "The extract fragment flag has been given\n";
    }

    #($help, $fastq1, $fastq2, $query, $outDir, $outName, $length, $identity, $gSize, $frag)
    return ($help, $fastq1, $fastq2, $query, $outDir, $outName, $length, $identity, $gSize, $frag);
}

sub help
{

die <<EOF

USAGE
emm_typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -q <query sequence file: file path> -o <output directory name: string> -n <output name prefix: string> -S <genome size>  [OPTIONS]
              
    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -q   query reference sequence file (including full path)
    -o   output directory
    -n   output name prefix
    -L   alignment length threshold (default is 0.5 (50%))
    -I   alignment identity threshold (default is 0.5 (50%))
    -S   organism genome size in bp (can use k/M/G suffix)
    -f   extract query fragment flag

EOF
}

my ($help, $fastq1, $fastq2, $query, $outDir, $outName, $length, $identity, $gSize, $frag) = checkOptions( @ARGV );




###Subroutines###

sub fasta_seq_length {
    my ($seq) = @_;
    #open ( my $q_seq, "<", $seq ) or die "Could not open file '$seq': $!";    
    my @lines = split /\n/, $seq;
    my $final_line;
    foreach my $line (@lines) {
	chomp($line);
	if ($line =~ /^>/) {
	    next;
	} else {
	    #print "line: $line\n";
	    $final_line = $final_line.$line;
        }
    }
    return length($final_line);
}

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




##Start Doing Stuff##
chdir "$outDir";
###Preprocess with Cutadapt###
my $fastq1_trimd = "cutadapt_".$outName."_S1_L001_R1_001.fastq";
my $fastq2_trimd = "cutadapt_".$outName."_S1_L001_R2_001.fastq";
if( -e $fastq1_trimd) {
    print "Fastq files have already been preprocessed\n";
} else {
    print "Beginning cutadapt\n";
    system("cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq $fastq1 $fastq2");
    system("cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output $fastq1_trimd -o $fastq2_trimd temp2.fastq temp1.fastq");
    my $tempDel_1 = "temp1.fastq";
    my $tempDel_2 = "temp2.fastq";
    unlink $tempDel_1;
    unlink $tempDel_2;
}

if( -d "./velvet_output") {
    print "Velvet assembly has already been completed\n";
} else {
    print "Beginning Velvet\n";
    my $velvetK_val = `velvetk.pl --best --size "$gSize" "$fastq1_trimd" "$fastq2_trimd"`;
    `VelvetOptimiser.pl -s "$velvetK_val" -e "$velvetK_val" -o "-scaffolding no" -f "-shortPaired -separate -fastq $fastq1_trimd $fastq2_trimd" -d velvet_output`;
}

print "Beginning Prodigal\n";
if (glob("prodigal_$outName*")) {
    print "Gene prediction has already been completed\n";
} else {
    system("prodigal -c -f gff -i ./velvet_output/contigs.fa -a PRE_$outName.faa -o prodigal_$outName.gff -d PRE_$outName.fasta");
    `cat PRE_"$outName".faa | sed 's/ # .*//g' > prodigal_"$outName".faa`;
    `cat PRE_"$outName".fasta | sed 's/ # .*//g' > prodigal_"$outName".fna`;
    unlink("PRE_$outName.faa");
    unlink("PRE_$outName.fasta");
}

#print "Create a blast database using the sequences in the query fasta file\n";
#(my $query_headr = $query) =~ s/\.f.*$//g;
#my $query_dir = dirname($query);
#my $blast_DB = "TEMP_".$query_headr."_nucl_blast_db";
#if (glob("$blast_DB*")) {
#    print "Gene blast database has already been created\n";
#} else {
#    system("makeblastdb -in $query -dbtype nucl -out $query_dir/$blast_DB");
#}

print "Create a blast database using the predicted genes obtained from Prodigal\n";
if (glob("TEMP_prod_nucl_blast_db*")) {
    print "Contig blast database has already been created\n";
} else {
    system("makeblastdb -in prodigal_$outName.fna -dbtype nucl -out TEMP_prod_nucl_blast_db");
}

print "Create a blast database using the assembled contigs obtained from Velvet\n";
if (glob("TEMP_velvet_nucl_blast_db*")) {
    print "Contig blast database has already been created\n";
} else {
    system("makeblastdb -in ./velvet_output/contigs.fa -dbtype nucl -out TEMP_velvet_nucl_blast_db");
}

###Blast each sequence given in the query fasta file against the blast nucleotide database.###
my @query_names;
open ( my $q_seq, "<", $query ) or die "Could not open file '$query': $!";
while ( my $line = <$q_seq> ) {
    if ($line =~ />.*/) {
        $line =~ s/>//g;
        chomp($line);
	push(@query_names,$line);
    }
}
close $q_seq;
#print "query names:\n$query_names[0]";

foreach (@query_names) {
    my $query_name = $_;
    my $extract_out = "EXTRACT_".$query_name."_target.fasta";
    my $error_out = "ERROR_".$query_name."_target.fasta";

    ###OPEN 'EXTRACT_$query_name_target.fasta' FOR APPENDING (CHECK FOR FAILURES)
    open ( my $exOUT, ">>", $extract_out ) or die "Could not open file $extract_out: $!";
    print $exOUT "Complete $query_name Gene Sequence:\n";

    my $query_seq = extractFastaByID($query_name,$query);
    my $query_length = fasta_seq_length($query_seq);
    open ( my $qOUT, ">", 'TEMP_query_sequence.fna' ) or die "Could not open file TEMP_query_sequence.fna: $!";
    print $qOUT $query_seq;
    close $qOUT;
    system("blastn -db TEMP_prod_nucl_blast_db -query TEMP_query_sequence.fna -outfmt 6 -word_size 7 -out TEMP_prod-vs-query_blast.txt");

    ###Get the best blast hit by sorting the blast output by bit score, then % ID, then alignment length and select the first hit as the best match.###
    my $bestHit = `cat TEMP_prod-vs-query_blast.txt | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
    my @bestArray = split('\t',$bestHit);
    my $best_name = $bestArray[1];
    my $best_iden = $bestArray[2];
    my $best_len = $bestArray[3];
    #my $query_strt = $bestArray[6];
    #my $query_end = $bestArray[7];
    my $match_len = $length * $query_length;

    print "\nprodigal gene name of best hit against the query sequence: $best_name\n";
    print "% identity of best hit against the query sequence: $best_iden\n";
    print "length of best hit against the query sequence: $best_len\n";

    if ($best_iden >= $identity && $best_len >= $match_len) {
	my $prodigal_fna = `extractFastaByID.pl $best_name < prodigal_"$outName".fna`;
	my $prodigal_faa = `extractFastaByID.pl $best_name < prodigal_"$outName".faa`;
	print $exOUT "$prodigal_fna\n$prodigal_faa\n\n";
    } else {
	open ( my $errOUT, ">>", $error_out ) or die "Could not open file $error_out: $!";
	print $errOUT "The best blast hit ($best_name) for $query_name didn't meet minimum criteria of length and identity to call a true match\n\n";
	close $errOUT;
	#next;
    }

    ###If the '-f' flag is called, then the script will also extract just the section of the target sequence that corresponds to the query fragment.###
    ###If the best blast hit didn't include the entire query fragment, then the code will calculate the expected start/end coordinates of the complete fragment
    ###and will attempt to extract the full fragment from the predicted gene sequence.###
    ###The number of non-aligning bases at each end of the matching target sequence will recorded in the header name.###
    if ($frag) {
	print "Extracting target fragment\n";
	#my $velvet_blast = "velvet-vs-query_".$query_name."_blast.txt";
	system("blastn -db TEMP_velvet_nucl_blast_db -query TEMP_query_sequence.fna -outfmt 6 -word_size 7 -out TEMP_velvet-vs-query_blast.txt");
	my $bestHit = `cat TEMP_velvet-vs-query_blast.txt | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
	my @bestArray = split('\t',$bestHit);
	my $best_name = $bestArray[1];
	my $best_iden = $bestArray[2];
	my $best_len = $bestArray[3];
	my $query_strt = $bestArray[6];
	my $query_end = $bestArray[7];
	my $frag_length = $best_len / $query_length;
	#print "best hit: $bestHit || $frag_length\n";
	
	print "\ncontig name of best hit against the query sequence: $best_name\n";
	print "% identity of best hit against the query sequence: $best_iden\n";
	print "length of best hit against the query sequence: $best_len\n";
	
	if ($best_iden >= 50 && $frag_length >= 0.50) {
	    print $exOUT "\n$query_name Query Fragment Sequence:\n";
	    if ($bestArray[9] > $bestArray[8]) {
		#my $frag_start = $bestArray[8] - 1;
		my $blast_endDiff = $query_length - $bestArray[7];
	        my $frag_start = $bestArray[8] - $bestArray[6];
		my $frag_end = $blast_endDiff + $bestArray[9];
		open(my $fh, '>', 'TEMP_frwd_extract.bed') or die "Could not open file 'TEMP_frwd_extract.bed' $!";
		print $fh "$best_name\t$frag_start\t$frag_end\n";
		close $fh;
		my $extract_frag_frwd = `bedtools getfasta -fi ./velvet_output/contigs.fa -bed TEMP_frwd_extract.bed -fo stdout`;
		print $exOUT "$extract_frag_frwd\n";
	    } elsif ($bestArray[9] < $bestArray[8]) {
		#my $query_extract = $query_strt - 500;
		my $blast_endDiff = $query_length - $bestArray[7];
		my $frag_start = $bestArray[9] - $bestArray[6];
		my $frag_end = $blast_endDiff + $bestArray[8];
		open(my $fh, '>', 'TEMP_rev_extract.bed');
		print $fh "$best_name\t$frag_start\t$frag_end\n";
		close $fh;

		my $extract_frag_rev = `bedtools getfasta -tab -fi ./velvet_output/contigs.fa -bed TEMP_rev_extract.bed -fo stdout`;
		print "extract frag is:\n$extract_frag_rev\n";
		my @rev_frag_array = split('\t',$extract_frag_rev);
		my $rev_comp_frag = reverse($rev_frag_array[1]);
		$rev_comp_frag =~ tr/ATGCatgc/TACGtacg/;

		print $exOUT ">$rev_frag_array[0]";
		print $exOUT "$rev_comp_frag\n";
	    }
	} else {
	    open ( my $errOUT, ">>", $error_out ) or die "Could not open file $error_out: $!";
	    print $errOUT "The best blast hit ($best_name) for $query_name fragment didn't meet minimum criteria of length and identity to call a true match\n\n";
	    close $errOUT;
	    next;
	}	
    }
close $exOUT;
#close $errOUT;
}

#rm TEMP*
