#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;

###MODULE LOAD###
#module load perl/5.12.3
#module load ncbi-blast+/2.2.29
#module load BEDTools/2.17.0
#module load Python/2.7

sub checkOptions {
    my %opts;
    getopts('h1:2:r:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $emm_DB, $outDir, $outName);

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

    if($opts{r}) {
        $emm_DB = $opts{r};
        if (-d $emm_DB) {
            print "Directory containing the emm reference sequence: $emm_DB\n";
        } else {
            print "The emm reference sequence directory location is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The emm reference sequence directory (including full path) has not been given.\n";
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

    return ($help, $fastq1, $fastq2, $emm_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
emm_typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <reference databases directory: file path> -o <output directory name: string> -n <output name prefix: string>  [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   emm reference sequence directory (including full path)
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $emm_DB, $outDir, $outName) = checkOptions( @ARGV );


##Start Doing Stuff##
chdir "$outDir";
my $emmType_output = "emm-Type_Results.txt";
open(my $fh,'>',$emmType_output) or die "Could not open file '$emmType_output' $!";
print $fh "Sample_Name\temm_Type\temm_Seq\tPercent_Identity\tMatch_Length\n";

###Preprocess with Cutadapt###
my $fastq1_trimd = "cutadapt_".$outName."_S1_L001_R1_001.fastq";
my $fastq2_trimd = "cutadapt_".$outName."_S1_L001_R2_001.fastq";
if( -e $fastq1_trimd) {
    print "Fastq files have already been preprocessed\n";
} else {
    print "Beginning cutadapt\n";
    system("cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq $fastq1 $fastq2");
    system("cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output $fastq1_trimd -o $fastq2_trimd temp2.fastq temp1.fastq");
    #`cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq $fastq1 $fastq2`;
    #`cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output $fastq1_trimd -o $fastq2_trimd temp2.fastq temp1.fastq`;
    #system("cutadapt --version > cut_temp.txt")
    #system("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 -o $fastq1_trimd -p $fastq2_trimd $fastq1 $fastq2");
    my $tempDel_1 = "temp1.fastq";
    my $tempDel_2 = "temp2.fastq";
    unlink $tempDel_1;
    unlink $tempDel_2;
}

if( -d "./velvet_output") {
    print "Velvet assembly has already been completed\n";
} else {
    print "Beginning Velvet\n";
    my $velvetK_val = `velvetk.pl --best --size 1.8M "$fastq1_trimd" "$fastq2_trimd"`;
    `VelvetOptimiser.pl -s "$velvetK_val" -e "$velvetK_val" -o "-scaffolding no" -f "-shortPaired -separate -fastq $fastq1_trimd $fastq2_trimd" -d velvet_output`;
}

print "Blast assembled contigs against the forward primer reference sequence\n";
if (glob("$emm_DB/blast_frwd_primr-nucl_DB*")) {
    system("blastn -db $emm_DB/blast_frwd_primr-nucl_DB -query ./velvet_output/contigs.fa -outfmt 6 -word_size 4 -out contig-vs-frwd_nucl.txt");
} else {
    system("makeblastdb -in $emm_DB/frwd_primr-DB_Final.fasta -dbtype nucl -out $emm_DB/blast_frwd_primr-nucl_DB");
    system("blastn -db $emm_DB/blast_frwd_primr-nucl_DB -query ./velvet_output/contigs.fa -outfmt 6 -word_size 4 -out contig-vs-frwd_nucl.txt");
}

###Get the best blast hit by sorting the blast output by bit score, then % ID, then alignment length and select the first hit as the best match.###
my $frwd_bestHit = `cat contig-vs-frwd_nucl.txt | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
my @frwd_bestArray = split('\t',$frwd_bestHit);
my $best_frwd_name = $frwd_bestArray[0];
my $best_frwd_len = $frwd_bestArray[3];
my $best_frwd_iden = $frwd_bestArray[2];
my $query_strt = $frwd_bestArray[6];
my $query_end = $frwd_bestArray[7];
#my $target_strt = $frwd_bestArray[8];
#my $target_end = $frwd_bestArray[9];

print "\nname of best hit against the emm forward primers: $best_frwd_name\n";
print "% identity of best hit against the emm foward primers: $best_frwd_iden\n";
print "length of best hit against the emm forward primer: $best_frwd_len\n";

if ($best_frwd_len == 19 && $best_frwd_iden >= 94.5) {
    if ($frwd_bestArray[9] > $frwd_bestArray[8]) {
	##Extract the +500 sequence with Bedtools##
	my $query_extract = $query_strt + 500;
	open(my $fh, '>', 'emm_region_extract.bed') or die "Could not open file 'emm_region_extract.bed' $!";
	print $fh "$best_frwd_name\t$query_strt\t$query_extract\n";
	close $fh;
	#print "done\n";
	system("bedtools getfasta -fi ./velvet_output/contigs.fa -bed emm_region_extract.bed -fo emm_region_extract.fasta");
    } else {
	print "extract from reverse strand\n";
	##Extract the -500 sequence with Bedtools##
        my $query_extract = $query_strt - 500;
        open(my $fh, '>', 'emm_region_extract.bed');
        print $fh "$best_frwd_name\t$query_extract\t$query_strt\n";
        close $fh;
	
	my $extract_emm1 = `bedtools getfasta -tab -fi ./velvet_output/contigs.fa -bed emm_region_extract.bed -fo stdout`;
        print "extract emm is:\n$extract_emm1\n";
	my @emm1_array = split('\t',$extract_emm1);
	my $rev_comp_emm1 = reverse($emm1_array[1]);
	$rev_comp_emm1 =~ tr/ATGCatgc/TACGtacg/;

        open(my $fh2, '>', 'emm_region_extract.fasta') or die "Could not open file 'emm_region_extract.fasta' $!";
        print $fh2 ">$emm1_array[0]";
	print $fh2 "$rev_comp_emm1\n";
        close $fh2;
    }
} else {
    print "The best blast hit ($best_frwd_name) obtained from querying the assembled contigs against the emm forward primers\ndidn't meet minimum criteria of length and identity to call a true match.\n";
    print $fh "$outName\tExtraction_Error\t--\t--\t--\n";
    my $old_name = "emm-Type_Results.txt";
    my $emm_out = $outName."__emm-Type__Results.txt";
    rename $old_name, $emm_out;
    exit
}

###Blast extracted emm region against database to find match###
if ( -s "emm_region_extract.fasta") {
    if (glob("$emm_DB/blast_emm_Gene-nucl_DB*")) {
	system("blastn -db $emm_DB/blast_emm_Gene-nucl_DB -query emm_region_extract.fasta -outfmt 6 -word_size 4 -out emm_vs_DB_nucl.txt");
    } else {
	system("makeblastdb -in $emm_DB/emm_Gene-DB_Final.fasta -dbtype nucl -out $emm_DB/blast_emm_Gene-nucl_DB");
	system("blastn", "-db", "$emm_DB/blast_emm_Gene-nucl_DB", "-query", "emm_region_extract.fasta", "-outfmt", "6", "-word_size", "4", "-out", "emm_vs_DB_nucl.txt");
    }
} else {
    print "Although the best blast hit found a true match against the 19bp primer, the matching contig didn't contain the\nfull 500bp region downstream of the primer that comprises the emm-typing region\n";
    print $fh "$outName\tExtraction_Error\t--\t--\t--\n";
    my $old_name = "emm-Type_Results.txt";
    my $emm_out = $outName."__emm-Type__Results.txt";
    rename $old_name, $emm_out;
    exit
}

my $emm_bestHit = `cat emm_vs_DB_nucl.txt | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
my @emm_bestArray = split('\t',$emm_bestHit);
my $best_emm_name = $emm_bestArray[1];
my $best_emm_len = $emm_bestArray[3];
my $best_emm_iden = $emm_bestArray[2];

#print "\nname of best hit in the emm database: $emm_bestArray[1]\n";
#print "identity of best hit in the emm database: $emm_bestArray[2]\n";
#print "length of best hit in the emm database: $emm_bestArray[3]\n";

$outName =~ /EMM_(.*)/;
my $finalName = $1;
if ($best_emm_iden == 100 && $best_emm_len == 180) {
    #$best_emm_name =~ /\d+__EMM(.*)__EMM.*__\d+/;
    $best_emm_name =~ /\d+__[A-Z]+(.*)__.*__\d+/;
    my $emmType = $1;
    print $fh "$finalName\t$emmType\t$best_emm_name\t$best_emm_iden\t$best_emm_len\n";
} else {
    #$best_emm_name =~ /\d+__EMM(.*)__EMM.*__\d+/;
    $best_emm_name =~ /\d+__[A-Z]+(.*)__.*__\d+/;
    my $emmType = $1;
    print $fh "$finalName\t$emmType*\t$best_emm_name\t$best_emm_iden\t$best_emm_len\n";
    ###OPEN 'Check_Target_Sequence.txt' FOR APPENDING (CHECK FOR FAILURES)
    open ( my $exOUT, ">>", 'Check_Target_Sequence.txt' ) or die "Could not open file 'Check_Target_Sequence.txt': $!";
    ###OPEN 'emm_region_extract.fasta' for READING (CHECK FOR FAILURES)
    open ( my $newEx, "<", 'emm_region_extract.fasta' ) or die "Could not open file 'emm_region_extract.fasta': $!";
    ###READ EACH LINE OF FILE B.txt (BAR) and add it to FILE A.txt (FOO)
    #print $exOUT "##############################--NEW EMM TYPE SEQUENCE--##############################\n";
    print $exOUT '#' x 65;
    print $exOUT "--NEW EMM TYPE SEQUENCE--";
    print $exOUT '#' x 65;
    print $exOUT "\n\n";    
    while ( my $line = <$newEx> ) {
	print $exOUT $line;
    }
    #print $exOUT "#####################################################################################\n";    
    print $exOUT '#' x 150;
    print $exOUT "\n\n";
    close $exOUT;
}
close $fh;

my $old_name = "emm-Type_Results.txt";
my $emm_out = $outName."__emm-Type__Results.txt";
rename $old_name, $emm_out;

###This will be part of --debug flag stuff###
#unlink($fastq1_trimd);
#unlink($fastq2_trimd);
#unlink("contig-vs-frwd_nucl.txt");
#unlink("emm_region_extract.bed");
#unlink("emm_vs_DB_nucl.txt");
#unlink("emm_region_extract.fasta");
