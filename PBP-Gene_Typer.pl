#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use File::Basename;
use File::Spec;

###MODULE LOAD###
#module load perl/5.12.3
#module load ncbi-blast+/2.2.29
#module load BEDTools/2.17.0
#module load Python/2.7

sub checkOptions {
    my %opts;
    getopts('h1:2:r:o:n:s:p:', \%opts);
    my ($help, $fastq1, $fastq2, $PBP_DB, $outDir, $outName, $pbp_genes, $species, @pbp_input);

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
        $PBP_DB = $opts{r};
        if (-e $PBP_DB) {
            print "The PBP reference database sequence: $PBP_DB\n";
        } else {
            print "The PBP reference sequence location is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/reference_DB).\n";
            help();
        }
    } else {
        print "The PBP reference sequence location (including full path) has not been given.\n";
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

    if($opts{s}) {
	if ($opts{s} eq "GAS" || $opts{s} eq "GBS" || $opts{s} eq "SPN") {
	    $species = $opts{s};
	    print "The species is: $species\n";
	} else {
	    print "The species argument needs to be one of the following: GAS, GBS or SPN\n";
	    help();
	}
    } else {
        print "The species argument hasn't been provided\n";
	help();
    }

    if($opts{p}) {
        $pbp_genes = $opts{p};
	@pbp_input = split(/,/,$pbp_genes);
	foreach my $pbp (@pbp_input) {
	    if ($pbp =~ /^1A$|^2B$|^2X$/) {
		print "The script will extract the following PBP gene: $pbp\n";
	    } else {
		print "The PBP gene <$pbp> doesn't exist or isn't in the right format.\n";
		print "Make sure you provide a comma delimited string containing 1A, 2B or 2X (e.g. '1A,2B,2X').\n";
		help();
	    }
        }
    } else {
        print "The PBP gene input argument has not been given.\n";
        help();
    }

    return ($help, $fastq1, $fastq2, $PBP_DB, $outDir, $outName, $species, @pbp_input);
}

sub help
{

die <<EOF

USAGE
GBS_PBP-Gene_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <reference databases directory: file path> -s <species name: string> -p <pbp genes: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   PBP reference sequence file path (including full path)
    -o   output directory
    -n   output name prefix
    -s   type of strep species. Input argument must be GAS, GBS or SPN.
    -p   comma delimited string that defines which PBP genes should be extracted (1A,2B,2X)

EOF
}

my ($help, $fastq1, $fastq2, $PBP_DB, $outDir, $outName, $species, @pbp_input) = checkOptions( @ARGV );




###SUBROUTINES###
sub PBP_blastTyper {
    my ($pbp_type,$pbp_name,$pbp_extract) = @_;
    chomp($pbp_type);
    chomp($pbp_name);
    chomp($pbp_extract);
    my $pbp_out;
    #print "type: $pbp_type || name: $pbp_name || LoTrac Extract: $pbp_extract\n";

    my $query_seq = extractFastaByID($pbp_name,$pbp_extract);
    #my $query_length = fasta_seq_length($query_seq);
    #open ( my $qOUT, ">", 'TEMP_query_sequence.fna' ) or die "Could not open file TEMP_query_sequence.fna: $!";
    #print $qOUT $query_seq;
    #close $qOUT;
    ##Translate TEMP_query_sequence.fna to .faa file
    #system("Translate_DNA-6Frame.pl -s TEMP_query_sequence.fna -f 1 -l 3 > TEMP_query_sequence.faa");
    #my $query_length = fasta_seq_length("TEMP_query_sequence.faa");
    my ($query_hedr, $query_aaSeq) = sixFrame_Translate($query_seq,1);
    `printf ">$query_hedr\n$query_aaSeq\n" > TEMP_query_sequence.faa`;
    my $query_length = fasta_seq_length($query_aaSeq);

    my $db_path = dirname($PBP_DB);
    my $blastDB_name = "Blast_bLactam_".$pbp_type."_prot_DB";
    my $blast_seq = $species."_bLactam_".$pbp_type."-DB.faa";
    my $blast_out = "TEMP_".$outName."_blast-out_".$pbp_type.".txt";
    print "Blast DB name: $db_path/$blastDB_name\n";
    if (!(glob("$db_path/$blastDB_name*"))) {
	print "Need to make a new Blast database\n";
	system("makeblastdb -in $db_path/$blast_seq -dbtype prot -out $db_path/$blastDB_name");
	system("blastp -db $db_path/$blastDB_name -query TEMP_query_sequence.faa -outfmt 6 -out $blast_out");
    } else {
	print "Blast database has already been created\n";
        system("blastp -db $db_path/$blastDB_name -query TEMP_query_sequence.faa -outfmt 6 -out $blast_out");
    }

    my $bestHit = `cat $blast_out | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
    print "best hit info: $bestHit";
    my @bestArray = split('\t',$bestHit);
    my $best_name = $bestArray[1];
    my $best_len = $bestArray[3];
    my $best_iden = $bestArray[2];
    my $frag_length = $best_len / $query_length;

    print "name of best hit in the PBP database: $bestArray[1]\n";
    print "identity of best hit in the PBP database: $bestArray[2]\n";
    print "length of best hit in the PBP database: $bestArray[3]\n";

    if ($best_iden == 100 && $best_len == $query_length) {
	print "Found a match\n\n";
	($pbp_out = $best_name) =~ s/^([0-9]+)\|\|.*/$1/g;
    } elsif (-s $blast_out && $best_iden >= 50 && $frag_length >= 0.50) {
	print "Didn't find match.  The sequence needs to be added to the database using 'bLactam-PBP_Updater.sh'\n\n";
	$pbp_out = "NEW";
	my $new_PBPseq = "PBP_".$pbp_type."_query_sequence.faa";
	rename("TEMP_query_sequence.faa",$new_PBPseq);
	my $fragPath = File::Spec->rel2abs("$new_PBPseq");
	open(my $f_new,'>>',"TEMP_newPBP_allele_info.txt") or die "Could not open file 'TEMP_newPBP_allele_info.txt' $!";
	print $f_new "$outName\t$fragPath\t$pbp_type\n";
	close $f_new;
    } else {
        $pbp_out = "ERROR";
    }
return $pbp_out;
}

sub sixFrame_Translate {
    my ($seq_input,$opt_f) = @_;

    sub codon2aa{
        my($codon)=@_;
        $codon=uc $codon;
        my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

        if(exists $g{$codon}){return $g{$codon};}
        elsif($codon=~/GC./i){return 'A';}
        elsif($codon=~/GG./i){return 'G';}
        elsif($codon=~/CC./i){return 'P';}
        elsif($codon=~/AC./i){return 'T';}
        elsif($codon=~/GT./i){return 'V';}
        elsif($codon=~/CG./i){return 'R';}
        elsif($codon=~/TC./i){return 'S';}
        else {
            return('x');
            print "Bad codon \"$codon\"!!\n";
        }
    }

    (my $DNAheader,my @DANseq)=split(/\n/,$seq_input);
    chomp $DNAheader;
    $DNAheader=~s/\s+$//g;
    my $DNAseq = join( '',@DANseq);
    $DNAseq =~ s/\s//g;
    $DNAheader=~s/>//g;
    $DNAseq=~s/>//g;
    my$DNA_length=length$DNAseq;
    #print "\nSeq:$DNAheader\t:$DNA_length nt\n\n";
    my $DNArevSeq = reverse($DNAseq);
    $DNArevSeq=~tr/ATGCatgc/TACGtacg/;
    #print "\nThe original DNA sequence is:\n$DNAseq \nThe reverse of DNA sequence is:\n$DNArevSeq\n";
    my @protein='';
    my @dna='';
    my $codon1;

    if ($opt_f == 1) {
        for(my $i=0;$i<(length($DNAseq)-2);$i+=3){
            $codon1=substr($DNAseq,$i,3);
            $protein[1].= codon2aa($codon1);
            #$dna[1].=codon2nt($codon1);
            }
        }
    if ($opt_f == 2) {
            my $codon2;
            for(my $i=1;$i<(length($DNAseq)-2);$i+=3){
                $codon2=substr($DNAseq,$i,3);
                $protein[2].= codon2aa($codon2);
                #$dna[2].=codon2nt($codon2);
                }
    }
    if ($opt_f == 3) {
            my $codon3;
            for(my $i=2;$i<(length($DNAseq)-2);$i+=3){
                $codon3=substr($DNAseq,$i,3);
                $protein[3].= codon2aa($codon3);
                #$dna[3].=codon2nt($codon3);
                }
    }
    if ($opt_f == 4) {
            my $codon4;
            for(my $i=0;$i<(length($DNArevSeq)-2);$i+=3){
                $codon4=substr($DNArevSeq,$i,3);
                $protein[4].= codon2aa($codon4);
                #$dna[4].=codon2nt($codon4);
                }
    }
    if ($opt_f == 5) {
            my $codon5;
            for(my $i=1;$i<(length($DNArevSeq)-2);$i+=3){
                $codon5=substr($DNArevSeq,$i,3);
                $protein[5].= codon2aa($codon5);
                #$dna[5].=codon2nt($codon5);
                }
    }
    if ($opt_f == 6) {
            my $codon6;
                for(my $i=2;$i<(length($DNArevSeq)-2);$i+=3){
                    $codon6=substr($DNArevSeq,$i,3);
                    $protein[6].= codon2aa($codon6);
                    #$dna[6].=codon2nt($codon6);
                    }
    }
#print "translate result\n$protein[$opt_f]\n";
return ($DNAheader,$protein[$opt_f]);
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

sub fasta_seq_length {
    my ($seq) = @_;
    #open ( my $q_seq, "<", $seq ) or die "Could not open file '$seq': $!";
    #my @lines = split /\n/, $q_seq;
    my @lines = split /\n/, $seq;
    my $final_line;
    foreach my $line (@lines) {
    #while (my $line = <$q_seq>) {
        chomp($line);
        if ($line =~ /^>/) {
            next;
        } else {
            #print "line: $line\n";
            $final_line = $final_line.$line;
        }
    }
    chomp($final_line);
    #print "final line: $final_line\n";
    return length($final_line);
}



#=pod
###Start Doing Stuff###
print "\n";
my $justName = `echo $outName | sed 's/PBP_//g'`;
chomp($justName);
system("LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q $PBP_DB -S 2.2M -f -n $justName -o $outDir");
chdir "$outDir";
my $PBP_output = "TEMP_pbpID_Results.txt";
open(my $fh,'>',$PBP_output) or die "Could not open file '$PBP_output' $!";
if (glob("ERROR_*1A*.fasta") || glob("ERROR_*2B*.fasta") || glob("ERROR_*2X*.fasta")) {
    print $fh "No PBP Type\n";
    exit
}

my ($code_1A, $code_2B, $code_2X);
foreach my $pbp (@pbp_input) {
    if ($pbp eq "1A") {
	my $pbp_1A = glob("EXTRACT_*1A*.fasta");
	my $pbp1A_fragName = `cat $pbp_1A | grep ">" | tail -n1 | sed 's/>//g'`;
	chomp($pbp1A_fragName);
	$code_1A = PBP_blastTyper("1A",$pbp1A_fragName,$pbp_1A);
        print "pbp extract file: $pbp_1A || pbp extraction name: $pbp1A_fragName || pbp ID is: $code_1A\n";
    } elsif ($pbp eq "2B") {
	my $pbp_2B = glob("EXTRACT_*2B*.fasta");
	my $pbp2B_fragName = `cat $pbp_2B | grep ">" | tail -n1 | sed 's/>//g'`;
	chomp($pbp2B_fragName);
	$code_2B = PBP_blastTyper("2B",$pbp2B_fragName,$pbp_2B);
    } elsif ($pbp eq "2X") {
	my $pbp_2X = glob("EXTRACT_*2X*.fasta");
	my $pbp2X_fragName = `cat $pbp_2X | grep ">" | tail -n1 | sed 's/>//g'`;
	chomp($pbp2X_fragName);
	$code_2X = PBP_blastTyper("2X",$pbp2X_fragName,$pbp_2X);
    }
}

if ($species eq "GBS") {
    print $fh "Sample_Name\tPBP_Code(1A:2X)\n";
    print $fh "$outName\t$code_1A:$code_2X\n";
} elsif ($species eq "SPN") {
    print $fh "Sample_Name\tPBP_Code(1A:2B:2X)\n";
    print $fh "$outName\t$code_1A:$code_2B:$code_2X\n";
} elsif ($species eq "GAS") {
    print $fh "Sample_Name\tPBP_Code(1A:2B:2X)\n";
    print $fh "$outName\t$code_2X\n";
}
close $fh;
#=cut







#my $pbp_1A = glob("EXTRACT_*1A*.fasta");
#my $pbp1A_fragName = `cat $pbp_1A | grep ">" | tail -n1 | sed 's/>//g'`;
#my $pbp_2B = glob("EXTRACT_*2B*.fasta");
#my $pbp2B_fragName = `cat $pbp_2B | grep ">" | tail -n1 | sed 's/>//g'`;
#my $pbp_2X = glob("EXTRACT_*2X*.fasta");
#my $pbp2X_fragName = `cat $pbp_2X | grep ">" | tail -n1 | sed 's/>//g'`;
#my $code_1A = PBP_blastTyper("1A",$pbp1A_fragName,$pbp_1A);
#my $code_2B = PBP_blastTyper("2B",$pbp2B_fragName,$pbp_2B);
#my $code_2X = PBP_blastTyper("2X",$pbp2X_fragName,$pbp_2X);
#print $fh "$outName\t$code_1A\t$code_2X\t$code_1A:$code_2X\n";
#print $fh "$outName\t$code_1A:$code_2X\n";
#close $fh;

###Delete Temp Files###
#unlink(TEMP*);
