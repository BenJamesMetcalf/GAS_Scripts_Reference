#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;
use File::Copy qw(copy);
use Env;
use lib $ENV{MODULESHOME}."/init";
use perl;

###MODULE LOAD###
#module load samtools/0.1.18
#module load bowtie2/2.1.0
#module load Python/2.7
#module load freebayes/0.9.21

sub checkOptions {
    my %opts;
    getopts('h1:2:d:r:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $ref_dir, $res_DB, $outDir, $outName);

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

    if($opts{d}) {
        $ref_dir = $opts{d};
        $ref_dir =~ s/\/$//g;
        if (-d $ref_dir) {
            print "Directory containing the  resistance reference sequences: $ref_dir\n";
        } else {
            print "The directory containing the resistance references is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/ref_dir).\n";
            help();
        }
    } else {
        print "The resistance reference sequence directory (including full path) has not been given.\n";
        help();
    }

    if($opts{r}) {
        $res_DB = "$ref_dir/$opts{r}";
        if ($res_DB) {
            print "The resistance reference sequence file: $opts{r}\n";
        } else {
            print "The resistance reference sequence file is not in the correct format or doesn't exist.\n";
            #print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The resistance reference sequence file has not been given.\n";
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

    return ($help, $fastq1, $fastq2, $ref_dir, $res_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
GAS_Res_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -d <ref directory: dir> -r <resistance seq: file> -o <output directory name: string> -n <output name prefix: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -d   reference sequence directory (including full path)
    -r   resistance reference sequence file
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $ref_dir, $res_DB, $outDir, $outName) = checkOptions( @ARGV );




###Subroutines###
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
return $protein[$opt_f];
}

sub extractSeqByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
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
            #$output = ">$id\n$seq\n";
            $output = $seq;
            last;
        }
    }
    return $output;
}

sub extractFastaByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
    local $/ = "\n>";  # read by FASTA record

    my $output;
    while (my $seq = <$fh>) {
        chomp $seq;
        #print "while seq:\n$seq\n";
        my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        if ($id eq $lookup) {
            $seq =~ s/^>*.+\n//;  # remove FASTA header
            $seq =~ s/\n//g;  # remove endlines
            #print ">$id\n";
            #print "$seq\n";
            $output = ">$id\n$seq\n";
            #$output = $seq;
            last;
        }
    }
    return $output;
}

sub freebayes_prior_fix {
    my ($bamFile, $refFile, $target) = @_;
    (my $samFile = $bamFile) =~ s/\.bam/\.sam/g;
    system("samtools view -h $bamFile > $samFile");
    system("cat $samFile | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > CHECK_target_seq.sam");
    system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $samFile >> CHECK_target_seq.sam");
    system("samtools view -bS CHECK_target_seq.sam > CHECK_target_seq.bam");
    system("samtools index CHECK_target_seq.bam CHECK_target_seq.bai");
    my $REF_seq = extractFastaByID("$target","$refFile");
    open(my $rf,'>',"CHECK_target_ref.fna");
    print $rf "$REF_seq\n";
    close $rf;
    system("freebayes -q 20 -p 1 -f CHECK_target_ref.fna CHECK_target_seq.bam -v CHECK_target_seq.vcf");
    system("bgzip CHECK_target_seq.vcf");
    system("tabix -p vcf CHECK_target_seq.vcf.gz");
    my $extractSeq = `echo "$REF_seq" | vcf-consensus CHECK_target_seq.vcf.gz`;
    chomp($extractSeq);
    #print "$target-----------------------------------\n";
    #system("cat CHECK_target_seq.sam");
    #system("zcat CHECK_target_seq.vcf.gz");
    #print "reference seq:\n$REF_seq\n";
    #print "extracted Seq:\n$extractSeq\n";
    #print "$target-----------------------------------\n";
    system("rm CHECK_target*");
    return $extractSeq;
}




###Start Doing Stuff###
chdir "$outDir";
my $Res_output = "TEMP_Res_Results.txt";
open(my $fh,'>',$Res_output) or die "Could not open file '$Res_output' $!";
my $BIN_res_out = "BIN_Res_Results.txt";
open(my $bh,'>',$BIN_res_out) or die "Could not open file '$BIN_res_out' $!";
my @Bin_Res_arr = (0) x 20;
#print $fh "Resistance_Group\tTarget\n";
#=pod
my $outNameRES = "RES_".$outName;
my $out_nameARG = "ARG_".$outName;
my $out_nameRESFI = "RESFI_".$outName;
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outNameRES --log --save_scores --min_coverage 99.9 --max_divergence 5 --gene_db $res_DB");
###Type ARG-ANNOT Resistance Genes###
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $out_nameARG --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db $ref_dir/ARGannot_r1.fasta");
###Type ResFinder Resistance Gene###
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $out_nameRESFI --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db $ref_dir/ResFinder.fasta");
#=cut

my @TEMP_RES_bam = glob("RES_*\.sorted\.bam");
my @TEMP_RES_fullgene = glob("RES_*__fullgenes__*__results\.txt");
my $RES_bam = $TEMP_RES_bam[0];
my $RES_full_name = $TEMP_RES_fullgene[0];
print "res bam is: $RES_bam || res full gene $RES_full_name\n";
(my $RES_vcf = $RES_bam) =~ s/\.bam/\.vcf/g;
(my $RES_bai = $RES_bam) =~ s/\.bam/\.bai/g;

my @TEMP_ARG_fullgene = glob("ARG_*__fullgenes__*__results\.txt");
my $ARG_full_name = $TEMP_ARG_fullgene[0];
my @TEMP_RESFI_fullgene = glob("RESFI_*__fullgenes__*__results\.txt");
my $RESFI_full_name = $TEMP_RESFI_fullgene[0];
my $merged_net = "ARG-RESFI_fullgenes_results.txt";
#copy $ARG_full_name, $merged_net;
system("tail -n+2 $ARG_full_name > $merged_net");
system("tail -n+2 $RESFI_full_name >> $merged_net");

my %drugRes_Col = (
    "ER_CL" => "neg",
    "TET" => "neg",
    #"RARE" => "neg",
    "GYRASPY_PARC" => "neg",
    #"FOLA_FOLP" => "neg",
    "OTHER" => "neg",
    );

my %Res_Targets = (
    "ERM" => "neg",
    "LNUB" => "neg",
    "LSA" => "neg",
    "MEF" => "neg",
    "TET" => "neg",
    "CAT" => "neg",
    "GYRA" => "neg",
    "PARC" => "neg",
    "23S1" => "neg",
    "23S3" => "neg",
    "FOLA" => "neg",
    "FOLP1" => "neg",
    "FOLP2" => "neg",
    "RPOB1" => "neg",
    "RPOBN" => "neg",
    );

###Type the Presence/Absence Targets###
open(MYINPUTFILE, "$RES_full_name");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @miscR_fullgene;
    @miscR_fullgene = split('\t',$line);
    if ($miscR_fullgene[5] >= 10) {
	if ($miscR_fullgene[3] =~ m/(ERM|LNUB|LSA|MEF)/) {
	    if ($drugRes_Col{"ER_CL"} eq "neg") {
		$drugRes_Col{"ER_CL"} = $miscR_fullgene[2];
	    } else {
		my $new_val = $drugRes_Col{"ER_CL"}.":".$miscR_fullgene[2];
		$drugRes_Col{"ER_CL"} = $new_val;
	    }
	}
	if ($miscR_fullgene[3] =~ m/TET/) {
	    if ($drugRes_Col{"TET"} eq "neg") {
                $drugRes_Col{"TET"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"TET"}.":".$miscR_fullgene[2];
                $drugRes_Col{"TET"} = $new_val;
            }
	}
	if ($miscR_fullgene[3] =~ m/CAT/) {
            if ($drugRes_Col{"OTHER"} eq "neg") {
                $drugRes_Col{"OTHER"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"OTHER"}.":".$miscR_fullgene[2];
                $drugRes_Col{"OTHER"} = $new_val;
            }
        }

	if ($miscR_fullgene[3] =~ m/ERM/) {
	    $Res_Targets{"ERM"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/LNUB/) {
	    $Res_Targets{"LNUB"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/LSA/) {
	    $Res_Targets{"LSA"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/MEF/) {
            $Res_Targets{"MEF"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/TET/) {
            $Res_Targets{"TET"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/CAT/) {
            $Res_Targets{"CAT"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/PARC/) {
	    $Res_Targets{"PARC"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/GYRA/) {
	    $Res_Targets{"GYRA"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/23S1/) {
            $Res_Targets{"23S1"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/23S3/) {
            $Res_Targets{"23S3"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/FOLA/) {
            $Res_Targets{"FOLA"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/FOLP1/) {
            $Res_Targets{"FOLP1"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/FOLP2/) {
            $Res_Targets{"FOLP2"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/RPOB1/) {
            $Res_Targets{"RPOB1"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/RPOBN/) {
            $Res_Targets{"RPOBN"} = "pos";
	}
    }
}

#while (my ($key, $val) = each %Res_Targets) {
#    my @val_arr = split(':',$val);
#    my @val_sort = sort(@val_arr);
#    my $val_out = join(':',@val_sort);
#    print "$key\t$val_out\n";
#}
#print "\n";

###############################################################################################
###ARG-ANNOT and ResFinder Safety Net###
open(MYINPUTFILE, "$merged_net");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @miscR_fullgene;
    @miscR_fullgene = split('\t',$line);
    if ($miscR_fullgene[5] >= 10) {
        if ($miscR_fullgene[3] =~ m/ERM/i) {#&& $Res_Targets{"ERM"} eq "neg") {
	    if ($Res_Targets{"ERM"} eq "neg") {
		if ($drugRes_Col{"ER_CL"} eq "neg") {
		    $drugRes_Col{"ER_CL"} = "ERM";
		} else {
		    my $new_val = $drugRes_Col{"ER_CL"}.":ERM";
		    $drugRes_Col{"ER_CL"} = $new_val;
		}
		$Res_Targets{"ERM"} = "pos";
	    }
	} elsif ($miscR_fullgene[3] =~ m/LNU/i) {#&& $Res_Targets{"LNUB"} eq "neg") {
	    if ($Res_Targets{"LNUB"} eq "neg") {
		if ($drugRes_Col{"ER_CL"} eq "neg") {
		    $drugRes_Col{"ER_CL"} = "LNU";
		} else {
		    my $new_val = $drugRes_Col{"ER_CL"}.":LNU";
		    $drugRes_Col{"ER_CL"} = $new_val;
		}
		$Res_Targets{"LNUB"} = "pos";
	    }
        } elsif ($miscR_fullgene[3] =~ m/LSA/i) {#&& $Res_Targets{"LSA"} eq "neg") {
	    if ($Res_Targets{"LSA"} eq "neg") {
		if ($drugRes_Col{"ER_CL"} eq "neg") {
		    $drugRes_Col{"ER_CL"} = "LSA";
		} else {
		    my $new_val = $drugRes_Col{"ER_CL"}.":LSA";
		    $drugRes_Col{"ER_CL"} = $new_val;
		}
		$Res_Targets{"LSA"} = "pos";
	    }
        } elsif ($miscR_fullgene[3] =~ m/MEF/i) {#&& $Res_Targets{"MEF"} eq "neg") {
	    if ($Res_Targets{"MEF"} eq "neg") {
		if ($drugRes_Col{"ER_CL"} eq "neg") {
		    $drugRes_Col{"ER_CL"} = "MEF";
		} else {
		    my $new_val = $drugRes_Col{"ER_CL"}.":MEF";
		    $drugRes_Col{"ER_CL"} = $new_val;
		}
		$Res_Targets{"MEF"} = "pos";
	    }
        } elsif ($miscR_fullgene[3] =~ m/TET/i) { #&& $Res_Targets{"TET"} eq "neg") {
	    if ($Res_Targets{"TET"} eq "neg") {
		if ($drugRes_Col{"TET"} eq "neg") {
		    $drugRes_Col{"TET"} = "TET";
		} else {
		    my $new_val = $drugRes_Col{"TET"}.":TET";
		    $drugRes_Col{"TET"} = $new_val;
		}
		$Res_Targets{"TET"} = "pos";
	    }
        } elsif ($miscR_fullgene[3] =~ m/CAT/i) {#&& $Res_Targets{"CAT"} eq "neg") {
	    if ($Res_Targets{"CAT"} eq "neg") {
		if ($drugRes_Col{"OTHER"} eq "neg") {
		    $drugRes_Col{"OTHER"} = "CAT";
		} else {
		    my $new_val = $drugRes_Col{"OTHER"}.":CAT";
		    $drugRes_Col{"OTHER"} = $new_val;
		}
		$Res_Targets{"CAT"} = "pos";
	    }
        } else {
	    if ($drugRes_Col{"OTHER"} eq "neg") {
                $drugRes_Col{"OTHER"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"OTHER"}.":".$miscR_fullgene[2];
                $drugRes_Col{"OTHER"} = $new_val;
            }
	}
    }
}
###############################################################################################


###Now Type the AA/nt Change Targets###
###############################################################################################
###Type the PARC and GYRA FLQ Targets###
if ($Res_Targets{"PARC"} eq "pos") {
    my @PARC_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $PARC_seq = freebayes_prior_fix($RES_bam, $res_DB,"15__PARC__PARC-1__15");
    my $PARC_aaSeq = sixFrame_Translate($PARC_seq,1);
    my $PARC_aaRef = "HPHGDSSIYDAMVRMSQ";
    print "PARCgbs sequence: $PARC_seq || $PARC_aaSeq\n";
    if ($PARC_aaSeq ne "HPHGDSSIYDAMVRMSQ") {
        my $mask = $PARC_aaSeq ^ $PARC_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($PARC_aaRef,$-[0],1), ' ', substr($PARC_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
	    my $diff_element = substr($PARC_aaRef,$-[0],1).($-[0]+1).substr($PARC_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "PARC seq: $PARC_seq\n";
        my $diff_output = join(',',@seq_diffs);
        my $bin_out = join(':',@seq_diffs);
        $Bin_Res_arr[14] = $bin_out;
	my $PARC_out = "PARC-".$diff_output;
	if ($drugRes_Col{"GYRASPY_PARC"} eq "neg") {
	    $drugRes_Col{"GYRASPY_PARC"} = $PARC_out;
	} else {
	    my $new_val = $drugRes_Col{"GYRASPY_PARC"}.":".$PARC_out;
	    $drugRes_Col{"GYRASPY_PARC"} = $new_val;
	}
    }
}

if ($Res_Targets{"GYRA"} eq "pos") {
    my @GYRA_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $GYRA_seq = freebayes_prior_fix($RES_bam, $res_DB,"10__GYRA__GYRA-1__10");
    my $GYRA_aaSeq = sixFrame_Translate($GYRA_seq,1);
    my $GYRA_aaRef = "VMGKYHPHGDSSIYEAMVRMAQWWSY";
    print "GYRA sequence: $GYRA_seq || $GYRA_aaSeq\n";
    if ($GYRA_aaSeq ne "VMGKYHPHGDSSIYEAMVRMAQWWSY") {
        my $mask = $GYRA_aaSeq ^ $GYRA_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($GYRA_aaRef,$-[0],1), ' ', substr($GYRA_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($GYRA_aaRef,$-[0],1).($-[0]+1).substr($GYRA_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "GYRA seq: $GYRA_seq\n";
        my $diff_output = join(',',@seq_diffs);
	my $bin_out = join(':',@seq_diffs);
	$Bin_Res_arr[9] = $bin_out;
        my $GYRA_out = "GYRA-".$diff_output;
        if ($drugRes_Col{"GYRASPY_PARC"} eq "neg") {
            $drugRes_Col{"GYRASPY_PARC"} = $GYRA_out;
        } else {
            my $new_val = $drugRes_Col{"GYRASPY_PARC"}.":".$GYRA_out;
            $drugRes_Col{"GYRASPY_PARC"} = $new_val;
        }
    }
}
###############################################################################################

###############################################################################################
###Type the 23S ribosomal RNA MLS resistance target###
if ($Res_Targets{"23S1"} eq "pos") {
    my @S123_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $PRE_S123_seq = freebayes_prior_fix($RES_bam, $res_DB,"1__23S1__23S1-1__1");
    my @ARR_S123_seq = split('\n',$PRE_S123_seq);
    my $S123_seq = $ARR_S123_seq[1];
    #my $23S1_aaSeq = sixFrame_Translate($GYRA_seq,1);
    my $S123_ntRef = "GTTACCCGCGACAGGACGGAAAGACCCCATGGAG";
    print "S123 sequence: $S123_seq\n";
    if ($S123_seq ne "GTTACCCGCGACAGGACGGAAAGACCCCATGGAG") {
        my $mask = $S123_seq ^ $S123_ntRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($S123_ntRef,$-[0],1), ' ', substr($S123_seq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($S123_ntRef,$-[0],1).($-[0]+1).substr($S123_seq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "23S1 seq: $S123_seq\n";
        my $diff_output = join(',',@seq_diffs);
	my $bin_out = join(':',@seq_diffs);
	$Bin_Res_arr[0] = $bin_out;
        my $S123_out = "23S1-".$diff_output;
        if ($drugRes_Col{"ER_CL"} eq "neg") {
            $drugRes_Col{"ER_CL"} = $S123_out;
        } else {
            my $new_val = $drugRes_Col{"ER_CL"}.":".$S123_out;
            $drugRes_Col{"ER_CL"} = $new_val;
        }
    }
}

if ($Res_Targets{"23S3"} eq "pos") {
    my @S323_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $PRE_S323_seq = freebayes_prior_fix($RES_bam, $res_DB,"2__23S3__23S3-1__2");
    my @ARR_S323_seq = split('\n',$PRE_S323_seq);
    my $S323_seq = $ARR_S323_seq[1];
    #my $23S1_aaSeq = sixFrame_Translate($GYRA_seq,1);
    my $S323_ntRef = "CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG";
    print "S323 sequence: $S323_seq\n";
    if ($S323_seq ne "CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG") {
        my $mask = $S323_seq ^ $S323_ntRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($S323_ntRef,$-[0],1), ' ', substr($S323_seq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($S323_ntRef,$-[0],1).($-[0]+1).substr($S323_seq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "23S3 seq: $S323_seq\n";
        my $diff_output = join(',',@seq_diffs);
	my $bin_out = join(':',@seq_diffs);
	$Bin_Res_arr[1] = $bin_out;
        my $S323_out = "23S3-".$diff_output;
        if ($drugRes_Col{"ER_CL"} eq "neg") {
            $drugRes_Col{"ER_CL"} = $S323_out;
        } else {
            my $new_val = $drugRes_Col{"ER_CL"}.":".$S323_out;
            $drugRes_Col{"ER_CL"} = $new_val;
        }
    }
}
###############################################################################################

###############################################################################################
###Type the FOLA-1 Trimethoprim targets###
###The FOLA-1 marker for tmR resistance may contain mosiac regions and must be extracted from the genome assembly###
my $REF_seq = extractFastaByID("7__FOLA__FOLA-1__7","$res_DB");
`echo "$REF_seq" > TEMP_FOLA_Ref.fna`;
module "unload perl/5.22.1";
module "load perl/5.16.1-MT";
system("LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q TEMP_FOLA_Ref.fna -S 2.2M -f -n $outName -o $outDir");
module "unload perl/5.16.1-MT";
module "load perl/5.22.1";
my $FOLA_file = glob("EXTRACT_*FOLA*.fasta");
my $FOLA_error = glob("ERROR_*FOLA*.fasta");
my @FOLA_output;
if ($FOLA_file && ! $FOLA_error) {
    my $FOLA_seq_hedr = `cat $FOLA_file | grep ">" | tail -n1 | sed 's/>//g'`;
    chomp($FOLA_seq_hedr);
    my $FOLA_seq = extractFastaByID("$FOLA_seq_hedr","$FOLA_file");
    chomp($FOLA_seq);
    #print "FOLA extract header: $FOLA_seq_hedr | FOLA extract file: $FOLA_file\n";
    my $FOLA_aaSeq = sixFrame_Translate($FOLA_seq,1);
    my $FOLA_aaRef = "TSIEKVLEWYHAQDKTLYII";
    print "FOLA Seq: $FOLA_aaRef || $FOLA_aaSeq || $FOLA_seq\n";
    if ($FOLA_aaSeq ne "TSIEKVLEWYHAQDKTLYII") {
	my $mask = $FOLA_aaSeq ^ $FOLA_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($FOLA_aaRef,$-[0],1), ' ', substr($FOLA_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($FOLA_aaRef,$-[0],1).($-[0]+1).substr($FOLA_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "FOLA seq: $FOLA_seq\n";
        my $diff_output = join(',',@seq_diffs);
	my $bin_out = join(':',@seq_diffs);
	$Bin_Res_arr[6] = $bin_out;
        my $FOLA_out = "FOLA-".$diff_output;
        if ($drugRes_Col{"OTHER"} eq "neg") {
            $drugRes_Col{"OTHER"} = $FOLA_out;
        } else {
            my $new_val = $drugRes_Col{"OTHER"}.":".$FOLA_out;
            $drugRes_Col{"OTHER"} = $new_val;
        }
    }
}

###Type the FOLP-1 Sulfonamide resistance target###
my $target = "8__FOLP1__FOLP1-1__8";
(my $RES_sam = $RES_bam) =~ s/\.bam/\.sam/g;
system("samtools view -h $RES_bam > $RES_sam");
system("cat $RES_sam | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > FOLP_target_seq.sam");
system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $RES_sam >> FOLP_target_seq.sam");
system("samtools view -bS FOLP_target_seq.sam > FOLP_target_seq.bam");
system("samtools index FOLP_target_seq.bam FOLP_target_seq.bai");
$REF_seq = extractFastaByID("$target","$res_DB");
open(my $rf,'>',"FOLP_target_ref.fna");
print $rf "$REF_seq\n";
close $rf;
system("freebayes -q 20 -p 1 -f FOLP_target_ref.fna FOLP_target_seq.bam -v FOLP_target_seq.vcf");
open(MYINPUTFILE, "FOLP_target_seq.vcf");
my %srst2_seroT;
while(<MYINPUTFILE>) {
    my $line = $_;
    chomp($line);
    if ($line =~ /^8__FOLP1__FOLP1-1__8/) {
        my @FOLP_line = split('\t', $line);
        my $ref_allele = $FOLP_line[3];
        #my $ref_len = length($ref_allele);
        my $alt_allele = $FOLP_line[4];
        #my $alt_len = length($alt_allele);
        #my $lociDP = ($FOLP_line[7] =~ /DP=([0-9]+);/);
        #$FOLP_line[7] =~ /DPB=([0-9]+);/;
        $FOLP_line[7] =~ /DPB=(\d+\.?\d*);/;
        print "FOLP DP: $1 | ref allele: $ref_allele | alt allele: $alt_allele\n";
        if (length($ref_allele) != length($alt_allele)) {
	    $Bin_Res_arr[7] = 1;
	    if ($drugRes_Col{"OTHER"} eq "neg") {
		$drugRes_Col{"OTHER"} = "FOLP1";
	    } else {
		my $new_val = $drugRes_Col{"OTHER"}.":FOLP1";
		$drugRes_Col{"OTHER"} = $new_val;
	    }
        }
    }
}

if ($Res_Targets{"FOLP2"} eq "pos") {
    my @FOLP2_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $FOLP2_seq = freebayes_prior_fix($RES_bam, $res_DB,"9__FOLP2__FOLP2-1__9");
    my $FOLP2_aaSeq = sixFrame_Translate($FOLP2_seq,1);
    print "FOLP2 sequence: $FOLP2_seq || $FOLP2_aaSeq\n";
    my $FOLP2_p13 = substr($FOLP2_aaSeq,12,1);
    if ($FOLP2_p13 eq "G") {
	$Bin_Res_arr[8] = 1;
	if ($drugRes_Col{"OTHER"} eq "neg") {
	    $drugRes_Col{"OTHER"} = "FOLP2";
	} else {
	    my $new_val = $drugRes_Col{"OTHER"}.":FOLP2";
	    $drugRes_Col{"OTHER"} = $new_val;
	}
    }
}
###############################################################################################

###############################################################################################
###Type Rifampicin Resistance Using 4 RPOB Targets###
if ($Res_Targets{"RPOB1"} eq "pos") {
    my @RPOB1_output;
    my $RPOB1_seq = freebayes_prior_fix($RES_bam, $res_DB,"16__RPOB1__RPOB1-1__16");
    my $RPOB1_aaSeq = sixFrame_Translate($RPOB1_seq,1);
    my $RPOB1_aaRef = "FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL";
    print "RPOB1 sequence: $RPOB1_seq || $RPOB1_aaSeq\n";
    if ($RPOB1_aaSeq ne "FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL") {
        my $mask = $RPOB1_aaSeq ^ $RPOB1_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($RPOB1_aaRef,$-[0],1), ' ', substr($RPOB1_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($RPOB1_aaRef,$-[0],1).($-[0]+1).substr($RPOB1_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "RPOB1 seq: $RPOB1_seq\n";
        my $diff_output = join(',',@seq_diffs);
	my $bin_out = join(':',@seq_diffs);
	$Bin_Res_arr[15] = $bin_out;
        my $RPOB1_out = "RPOB1-".$diff_output;
        if ($drugRes_Col{"OTHER"} eq "neg") {
            $drugRes_Col{"OTHER"} = $RPOB1_out;
        } else {
            my $new_val = $drugRes_Col{"OTHER"}.":".$RPOB1_out;
            $drugRes_Col{"OTHER"} = $new_val;
        }
    }
}

if ($Res_Targets{"RPOBN"} eq "pos") {
    my @RPOBN_output;
    my $RPOBN_seq = freebayes_prior_fix($RES_bam, $res_DB,"17__RPOBN__RPOBN-1__17");
    my $RPOBN_aaSeq = sixFrame_Translate($RPOBN_seq,1);
    my $RPOBN_aaRef = "RIIVSQLVRSPGV";
    print "RPOBN sequence: $RPOBN_seq || $RPOBN_aaSeq\n";
    if ($RPOBN_aaSeq ne "RIIVSQLVRSPGV") {
        my $mask = $RPOBN_aaSeq ^ $RPOBN_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($RPOBN_aaRef,$-[0],1), ' ', substr($RPOBN_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($RPOBN_aaRef,$-[0],1).($-[0]+1).substr($RPOBN_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "RPOBN seq: $RPOBN_seq\n";
        my $diff_output = join(',',@seq_diffs);
	my $bin_out = join(':',@seq_diffs);
	$Bin_Res_arr[16] = $bin_out;
        my $RPOBN_out = "RPOBN-".$diff_output;
        if ($drugRes_Col{"OTHER"} eq "neg") {
            $drugRes_Col{"OTHER"} = $RPOBN_out;
        } else {
            my $new_val = $drugRes_Col{"OTHER"}.":".$RPOBN_out;
            $drugRes_Col{"OTHER"} = $new_val;
        }
    }
}
###############################################################################################

###############################################################################################
###Make Binary Output Table###
#my $BIN_res_out = "BIN_Res_Results.txt";
#open(my $bh,'>',$BIN_res_out) or die "Could not open file '$BIN_res_out' $!";
#my @Bin_Res_arr = (0) x 19;
open(MYINPUTFILE, "$RES_full_name");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @feat_fullgene;
    @feat_fullgene = split('\t',$line);
    if ($feat_fullgene[5] >= 10) {
        if ($feat_fullgene[3] =~ m/CAT/) {
            $Bin_Res_arr[2] = 1;
        }
        if ($feat_fullgene[3] =~ m/ERMB/) {
            $Bin_Res_arr[3] = 1;
        }
        if ($feat_fullgene[3] =~ m/ERMT/) {
            $Bin_Res_arr[4] = 1;
        }
        if ($feat_fullgene[3] =~ m/ERMA/) {
            $Bin_Res_arr[5] = 1;
        }
        if ($feat_fullgene[3] =~ m/LNUB/) {
            $Bin_Res_arr[10] = 1;
        }
        if ($feat_fullgene[3] =~ m/LSAC/) {
            $Bin_Res_arr[11] = 1;
        }
        if ($feat_fullgene[3] =~ m/LSAE/) {
            $Bin_Res_arr[12] = 1;
        }
        if ($feat_fullgene[3] =~ m/MEF/) {
            $Bin_Res_arr[13] = 1;
        }
        if ($feat_fullgene[3] =~ m/TETL/) {
            $Bin_Res_arr[17] = 1;
        }
        if ($feat_fullgene[3] =~ m/TETM/) {
            $Bin_Res_arr[18] = 1;
        }
        if ($feat_fullgene[3] =~ m/TETO/) {
            $Bin_Res_arr[19] = 1;
        }
    }
}

###Print GAS Binary Output###
print $bh join(',',@Bin_Res_arr);
###############################################################################################

###############################################################################################
###Print Drug Resistance Output###
while (my ($key, $val) = each %drugRes_Col) {
    my @val_arr = split(':',$val);
    #print "@val_arr\n";
    my @val_sort = sort { "\L$a" cmp "\L$b" } @val_arr;
    #print "@val_sort\n";
    my $val_out = join(':',@val_sort);
    print "$key\t$val_out\n";
    $drugRes_Col{"$key"} = $val_out;
    #print $fh "$key\t$val_out\n";
}

print $fh "ER_CL\t$drugRes_Col{'ER_CL'}\n";
print $fh "TET\t$drugRes_Col{'TET'}\n";
print $fh "GYRA_PARC\t$drugRes_Col{'GYRASPY_PARC'}\n";
#print $fh "FOLA_FOLP\t$drugRes_Col{'FOLA_FOLP'}\n";
#print $fh "RARE\t$drugRes_Col{'RARE'}\n";
print $fh "OTHER\t$drugRes_Col{'OTHER'}\n";
###############################################################################################
