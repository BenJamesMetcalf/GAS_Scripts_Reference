#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;

###MODULE LOAD###
#module load samtools/0.1.18
#module load bowtie2/2.1.0
#module load Python/2.7
#module load freebayes/0.9.21

sub checkOptions {
    my %opts;
    getopts('h1:2:d:f:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $ref_dir, $feat_DB, $outDir, $outName);

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
            print "Directory containing the surface and secretory protein reference sequences: $ref_dir\n";
        } else {
            print "The directory containing the surface and secretory protein references is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The surface protein reference sequence directory (including full path) has not been given.\n";
        help();
    }

    if($opts{f}) {
        $feat_DB = "$ref_dir/$opts{f}";
        if ($feat_DB) {
            print "The surface and secretory protein reference sequence file: $opts{f}\n";
        } else {
            print "The surface and secretory protein reference sequence file is not in the correct format or doesn't exist.\n";
            #print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The surface and secretory protein reference sequence file has not been given.\n";
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
        $outName = `echo "$fastq1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\\+_L[0-9]\\+_R[0-9]\\+.*//g'`;
        print "The default output file name prefix is: $outName";
    }

    return ($help, $fastq1, $fastq2, $ref_dir, $feat_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
GAS_Features_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <ref directory: dir> -p <feature seq: file> -o <output directory name: string> -n <output name prefix: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -d   reference sequence directory (including full path)
    -f   features sequence reference file
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $ref_dir, $feat_DB, $outDir, $outName) = checkOptions( @ARGV );





###Subroutines###
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




###Start Doing Stuff###
chdir "$outDir";
my $feat_out = "TEMP_protein_Results.txt";
open(my $fh,'>',$feat_out) or die "Could not open file '$feat_out' $!";
my $BIN_feat_out = "BIN_Features_Results.txt";
open(my $bh,'>',$BIN_feat_out) or die "Could not open file '$BIN_feat_out' $!";
my @Bin_Feat_arr = (0) x 26;
#print $fh "Feature_Group\tTarget\n";
my $outNameFEAT = "FEAT_".$outName;
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outNameFEAT --log --save_scores --min_coverage 99.9 --max_divergence 8 --gene_db $feat_DB");
my @TEMP_FEAT_fullgene = glob("FEAT_*__fullgenes__*__results\.txt");
my $FEAT_full_name = $TEMP_FEAT_fullgene[0];

my %Feat_Col = (
    "T_Type" => "neg",
    "GACI" => "neg",
    "EMM_Family" => "neg",
    "ECM" => "neg",
    "HASA" => "neg",
    "SDA1" => "neg",
    "SIC" => "neg",
    "ROCA" => "neg",
    "PNGA" => "neg",
    "SLO-G" => "neg",
    "Exotoxins" => "neg",
    "SLAA" => "neg",
    );

###Type the Presence/Absence Targets###
###############################################################################################
open(MYINPUTFILE, "$FEAT_full_name");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    print "$line\n";
    my @feat_fullgene;
    @feat_fullgene = split('\t',$line);
    if ($feat_fullgene[5] >= 10) {
        if ($feat_fullgene[3] =~ m/^T[0-9]+/) {
            if ($Feat_Col{"T_Type"} eq "neg") {
                $Feat_Col{"T_Type"} = $feat_fullgene[2];
		$Bin_Feat_arr[0] = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"T_Type"}.":".$feat_fullgene[2];
                $Feat_Col{"T_Type"} = $new_val;
                $Bin_Feat_arr[0] = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/GACI/) {
	    #$Feat_Col{"GACI"} = $feat_fullgene[2];
	    $Feat_Col{"GACI"} = "pos";
        }
        if ($feat_fullgene[3] =~ m/(MRP|ENN)/) {
            if ($Feat_Col{"EMM_Family"} eq "neg") {
                $Feat_Col{"EMM_Family"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"EMM_Family"}.":".$feat_fullgene[2];
                $Feat_Col{"EMM_Family"} = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/(FBAA|PRTF2|SFB1|R28|SOF)/) {
	    print "FBAA has been found: $feat_fullgene[3]\n";
            if ($Feat_Col{"ECM"} eq "neg") {
                $Feat_Col{"ECM"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"ECM"}.":".$feat_fullgene[2];
                $Feat_Col{"ECM"} = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/HASA/) {
	    $Feat_Col{"HASA"} = $feat_fullgene[2];
	    my $target = "8__HASA__HASA-1__8";
	    my @TEMP_FEAT_bam = glob("FEAT_*\.sorted\.bam");
	    my $bamFile = $TEMP_FEAT_bam[0];
	    (my $samFile = $bamFile) =~ s/\.bam/\.sam/g;
	    system("samtools view -h $bamFile > $samFile");
	    system("cat $samFile | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > HASA_target_seq.sam");
	    system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $samFile >> HASA_target_seq.sam");
	    system("samtools view -bS HASA_target_seq.sam > HASA_target_seq.bam");
	    system("samtools index HASA_target_seq.bam HASA_target_seq.bai");
	    my $REF_seq = extractFastaByID("$target","$feat_DB");
	    open(my $rf,'>',"HASA_target_ref.fna");
	    print $rf "$REF_seq\n";
	    close $rf;
	    system("freebayes -q 20 -p 1 -f HASA_target_ref.fna HASA_target_seq.bam -v HASA_target_seq.vcf");
	    open(MYINPUTFILE2, "HASA_target_seq.vcf");
	    my %srst2_seroT;
	    while(<MYINPUTFILE2>) {
		my $line = $_;
		chomp($line);
		if ($line =~ /^8__HASA__HASA-1__8/) {
		    my @HASA_line = split('\t', $line);
		    my $ref_allele = $HASA_line[3];
		    my $alt_allele = $HASA_line[4];
		    $HASA_line[7] =~ /DPB=(\d+\.?\d*);/;
		    print "HASA DP: $1 | ref allele: $ref_allele | alt allele: $alt_allele\n";
		    my $HASA_dp = $1;
		    my $HASA_loc = $HASA_line[1];
		    if (length($ref_allele) != length($alt_allele) && $HASA_dp >= 2) {
			$Feat_Col{"HASA"} = "neg";
		    }
		}
	    }
	}
        if ($feat_fullgene[3] =~ m/SDA1/) {
            $Feat_Col{"SDA1"} = $feat_fullgene[2];
        }
	if ($feat_fullgene[3] =~ m/SLAA/) {
            $Feat_Col{"SLAA"} = $feat_fullgene[2];
        }
        if ($feat_fullgene[3] =~ m/SIC/) {
            $Feat_Col{"SIC"} = $feat_fullgene[2];
        }
        if ($feat_fullgene[3] =~ m/(^A-1|^C-1|^G-1|^H-1|^I-1|^J-1|^K-1|^L-1|^M-1|^S-1|^Z-1)/) {
            if ($Feat_Col{"Exotoxins"} eq "neg") {
                $Feat_Col{"Exotoxins"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"Exotoxins"}.":".$feat_fullgene[2];
                $Feat_Col{"Exotoxins"} = $new_val;
            }
        }


        if ($feat_fullgene[3] =~ m/(ROCAM18|ROCAM3)/ && !$feat_fullgene[6]) {
	    if ($feat_fullgene[3] =~ m/ROCAM18/) {
		$Bin_Feat_arr[12] = 1;
	    }
	    if ($feat_fullgene[3] =~ m/ROCAM3/) {
		$Bin_Feat_arr[13] = 1;
	    }

            if ($Feat_Col{"ROCA"} eq "neg") {
                $Feat_Col{"ROCA"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"ROCA"}.":".$feat_fullgene[2];
                $Feat_Col{"ROCA"} = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/PNGA/ && !$feat_fullgene[6]) {
	    $Bin_Feat_arr[13] = 1;
            #$Feat_Col{"PNGA"} = $feat_fullgene[2];
	    $Feat_Col{"PNGA"} = "PNGA3";
        }
    }

    if ($feat_fullgene[3] =~ m/SLOG/ && !$feat_fullgene[6] && $feat_fullgene[5] >= 7) {
	$Bin_Feat_arr[14] = 1;
	#$Feat_Col{"SLO-G"} = $feat_fullgene[2];
	$Feat_Col{"SLO-G"} = "330G";
    }
}
###############################################################################################

###############################################################################################
###Make Binary Output Table###
open(MYINPUTFILE, "$FEAT_full_name");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @feat_fullgene;
    @feat_fullgene = split('\t',$line);
    if ($feat_fullgene[5] >= 10) {
	if ($feat_fullgene[3] =~ m/MRP/) {
	    $Bin_Feat_arr[1] = 1;
	}
        if ($feat_fullgene[3] =~ m/ENN/) {
            $Bin_Feat_arr[2] = 1;
        }
        if ($feat_fullgene[3] =~ m/FBAA/) {
            $Bin_Feat_arr[3] = 1;
        }
        if ($feat_fullgene[3] =~ m/PRTF2/) {
            $Bin_Feat_arr[4] = 1;
        }
        if ($feat_fullgene[3] =~ m/SFB1/) {
            $Bin_Feat_arr[5] = 1;
        }
        if ($feat_fullgene[3] =~ m/R28/) {
            $Bin_Feat_arr[6] = 1;
        }
        if ($feat_fullgene[3] =~ m/SOF/) {
            $Bin_Feat_arr[7] = 1;
        }
        if ($feat_fullgene[3] =~ m/HASA/) {
            $Bin_Feat_arr[8] = 1;
        }
        if ($feat_fullgene[3] =~ m/SDA1/) {
            $Bin_Feat_arr[9] = 1;
        }
        if ($feat_fullgene[3] =~ m/SIC/) {
            $Bin_Feat_arr[10] = 1;
        }
        #if ($feat_fullgene[3] =~ m/ROCAM3/) {
        #    $Bin_Feat_arr[11] = 1;
        #}
        #if ($feat_fullgene[3] =~ m/ROCAM18/) {
        #    $Bin_Feat_arr[12] = 1;
        #}
        #if ($feat_fullgene[3] =~ m/PNGA/) {
        #    $Bin_Feat_arr[13] = 1;
        #}
        #if ($feat_fullgene[3] =~ m/SLOG/) {
        #    $Bin_Feat_arr[14] = 1;
        #}
        if ($feat_fullgene[3] =~ m/^A-1/) {
            $Bin_Feat_arr[15] = 1;
        }
        if ($feat_fullgene[3] =~ m/^C-1/) {
            $Bin_Feat_arr[16] = 1;
        }
        if ($feat_fullgene[3] =~ m/^G-1/) {
            $Bin_Feat_arr[17] = 1;
        }
        if ($feat_fullgene[3] =~ m/^H-1/) {
            $Bin_Feat_arr[18] = 1;
        }
        if ($feat_fullgene[3] =~ m/^I-1/) {
            $Bin_Feat_arr[19] = 1;
        }
        if ($feat_fullgene[3] =~ m/^J-1/) {
            $Bin_Feat_arr[20] = 1;
        }
        if ($feat_fullgene[3] =~ m/^K-1/) {
            $Bin_Feat_arr[21] = 1;
        }
        if ($feat_fullgene[3] =~ m/^L-1/) {
            $Bin_Feat_arr[22] = 1;
        }
        if ($feat_fullgene[3] =~ m/^M-1/) {
            $Bin_Feat_arr[23] = 1;
        }
        if ($feat_fullgene[3] =~ m/^S-1/) {
            $Bin_Feat_arr[24] = 1;
        }
        if ($feat_fullgene[3] =~ m/^Z-1/) {
            $Bin_Feat_arr[25] = 1;
        }
    }
}

###Print GAS Binary Output###
print $bh join(',',@Bin_Feat_arr);

###############################################################################################
###Print GAS Features Output###
while (my ($key, $val) = each %Feat_Col) {
    my @val_arr = split(':',$val);
    #print "@val_arr\n";
    my @val_sort = sort { "\L$a" cmp "\L$b" } @val_arr;
    #print "@val_sort\n";
    my $val_out = join(':',@val_sort);
    print "$key\t$val_out\n";
    $Feat_Col{"$key"} = $val_out;
    #print $fh "$key\t$val_out\n";
}

print $fh "T_Type\t$Feat_Col{'T_Type'}\n";
print $fh "GACI\t$Feat_Col{'GACI'}\n";
print $fh "EMM_Family\t$Feat_Col{'EMM_Family'}\n";
print $fh "ECM\t$Feat_Col{'ECM'}\n";
print $fh "HASA\t$Feat_Col{'HASA'}\n";
print $fh "SDA1\t$Feat_Col{'SDA1'}\n";
print $fh "SLAA\t$Feat_Col{'SLAA'}\n";
print $fh "SIC\t$Feat_Col{'SIC'}\n";
print $fh "ROCA\t$Feat_Col{'ROCA'}\n";
print $fh "PNGA\t$Feat_Col{'PNGA'}\n";
print $fh "SLO-G\t$Feat_Col{'SLO-G'}\n";
print $fh "Exotoxins\t$Feat_Col{'Exotoxins'}\n";
###############################################################################################
