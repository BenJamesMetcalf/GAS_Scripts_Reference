#!/usr/bin/perl -w

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
    getopts('h1:2:r:p:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $ref_dir, $protein_DB, $outDir, $outName);

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
        $ref_dir = $opts{r};
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

    if($opts{p}) {
        $protein_DB = "$ref_dir/$opts{p}";
        if ($protein_DB) {
            print "The surface and secretory protein reference sequence file: $opts{p}\n";
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

    return ($help, $fastq1, $fastq2, $ref_dir, $protein_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
GBS_miscRes_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <ref directory: dir> -m <misc. resistance seq: file> -v <vanc. resistance seq: file> -o <output directory name: string> -n <output name prefix: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   reference sequence directory (including full path)
    -p   surface and secretory protein reference sequence file
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $ref_dir, $protein_DB, $outDir, $outName) = checkOptions( @ARGV );





###Subroutines###




###Start Doing Stuff###
chdir "$outDir";
my $surface_output = "TEMP_protein_Results.txt";
open(my $fh,'>',$surface_output) or die "Could not open file '$surface_output' $!";
print $fh "Target\tMatch_Type\tProtein_Type\tCoverage\n";

#my $outName_PROT = "PROT_".$outName;
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outName --log --save_scores --min_coverage 99.9 --max_divergence 8 --gene_db $protein_DB");
my @TEMP_PROT_fullgene = glob("PROT_*__fullgenes__*__results\.txt");
my $PROT_full_name = $TEMP_PROT_fullgene[0];

my %surface_Targets = (
    "ENN-1" => "M_Like",
    "MRP-1" => "M_Like",
    "FBAA-1" => "Surface",
    "PRFT2-1" => "Surface",
    "R28-1" => "Surface",
    "SFB1-2" => "Surface",
    "SOF-1" => "Surface",
    "SMEZ-1" => "Exotoxin",
    "SPEA-1" => "Exotoxin",
    "SPEC-1" => "Exotoxin",
    "SPEG-1" => "Exotoxin",
    "SPEH-1" => "Exotoxin",
    "SPEI-1" => "Exotoxin",
    "SPEJ-1" => "Exotoxin",
    "SPEK-1" => "Exotoxin",
    "SPEL-1" => "Exotoxin",
    "SPEM-1" => "Exotoxin",
    "SSA-1" => "Exotoxin",
    "T1-1" => "T_Antigen",
    "T11-1" => "T_Antigen",
    "T12-1" => "T_Antigen",
    "T2-1" => "T_Antigen",
    "T25-1" => "T_Antigen",
    "T28-1" => "T_Antigen",
    "T3-1" => "T_Antigen",
    "T4-1" => "T_Antigen",
    "T5-1" => "T_Antigen",
    "T6-1" => "T_Antigen",
    "T9-1" => "T_Antigen",
    "HASA-1" => "Capsule",
    "SDA1-1" => "Exotoxin",
    "SIC-1" => "Exotoxin",
    );

my $isNotResistant = "yes";
###Type the non-vancomycin resistance targets###
my %protein_extract;
open(MYINPUTFILE, "$PROT_full_name");
my %protein_Type;
my @protein_fullgene;
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    @protein_fullgene = split('\t',$line);
    my $isMismatch = "no";
    if ($protein_fullgene[6]) {
        $isMismatch = "yes";
    }
    my $value_out = "$protein_fullgene[5]:$isMismatch";
    $protein_Type{$protein_fullgene[3]} = $value_out;
}
print Dumper(\%protein_Type);

foreach my $key (keys(%protein_Type)) {
    #print "$miscR_Type{$key}\n";
    if ($surface_Targets{$key}) {
        my @protein_value = split(':',$protein_Type{$key});
        my $status = "identical";
        if ($protein_value[1] eq "yes") {
            $status = "imperfect";
        }
        print $fh "$key\t$status\t$surface_Targets{$key}\t$protein_value[0]\n";
        $isNotResistant = "no";
    }
}

if ($protein_Type{"PNGA-1"}) {
    my @protein_value = split(':',$protein_Type{"PNGA-1"});
    if ($protein_value[1] eq "no") {
        print $fh "PNGA-1\tidentical\tPromoter\t$protein_value[0]\n";
    }
}

if ($protein_Type{"ROCAM18-1"}) {
    my @protein_value = split(':',$protein_Type{"ROCAM18-1"});
    if ($protein_value[1] eq "no") {
        print $fh "ROCAM18-1\tidentical\tRocA\t$protein_value[0]\n";
    }
}

if ($protein_Type{"ROCAM3-1"}) {
    my @protein_value = split(':',$protein_Type{"ROCAM3-1"});
    if ($protein_value[1] eq "no") {
        print $fh "ROCAM3-1\tidentical\tRocA\t$protein_value[0]\n";
    }
}

###Delete TEMP files and close open files###
#unlink(TEMP*);
close $fh;
