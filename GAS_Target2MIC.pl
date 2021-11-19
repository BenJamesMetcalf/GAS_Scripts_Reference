#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;
use File::Copy qw(copy);
use Env;
#use lib $ENV{MODULESHOME}."/init";
use lib "/usr/share/Modules/init/";
use perl;

###Start Doing Stuff###
my $Res_output = "RES-MIC_".$ARGV[1];
open(my $fh,'>',$Res_output) or die "Could not open file '$Res_output' $!";
#print Dumper \@ARGV;
print "Output file name is: $Res_output\n";

my %Res_hash;
my $RES_full_name = $ARGV[0];
open(MYINPUTFILE, "$RES_full_name");
while(<MYINPUTFILE>) {
    #next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @res_arr;
    @res_arr = split('\t',$line);
    $Res_hash{$res_arr[0]} = $res_arr[1];
}
close MYINPUTFILE;

my $PBP_full_name = $ARGV[2];
my $pbp_type = "NA";
open (MYINPUTFILE, "$PBP_full_name");
while (<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @pbp_arr;
    @pbp_arr = split('\t',$line);
    $pbp_type = $pbp_arr[1];
    print "PBP Type: $pbp_type\n";
}
close MYINPUTFILE;

#while (my ($key, $val) = each %Res_hash) {
#    my @val_arr = split(':',$val);
#    my @val_sort = sort(@val_arr);
#    my $val_out = join(':',@val_sort);
#    print "$key\t$val_out\n";
#}
#print "\n";
print Dumper \%Res_hash;

my %drug;
my %Out_hash;
###TET Category###
if ($Res_hash{"TET"} eq "neg") {
    printf "TET,neg,<=,2,S\n";
    $Out_hash{TET} = "neg,<=,2,S";
} else {
    my @Res_targs = split(':',$Res_hash{TET});
    if ( grep( /TET/i, @Res_targs ) ) {
	printf "TET,$Res_hash{TET},>=,8,R\n";
	$Out_hash{TET} = "$Res_hash{TET},>=,8,R";
    }
}

###PBP Category###
$drug{ZOX} = "Flag,Flag,Flag";
$drug{FOX} = "Flag,Flag,Flag";
$drug{TAX} = "Flag,Flag,Flag";
$drug{CZL} = "Flag,Flag,Flag";
$drug{CFT} = "Flag,Flag,Flag";
$drug{CPT} = "Flag,Flag,Flag";
$drug{AMP} = "Flag,Flag,Flag";
$drug{PEN} = "Flag,Flag,Flag";
$drug{MER} = "Flag,Flag,Flag";
if ($pbp_type <= 18 && $pbp_type != [11|17] && $pbp_type !~ /[A-Za-z]/) {
    $drug{ZOX} = "<=,0.12,U";
    $drug{FOX} = "<=,2,U";
    $drug{TAX} = "<=,0.06,S";
    $drug{CFT} = "<=,0.06,S";
    $drug{CPT} = "<=,0.06,S";
    $drug{CZL} = "<=,0.5,U";
    $drug{AMP} = "<=,0.06,S";
    $drug{PEN} = "<=,0.06,S";
    $drug{MER} = "<=,0.06,S";
}
print "PBP,$drug{ZOX},$drug{FOX},$drug{TAX},$drug{CFT},$drug{CPT},$drug{CZL},$drug{AMP},$drug{PEN},$drug{MER}\n";
$Out_hash{PBP} = "$drug{ZOX},$drug{FOX},$drug{TAX},$drug{CFT},$drug{CPT},$drug{CZL},$drug{AMP},$drug{PEN},$drug{MER}";


###ER_CL Category###
$drug{ERY} = "<=,0.25,S";
$drug{CLI} = "<=,0.25,S";
$drug{LZO} = "<=,2.0,S";
$drug{SYN} = "<=,1,S";
$drug{ERY_SYN} = "neg";
if ($Res_hash{"ER_CL"} eq "neg") {
    print "ER_CL,$Res_hash{ER_CL},".$drug{ERY}.",".$drug{CLI}.",".$drug{LZO}.",".$drug{SYN}.",".$drug{ERY_SYN}."\n";
    $Out_hash{ER_CL} = "$Res_hash{ER_CL},$drug{ERY},$drug{CLI},$drug{LZO},$drug{SYN},$drug{ERY_SYN}";
} else {
    my @Res_targs = split(':',$Res_hash{ER_CL});
    if (grep(/23S/i,@Res_targs) || grep(/RPLD/i,@Res_targs)) {
	$drug{ERY} = "Flag,Flag,Flag";
	$drug{CLI} = "Flag,Flag,Flag";
	$drug{LZO} = "Flag,Flag,Flag";
	$drug{SYN} = "Flag,Flag,Flag";
	$drug{ERY_SYN} = "Flag";
    }
    if (grep(/MEF/i,@Res_targs)) {
        print "Found MEF\n";
	$drug{ERY} = ">=,1,R";
    }
    if (grep(/LSA/i,@Res_targs)) {
	print "Found LSA\n";
	$drug{CLI} = "Flag,Flag,Flag";
    }
    if (grep(/ERM/i,@Res_targs)) {
	print "Found ERM\n";
	$drug{ERY} = ">=,1,R";
	$drug{CLI} = ">=,1,R";
	$drug{ERY_SYN} = "pos";
    }
    if (grep(/ERM/i,@Res_targs) && grep(/LSA/i,@Res_targs)) {
	print "Found ERM and LSAC\n";
        $drug{ERY} = ">=,1,R";
        $drug{CLI} = ">=,1,R";
	$drug{SYN} = "Flag,Flag,Flag";
	$drug{ERY_SYN} = "pos";
    }
    print "ER_CL,$Res_hash{ER_CL},".$drug{ERY}.",".$drug{CLI}.",".$drug{LZO}.",".$drug{SYN}.",".$drug{ERY_SYN}."\n";
    $Out_hash{ER_CL} = "$Res_hash{ER_CL},$drug{ERY},$drug{CLI},$drug{LZO},$drug{SYN},$drug{ERY_SYN}";
}

###GYRA_PARC Category###
$drug{LFX} = "<=,2,S";
$drug{CIP} = "NA,NA,NA";
if ($Res_hash{"GYRA_PARC"} eq "neg") {
    print "GYRA_PARC,$Res_hash{GYRA_PARC},".$drug{LFX}."\n";
    $Out_hash{"GYRA_PARC"} = "$Res_hash{GYRA_PARC},$drug{CIP},$drug{LFX}";
} else {
    my @Res_targs = split(':',$Res_hash{GYRA_PARC});
    if ((grep/GYRA-S11F/i,@Res_targs) && (grep/PARC-S6[F|Y]/i,@Res_targs)) {
	print "Found GYRA-S11F:PARC-S6[F|Y]\n";
	$drug{LFX} = ">=,8,R";
    } elsif ((grep/PARC-S6F/i,@Res_targs)) {
	print "Found PARC-S6F\n";
	$drug{LFX} = "=,4,I";
    } elsif ((grep/PARC-D10[G|N]/i,@Res_targs) || (grep/PARC-D5N/i,@Res_targs) || (grep/PARC-S6A/i,@Res_targs)) {
	print "Found PARC-D10[G|N] or PARC-D5N or PARC-S6A\n";
	$drug{LFX} = "=,2,S";
    } else {
	$drug{LFX} = "Flag,Flag,Flag";
    }
    print "GYRA_PARC,$Res_hash{GYRA_PARC},".$drug{LFX}."\n";
    $Out_hash{"GYRA_PARC"} = "$Res_hash{GYRA_PARC},$drug{CIP},$drug{LFX}";
}

###Other Category###
$drug{DAP} = "<=,1,S";
$drug{VAN} = "<=,1,S";
$drug{RIF} = "<=,1,U";
$drug{CHL} = "<=,4,S";
$drug{SXT} = "<=,0.5,U";
if ($Res_hash{"OTHER"} eq "neg") {
    print "OTHER,$Res_hash{OTHER},".$drug{DAP}.",".$drug{VAN}.",".$drug{RIF}.",".$drug{CHL}.",".$drug{SXT}."\n";
    $Out_hash{"OTHER"} = "$Res_hash{OTHER},$drug{DAP},$drug{VAN},$drug{RIF},$drug{CHL},$drug{SXT}";
} else {
    my @Res_targs = split(':',$Res_hash{OTHER});
    if ($Res_hash{OTHER} ne "ant(6)-Ia:Ant6-Ia_AGly:aph(3')-III:Aph3-III_AGly:Sat4A_Agly" && $Res_hash{OTHER} ne "aph(3')-III:Aph3-III_AGly:Sat4A_AGly") {
	print "I'm in the wrong loop\n";
	foreach my $target (@Res_targs) {
	    if ($target !~ m/CAT|FOLA|FOLP|RPOB|MSR/i) { #&& $Res_hash{"OTHER"} !~  m/MSRD/i) {
		print "Found an ARGANNOT/RESFINDER target. Flag everything\n";
		$drug{DAP} = "Flag,Flag,Flag";
		$drug{VAN} = "Flag,Flag,Flag";
		$drug{RIF} = "Flag,Flag,Flag";
		$drug{CHL} = "Flag,Flag,Flag";
		$drug{SXT} = "Flag,Flag,Flag";
		last;
	    }
	}
    }
	foreach my $target (@Res_targs) {
	    if ($target =~ m/CAT/i) {
		print "Found CAT\n";
		$drug{CHL} = ">=,16,R";
	    } elsif ($target =~ m/FOLA/i && $target =~ m/FOLP/i) {
		print "Found FOLA and FOLP\n";
		$drug{SXT} = ">=,4,U";
	    } elsif ($target =~ m/FOLA/i || $target =~ m/FOLP/i) {
		print "Found FOLA or FOLP\n";
		$drug{SXT} = "=,2,U";
	    } elsif ($target =~ m/RPOB/i) {
		print "Found RPOB\n";
		$drug{RIF} = "Flag,Flag,Flag";
	    }
	}
    print "OTHER,$Res_hash{OTHER},".$drug{DAP}.",".$drug{VAN}.",".$drug{RIF}.",".$drug{CHL}."\n";
    $Out_hash{"OTHER"} = "$Res_hash{OTHER},$drug{DAP},$drug{VAN},$drug{RIF},$drug{CHL},$drug{SXT}";
}

print $fh $Out_hash{PBP}.",".$Out_hash{ER_CL}.",".$Out_hash{TET}.",".$Out_hash{GYRA_PARC}.",".$Out_hash{OTHER}."\n";
print "$Out_hash{PBP}||$Out_hash{ER_CL}||$Out_hash{TET}||$Out_hash{GYRA_PARC}||$Out_hash{OTHER}\n";
