#!/usr/bin/perl 
# This script will convert your DNA sequence to PROTEIN Sequence in 6 frames. Each frame will be saved seperately in a file.
# Provide sequence file containing sequences in fasta format at the command line or when asked.
# Author: Ratnesh Singh
# Report bugs/errors at : ratnesh@hawaii.edu

use strict;
use warnings;
use Getopt::Std;

our($opt_s,$opt_r,$opt_f,$opt_l,$opt_c);
#print "\n\n\t\#################### DNA 2 PROTEIN #################### \n\n";

$opt_r="o";
$opt_f=1;
$opt_l=5;
$opt_c='no';

getopt('sfrlc');

my$usage="This script will convert your DNA sequence to PROTEIN Sequence in 6 frames\n
Usage: perl script -options
-ssequence_file
-fframe [1]
-r[o]nly || [r]ange [o]
-lsmallest_protein-allowed [5]
-cyes|no Clean sequence of stop codon in frame 1.[no]
";

die "$usage\nSequence file not found\n" if !$opt_s;

open(DNAFILE, $opt_s) or die print "Cannot open file \"$opt_s\"\n\n";
open OUT,">$opt_s.TranslatedProtein.fasta";
#open OUT2,">$opt_s.ORF.out.fasta";
#open OUT3,">$opt_s.dna.fasta";
my $range=uc$opt_r;
chomp $range;
$/="\n>";

while(<DNAFILE>){

    (my $DNAheader,my @DANseq)=split(/\n/,$_);
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

    #####################################
    my @protein='';
    my @dna='';
    my $codon1;

    for(my $i=0;$i<(length($DNAseq)-2);$i+=3){
	$codon1=substr($DNAseq,$i,3);
	$protein[1].= codon2aa($codon1);
	$dna[1].=codon2nt($codon1);
    }
    #####################################
    #my $protein2='';
    my $codon2;
    for(my $i=1;$i<(length($DNAseq)-2);$i+=3){
	$codon2=substr($DNAseq,$i,3);
	$protein[2].= codon2aa($codon2);
	$dna[2].=codon2nt($codon2);
    }
    #####################################
    #my $protein3='';
    my $codon3;
    for(my $i=2;$i<(length($DNAseq)-2);$i+=3){
	$codon3=substr($DNAseq,$i,3);
	$protein[3].= codon2aa($codon3);
	$dna[3].=codon2nt($codon3);
    }
    #####################################
    #my $protein4='';
    my $codon4;
    for(my $i=0;$i<(length($DNArevSeq)-2);$i+=3){
	$codon4=substr($DNArevSeq,$i,3);
	$protein[4].= codon2aa($codon4);
	$dna[4].=codon2nt($codon4);
    }
    #####################################
    #my $protein5='';
    my $codon5;
    for(my $i=1;$i<(length($DNArevSeq)-2);$i+=3){
	$codon5=substr($DNArevSeq,$i,3);
	$protein[5].= codon2aa($codon5);
	$dna[5].=codon2nt($codon5);
    }
    #####################################
    #my $protein6='';
    my $codon6;
    for(my $i=2;$i<(length($DNArevSeq)-2);$i+=3){
	$codon6=substr($DNArevSeq,$i,3);
	$protein[6].= codon2aa($codon6);
	$dna[6].=codon2nt($codon6);
    }
    #####################################
    #print protein in requested frames

    if(uc$range eq 'O'){
	my $newframe;
	if($opt_f==1){$newframe='+1';}
	if($opt_f==2){$newframe='+2';}
	if($opt_f==3){$newframe='+3';}
	if($opt_f==4){$newframe='-1';}
	if($opt_f==5){$newframe='-2';}
	if($opt_f==6){$newframe='-3';}
	#print OUT ">$DNAheader".'_'."$newframe\n$protein[$opt_f]\n";
	#print ">$DNAheader".'_'."$newframe\n$protein[$opt_f]\n";

        print OUT ">$DNAheader\n$protein[$opt_f]\n";
        print ">$DNAheader\n$protein[$opt_f]\n";	
	
	#print OUT3 ">$DNAheader".'_'."$newframe\n$dna[$opt_f]\n";
	#print ">$DNAheader".'_'."$newframe\n$dna[$opt_f]\n";
	
    }
    else{for(my $i=1;$i<=$opt_f;$i++){
	my $newframe;
	if($i==1){$newframe='+1';}
	if($i==2){$newframe='+2';}
	if($i==3){$newframe='+3';}
	if($i==4){$newframe='-1';}
	if($i==5){$newframe='-2';}
	if($i==6){$newframe='-3';}

	#print OUT">$DNAheader".'_'."$newframe\n$protein[$i]\n";
	#print ">$DNAheader".'_'."$newframe\n$protein[$i]\n";
	
	#print OUT3 ">$DNAheader".'_'."$newframe\n$dna[$i]\n";
	#print ">$DNAheader".'_'."$newframe\n$dna[$i]\n";

	# Selecting ORFs from translated protein
	my (@protein_ORF)=split(/[\*]+/,$protein[$i]);
	my $protein_ORF=@protein_ORF;
	my $k=0;
	for(my$j=0;$j<=$protein_ORF;$j++){
	    next if !defined $protein_ORF[$j] ;
	    $protein_ORF[$j]=~s/\s*//g;
	    if ($protein_ORF[$j]=~/M/){$protein_ORF[$j]=~s/(\w*?M)/M/;} 
	    else{next;}
	    if(length$protein_ORF[$j]<$opt_l){next;}
	    $k++;
	    my $prot_length=length$protein_ORF[$j];
	    #print OUT2">$DNAheader".'_'."$newframe\.$k\t$prot_length aa\n$protein_ORF[$j]\n" if defined $protein_ORF[$j];
	    #print ">$DNAheader".'_'."$newframe\.$k\t$prot_length aa\n$protein_ORF[$j]\n" if defined $protein_ORF[$j];
	}
	 }
    }
}

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
    else{ return('x');
	  
	  print "Bad codon \"$codon\"!!\n";
    }
}
                             
sub codon2nt{
    my($codon)=@_;
    $codon=uc $codon;
    my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
    
    if(lc$opt_c eq 'yes'){
	if($g{$codon} ne '*'){return $codon;}
	else{ return('---');}
    }
    else{return $codon}
}
