#!/bin/bash -l

temp_path=$(pwd)
export PATH=$PATH:$temp_path

## -- begin embedded SGE options --
read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1/job-control.txt)
## -- end embedded SGE options --

###Load Modules###
#. /usr/share/Modules/init/bash
module load perl/5.22.1
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load freebayes/0.9.21
module load prodigal/2.60
#module load cutadapt/1.8
module load cutadapt/1.8.3
module load srst2/0.1.7

###This script is called for each job in the qsub array. The purpose of this code is to read in and parse a line of the job-control.txt file
###created by 'StrepLab-JanOw_GAS-wrapr.sh' and pass that information, as arguments, to other programs responsible for various parts of strain
###characterization (MLST, emm type and antibiotic drug resistance prediction).

readPair_1=${PARAM[0]}
readPair_2=${PARAM[1]}
allDB_dir=${PARAM[2]}
batch_out=${PARAM[3]}
sampl_out=${PARAM[4]}
delete='true'


###Start Doing Stuff###
cd "$sampl_out"
batch_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)}')
out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')  ###Use This For Batches off the MiSeq###
#out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-1)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')   ###Otherwise Use This###
just_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')
out_nameMLST=MLST_"$just_name"
#out_nameEMM=EMM_"$just_name"
#out_namePBP=PBP_"$just_name"

###Pre-Process Paired-end Reads###
fastq1_trimd=cutadapt_"$just_name"_S1_L001_R1_001.fastq
fastq2_trimd=cutadapt_"$just_name"_S1_L001_R2_001.fastq
cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq $readPair_1 $readPair_2
cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output $fastq1_trimd -o $fastq2_trimd temp2.fastq temp1.fastq
rm temp1.fastq
rm temp2.fastq
module load fastqc/0.11.5
mkdir "$just_name"_R1_cut
mkdir "$just_name"_R2_cut
fastqc "$fastq1_trimd" --outdir=./"$just_name"_R1_cut
fastqc "$fastq2_trimd" --outdir=./"$just_name"_R2_cut
module unload fastqc/0.11.5

###Call MLST###
srst2 --samtools_args "\-A" --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_pyogenes.fasta" --mlst_definitions "$allDB_dir/spyogenes.txt" --min_coverage 99.999
###Check and extract new MLST alleles###
MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt "$out_nameMLST"__*.Streptococcus_pyogenes.sorted.bam "$allDB_dir/Streptococcus_pyogenes.fasta"

###Call emm Type###
module unload perl/5.22.1
module load perl/5.16.1-MT
emm_typer.pl -1 "$readPair_1" -2 "$readPair_1" -r "$allDB_dir" -n "$just_name"
PBP-Gene_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/GAS_bLactam_Ref.fasta" -n "$just_name" -s GAS -p 2X
module unload perl/5.16.1-MT
module load perl/5.22.1

###Call GAS Misc Resistance###
GAS_Res_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -d "$allDB_dir" -r GAS_Res_Gene-DB_Final.fasta -n "$just_name"
GAS_Target2MIC.pl TEMP_Res_Results.txt "$just_name" TEMP_pbpID_Results.txt

###Type Surface and Secretory Proteins###
GAS_Features_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -d "$allDB_dir" -f GAS_features_Gene-DB_Final.fasta -n "$just_name"


###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
tabl_out="TABLE_Isolate_Typing_results.txt"
bin_out="BIN_Isolate_Typing_results.txt"
#contamination_level=10
printf "$just_name\t" >> "$tabl_out"
printf "$just_name," >> "$bin_out"
###EMM TYPE OUTPUT###
emm_out="NF"
while read -r line
do
    if [[ -n "$line" ]]
    then
        justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
        if [[ "$emm_out" == "NF" ]]
        then
            emm_out="$justTarget"
        else
            emm_out="$emm_out;$justTarget"
        fi
    fi
done <<< "$(sed 1d *__emm-Type__Results.txt)"
printf "$emm_out\t" >> "$tabl_out"
printf "$emm_out," >> "$bin_out"
###MLST OUTPUT###
sed 1d "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt | while read -r line
do
    MLST_tabl=$(echo "$line" | cut -f2-9)
    echo "MLST line: $MLST_tabl\n";
    printf "$MLST_tabl\t" >> "$tabl_out"
    MLST_val=$(echo "$line" | awk -F" " '{print $2}')
    printf "$MLST_val," >> "$bin_out"
done #< "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt
#tail -n+2 "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt | cut -f2-9 >> "$tabl_out"

###Features Targets###
while read -r line
do
    FEAT_targ=$(echo "$line" | cut -f2)
    printf "$FEAT_targ\t" >> "$tabl_out"
done < TEMP_protein_Results.txt

###PBP_ID Output###
justPBPs="NF"
sed 1d TEMP_pbpID_Results.txt | while read -r line
do
    if [[ -n "$line" ]]
    then
	justPBPs=$(echo "$line" | awk -F"\t" '{print $2}')
    fi
    printf "$justPBPs\t" >> "$tabl_out"
done

###Resistance Targets###
while read -r line
do
    #RES_targ=$(echo "$line" | cut -f2)
    #printf "$RES_targ\t" >> "$tabl_out"
    printf "$line\t" | tr ',' '\t' >> "$tabl_out"
done < RES-MIC_"$just_name"

#if test -n "$(shopt -s nullglob; echo ./velvet_output/*_Logfile.txt)"
if [[ -e $(echo ./velvet_output/*_Logfile.txt) ]]
then
    vel_metrics=$(echo ./velvet_output/*_Logfile.txt)
    print "velvet metrics file: $vel_metrics\n";
    velvetMetrics.pl -i "$vel_metrics";
    line=$(cat velvet_qual_metrics.txt | tr ',' '\t') 
    printf "$line\t" >> "$tabl_out"

    printf "$readPair_1\t" >> "$tabl_out";
    pwd | xargs -I{} echo {}"/velvet_output/contigs.fa" >> "$tabl_out"
else
    printf "NA\tNA\tNA\tNA\t$readPair_1\tNA\n" >> "$tabl_out"
fi
#printf "\n" >> "$tabl_out"

cat BIN_Features_Results.txt | sed 's/$/,/g' >> "$bin_out"
cat BIN_Res_Results.txt >> "$bin_out"
printf "\n" >> "$bin_out"

###Remove Unneeded Files###
if ${delete};
then
    echo "The delete flag is true. Will delete all unneeded files"
    rm *.pileup
    rm *.*sorted.*am
    rm *.scores
    rm *ARGannot_r1.*
    rm *__genes__*
    rm *log
    rm BIN*
    rm *Final.sam*
    rm *Results.txt
    rm emm_region_extract.bed ARG-RESFI_fullgenes_results.txt contig-vs-frwd_nucl.txt velvet_qual_metrics.txt ./velvet_output/Sequences ./velvet_output/stats.txt ./velvet_output/Log
    rm RES-MIC*
    rm emm_v*
    rm cutadapt_*_R*
    rm FOLP_target_*
    rm HASA_target_ref.fna*
    rm HASA_target_seq.[sb]*
    rm TEMP_*
    rm velvet_output/*Graph*
    rm -r *_R1_cut
    rm -r *_R2_cut
else
    echo "The delete flag is false. Will keep all unneeded files"
fi



###Unload Modules###
module unload perl/5.22.1
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload freebayes/0.9.21
module unload prodigal/2.60
#module unload cutadapt/1.8
module unload cutadapt/1.8.3
module unload srst2/0.1.7
