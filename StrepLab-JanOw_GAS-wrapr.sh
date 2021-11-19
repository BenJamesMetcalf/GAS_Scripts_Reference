#!/bin/bash -l

#. /usr/share/Modules/init/bash
###This wrapper script validates the input arguments and creates the job-control.txt file which is needed to submit the qsub array job to the cluster.###

while getopts :s:r:o: option
do
    case $option in
        s) batch_dir=$OPTARG;;
        r) allDB_dir=$OPTARG;;
        o) output_dir=$OPTARG;;
    esac
done

###Check if batch directory and reference database directory arguments were given and if they exist###
if [[ ! -z "$batch_dir" ]]
then
    if [[ -d "$batch_dir" ]]
    then
	batch_dir=$(echo "$batch_dir" | sed 's/\/$//g')
	echo "The sequence directory is in the following location: $batch_dir"
    else
	echo "This sequence directory is not in the correct format or doesn't exist."
	echo "Make sure you provide the full directory path (/root/path/sequence_directory)."
	exit 1
    fi
else
    echo "No sequence data directory path argument given."
    exit 1
fi

if [[ ! -z "$allDB_dir" ]]
then
    if [[ -d "$allDB_dir" ]]
    then
	allDB_dir=$(echo "$allDB_dir" | sed 's/\/$//g')
	echo "The references directory is in the following location: $allDB_dir"
    else
	echo "This reference directory is not in the correct format or doesn't exist."
	echo "Make sure you provide the full directory path (/root/path/reference_directory)."
	exit 1
    fi
else
    echo "No reference database directory path argument given."
    exit 1
fi

###Check if the output directory argument has been given. If yes, create the 'GAS_Typing_Output' and 'qsub_files' folders within the output dir###
###If no, output the results into a subdirectory of '~/GAS_Typing_Analysis'. The subdirectory name is extracted from the batch sequence full path###
if [[ -z "$output_dir" ]]
then
    echo "The files will be output into the default directory 'GAS_Typing_Analysis'."
    if [[ ! -d ~/GAS_Typing_Analysis ]]
    then
        mkdir ~/GAS_Typing_Analysis
        out_dir="~/GAS_Typing_Analysis"
        eval out_dir=$out_dir
        echo "The output directory has been created: $out_dir"
    else
        out_dir="~/GAS_Typing_Analysis"
        eval out_dir=$out_dir
    fi
    batch_name=$(echo "$batch_dir" | awk -F"/" '{print $(NF-3)}')
    out_analysis="${out_dir}"/"${batch_name}"/GAS_Typing_Output
    out_qsub="${out_dir}"/"${batch_name}"/qsub_files/
    out_jobCntrl="${out_dir}/${batch_name}/"
    eval out_analysis=$out_analysis
    eval out_qsub=$out_qsub
    eval out_jobCntrl=$out_jobCntrl
    mkdir -p "$out_analysis"
    mkdir -p "$out_qsub"
elif [[ ! -d "$output_dir" ]]
then
    output_dir=$(echo "$output_dir" | sed 's/\/$//g')
    mkdir "$output_dir"
    out_dir="$output_dir"
    eval out_dir=$out_dir
    echo "The output directory has been created: $out_dir"
    out_analysis="${out_dir}"/GAS_Typing_Output
    out_qsub="${out_dir}"/qsub_files/
    out_jobCntrl="${out_dir}/"
    eval out_analysis=$out_analysis
    eval out_qsub=$out_qsub
    eval out_jobCntrl=$out_jobCntrl
    mkdir -p "$out_analysis"
    mkdir -p "$out_qsub"
else
    output_dir=$(echo "$output_dir" | sed 's/\/$//g')
    out_dir="$output_dir"
    eval out_dir=$out_dir
    out_analysis="${out_dir}"/GAS_Typing_Output
    out_qsub="${out_dir}"/qsub_files/
    out_jobCntrl="${out_dir}/"
    eval out_analysis=$out_analysis
    eval out_qsub=$out_qsub
    eval out_jobCntrl=$out_jobCntrl
    mkdir -p "$out_analysis"
    mkdir -p "$out_qsub"
fi

###Create the batch output files###
batch_name=$(echo "$batch_dir" | awk -F"/" '{print $(NF-3)}')
#printf "Sample_Name\temm_Type\temm_Seq\t%_identity\tmatch_length\n" >> "$out_analysis"/JanOw_"$batch_name"_emmType_results.txt
printf "Sample\temm_Type\tST\tgki\tgtr\tmurI\tmutS\trecP\txpt\tyqiL\tT_Type\tGroup_A\tEMM_Family\tOther_Surface_Proteins\tCapsule\tSDA1\tSLAA\tSIC\tROCA\tPNGA3\tNADase_D330G\tExotoxins\tPBP_ID\tWGS_ZOX_SIGN\tWGS_ZOX\tWGS_ZOX_SIR\tWGS_FOX_SIGN\tWGS_FOX\tWGS_FOX_SIR\tWGS_TAX_SIGN\tWGS_TAX\tWGS_TAX_SIR\tWGS_CFT_SIGN\tWGS_CFT\tWGS_CFT_SIR\tWGS_CPT_SIGN\tWGS_CPT\tWGS_CPT_SIR\tWGS_CZL_SIGN\tWGS_CZL\tWGS_CZL_SIR\tWGS_AMP_SIGN\tWGS_AMP\tWGS_AMP_SIR\tWGS_PEN_SIGN\tWGS_PEN\tWGS_PEN_SIR\tWGS_MER_SIGN\tWGS_MER\tWGS_MER_SIR\tER_CL\tWGS_ERY_SIGN\tWGS_ERY\tWGS_ERY_SIR\tWGS_CLI_SIGN\tWGS_CLI\tWGS_CLI_SIR\tWGS_LZO_SIGN\tWGS_LZO\tWGS_LZO_SIR\tWGS_SYN_SIGN\tWGS_SYN\tWGS_SYN_SIR\tWGS_ERY/CLI\tTET\tWGS_TET_SIGN\tWGS_TET\tWGS_TET_SIR\tGYRA_PARC\tWGS_CIP_SIGN\tWGS_CIP\tWGS_CIP_SIR\tWGS_LFX_SIGN\tWGS_LFX\tWGS_LFX_SIR\tOTHER\tWGS_DAP_SIGN\tWGS_DAP\tWGS_DAP_SIR\tWGS_VAN_SIGN\tWGS_VAN\tWGS_VAN_SIR\tWGS_RIF_SIGN\tWGS_RIF\tWGS_RIF_SIR\tWGS_CHL_SIGN\tWGS_CHL\tWGS_CHL_SIR\tWGS_SXT_SIGN\tWGS_SXT\tWGS_SXT_SIR\tContig_Num\tN50\tLongest_Contig\tTotal_Bases\tReadPair_1\tContig_Path\n" >> "$out_analysis"/TABLE_GAS_"$batch_name"_Typing_Results.txt
#printf "Sample,MLST,emm_Type,T_Type,MRP,ENN,FBAA,PRTF2,SFB1,R28,SOF,HASA,SDA1,SIC,ROCAM3,ROCAM18,PNGA,SLOG,SpeA,SpeC,SpeG,SpeH,SpeI,SpeJ,SpeK,SpeL,SpeM,SSA,SMEZ,23S1,23S3,CAT,ERMB,ERMT,ERMA,FOLA,FOLP1,FOLP2,GYRA,LNUB,LSAC,LSAE,MEF,PARC,RPOB1,RPOBN,TETL,TETM,TETO\n" >> "$out_analysis"/BIN_GAS_"$batch_name"_Typing_Results.txt

###Will search thru every file in the batch directory and check if it matches the following regexs: _L.*_R1_001.fastq and _L.*_R2_001.fastq###
###If both paired end fastq files are found then the full paths of each file will be written to the 'job-control.txt' file###
batch_dir_star="${batch_dir}/*"
for sample in $batch_dir_star
do
    sampl_name=$(echo "$sample" | sed 's/^.*\///g' | sed 's/_S[0-9]\+_.*_001.fastq.gz//g')
    sampl_out="${out_analysis}"/"${sampl_name}"
    eval sampl_out=$sampl_out
    echo The sample file is: $sample
    if [[ $sampl_name =~ ^Undetermined ]]
    then
	echo "Skipping the 'Undetermined' fastq files"
	continue
    fi

    if [[ $sample =~ _L.*_R1_001.fastq && ! $sample =~ S[0-9]+ ]]
    then
        readPair_1=$(echo "$sample" | sed 's/_L\([0-9]\+\)_R1/_S1_L\1_R1/g')
        mv $sample $readPair_1
    elif [[ $sample =~ _L.*_R1_001.fastq && $sample =~ S[0-9]+ ]]
    then
        readPair_1=$sample
    fi

    if [[ $sample =~ _L.*_R2_001.fastq && ! $sample =~ S[0-9]+ ]]
    then
        readPair_2=$(echo "$sample" | sed 's/_L\([0-9]\+\)_R2/_S1_L\1_R2/g')
        mv $sample $readPair_2
    elif [[ $sample =~ _L.*_R2_001.fastq && $sample =~ S[0-9]+ ]]
    then
        readPair_2=$sample
    fi

    if [ -n "$readPair_1" -a -n "$readPair_2" ]
    then
        if [[ ! -d "$sampl_out" ]]
        then
            mkdir "$sampl_out"
        fi
        echo "Both Forward and Reverse Read files exist."
        echo "Paired-end Read-1 is: $readPair_1"
        echo "Paired-end Read-2 is: $readPair_2"
        printf "\n"
        echo "$readPair_1 $readPair_2 $allDB_dir $out_analysis $sampl_out" >> $out_jobCntrl/job-control.txt
        ###Prepare script for next sample###
        readPair_1=""
        readPair_2=""
    fi
done

qsub -sync y -q short.q -t 1-$(cat $out_jobCntrl/job-control.txt | wc -l) -cwd -o "$out_qsub" -e "$out_qsub" ./StrepLab-JanOw_GAS-Typer.sh $out_jobCntrl

###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
while read -r line
do
    batch_name=$(echo $line | awk -F" " '{print $1}' | awk -F"/" '{print $(NF-4)}')
    final_outDir=$(echo $line | awk -F" " '{print $5}')
    final_result_Dir=$(echo $line | awk -F" " '{print $4}')
    #cat $final_outDir/SAMPLE_Isolate__Typing_Results.txt >> $final_result_Dir/SAMPL_GAS_"$batch_name"_Typing_Results.txt
    cat $final_outDir/TABLE_Isolate_Typing_results.txt >> $final_result_Dir/TABLE_GAS_"$batch_name"_Typing_Results.txt
    rm $final_outDir/TABLE_Isolate_Typing_results.txt
    #cat $final_outDir/BIN_Isolate_Typing_results.txt >> $final_result_Dir/BIN_GAS_"$batch_name"_Typing_Results.txt
    #cat $final_outDir/TEMP_newPBP_allele_info.txt >> $final_result_Dir/UPDATR_GBS_"$batch_name"_Typing_Results.txt
    if [[ -e $final_outDir/newPBP_allele_info.txt ]]
    then
        cat $final_outDir/newPBP_allele_info.txt >> $final_result_Dir/UPDATR_GBS_"$batch_name"_Typing_Results.txt
    fi
done < $out_jobCntrl/job-control.txt
