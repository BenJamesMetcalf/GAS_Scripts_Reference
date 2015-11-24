#!/bin/bash -l

## -- begin embedded SGE options --
read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1/job-control.txt)
## -- end embedded SGE options --

###Load Modules###
#. /usr/share/Modules/init/bash
module load perl/5.12.3
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load Python/2.7
module load samtools/0.1.18
module load bowtie2/2.1.0
module load freebayes/0.9.21
module load prodigal/2.60

###This script is called for each job in the qsub array. The purpose of this code is to read in and parse a line of the job-control.txt file
###created by 'StrepLab-JanOw_GAS-wrapr.sh' and pass that information, as arguments, to other programs responsible for various parts of strain
###characterization (MLST, emm type and antibiotic drug resistance prediction).

readPair_1=${PARAM[0]}
readPair_2=${PARAM[1]}
allDB_dir=${PARAM[2]}
batch_out=${PARAM[3]}
sampl_out=${PARAM[4]}
#mlst_ref=${PARAM[5]}
#mlst_def=${PARAM[6]}




###Start Doing Stuff###
cd "$sampl_out"
batch_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)}')
out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')  ###Use This For Batches off the MiSeq###
#out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-1)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')   ###Otherwise Use This###
just_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')
out_nameMLST=MLST_"$just_name"
out_nameEMM=EMM_"$just_name"
out_namePROT=PROT_"$just_name"
out_nameMISC=MISC_"$just_name"

###Call MLST###
mod-srst2.py --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_pyogenes.fasta" --mlst_definitions "$allDB_dir/spyogenes.txt"
###Check and extract new MLST alleles###
perl ~/TEMP_GAS-Typing/MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt "$out_nameMLST"__*.Streptococcus_pyogenes.sorted.bam "$allDB_dir/Streptococcus_pyogenes.fasta"

###Call emm Type###
perl ~/TEMP_GAS-Typing/emm_typer.pl -1 "$readPair_1" -2 "$readPair_1" -r "$allDB_dir" -n "$out_nameEMM"

###Call GAS Misc Resistance###
perl ~/TEMP_GAS-Typing/GAS_miscRes_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -m GAS_miscR_Gene-DB_Final.fasta -v GAS_vancR_Gene-DB_Final.fasta -n "$just_name"

###Type Surface and Secretory Proteins###
perl ~/TEMP_GAS-Typing/GAS_surface-secretory_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -p GAS_protein_Gene-DB_Final.fasta -n "$out_namePROT"

#Add in perl script to find contamination threshold here 
contamination_level=10





###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
sampl_out="TEMP_GAS_Typing_Results.txt"

printf "$just_name\n" >> "$sampl_out"
printf "$just_name\t" >> TEMP_table_results.txt
###EMM TYPE OUTPUT###
printf "\tEmm_Type:\n" >> "$sampl_out"
lineNum=$(cat "$out_nameEMM"__emm-Type__Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no emm type was found
    firstLine=$(head -n1 Serotype_results.txt)
    printf "\t\t$firstLine\n\t\tNo_emm_Type\n" >> "$sampl_out"
    printf "No_emm_Type\t" >> TEMP_table_results.txt
else
    count=0
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
            printf "\t\t$line\n" >> "$sampl_out"
        else
            printf "\t\t$line\n" >> "$sampl_out"
            justTarget=$(echo "$line" | awk -F"\t" '{print $2}')
            printf "$justTarget;" >> TEMP_table_results.txt
        fi
    done < "$out_nameEMM"__emm-Type__Results.txt
fi
printf "\t" >> TEMP_table_results.txt

###T-ANTIGEN TYPE OUTPUT###
printf "\tT-Antigen:\n" >> "$sampl_out"
lineNum=$(cat TEMP_protein_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no T-antigens were typed
    firstLine=$(head -n1 TEMP_protein_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_T-Antigen_Targets\n" >> "$sampl_out"
    printf "No_T-Antigen_Targets" >> TEMP_table_results.txt
else
    count=0
    Tantigen_target=()
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
	    #print the header line to both output files
            printf "\t\t$line\n" >> "$sampl_out"
	    printf "$line\n" >> TEMP_no-Tantigen_Results.txt
        else
	    target_type=$(echo "$line" | awk -F"\t" '{print $3}')
	    if [[ "$target_type" == "T_Antigen" ]]
	    then
		#if target is a 'T_Antigen' and if the depth is greater than the contamination level 
	        #then add that T-type to the output array 'Tantigen_target'
		printf "\t\t$line\n" >> "$sampl_out"
		justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
		justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
		#printf "$justTarget-($justDepth);" >> TEMP_table_results.txt
		if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
		then
		    echo "Target $justTarget is a match"
		    Tantigen_target+=("$justTarget($justDepth)")
		fi
	    else
		#if target is not a 'T_Antigen' then print it out to the non-T Antigen protein typing file
		printf "$line\n" >> TEMP_no-Tantigen_Results.txt
	    fi
        fi
    done < TEMP_protein_Results.txt
    #if the output array 'Tantigen_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [[ ${#Tantigen_target[@]} -eq 0 ]] 
    then
	printf "No_T-Antigen_Targets" >> TEMP_table_results.txt
    else
	printf '%s\n' "${Tantigen_target[@]}" | sort | tr '\n' ';'
	printf '%s\n' "${Tantigen_target[@]}" | sort | tr '\n' ';' >> TEMP_table_results.txt
    fi
fi
printf "\t" >> TEMP_table_results.txt

###MLST OUTPUT###
printf "\tMLST:\n" >> "$sampl_out"
count=0
while read -r line
do
    count=$(( $count + 1 ))
    if [[ "$count" -eq 1 ]]
    then
        printf "\t\t$line\n" >> "$sampl_out"
    else
        printf "\t\t$line\n" >> "$sampl_out"
        MLST_tabl=$(echo "$line" | cut -f2-9)
        printf "$MLST_tabl\t" >> TEMP_table_results.txt
    fi
done < "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt

###MISC. RESISTANCE###
printf "\tMisc. Resistance:\n" >> "$sampl_out"
lineNum=$(cat TEMP_miscR_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no T-antigens were typed
    firstLine=$(head -n1 TEMP_miscR_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Resistance\n" >> "$sampl_out"
    printf "No_Resistance" >> TEMP_table_results.txt
else
    count=0
    misc_target=()
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
	    #print the header line to the 'SAMPL' output file
            printf "\t\t$line\n" >> "$sampl_out"
        else
	    #If misc. resistance target is greater than the contamination threshold then add that 
	    #misc. resistance target to the output array 'misc_target'
            printf "\t\t$line\n" >> "$sampl_out"
	    justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
	    justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
	    justMatchType=$(echo "$line" | awk -F"\t" '{print $2}')
	    if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
		echo "Target $justTarget is a match"
		#if [[ ! "$justMatchType" =~ [identical|imperfect] ]]
		#then
		#    misc_target+=("$justTarget($justDepth)")
		#else
		    misc_target+=("$justTarget($justDepth|$justMatchType)")
		#fi
	    fi
        fi
    done < TEMP_miscR_Results.txt
    #if the output array 'misc_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [ ${#misc_target[@]} -eq 0 ];
    then
        printf "No_Resistance" >> TEMP_table_results.txt
    else
	printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' >> TEMP_table_results.txt
    fi
fi
printf "\t" >> TEMP_table_results.txt

###Surface / Secretory Protein Output (Not including T-Antigens)###
printf "\tProtein Targets:\n" >> "$sampl_out"
lineNum=$(cat TEMP_protein_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no surface/secretory targets were typed
    firstLine=$(head -n1 TEMP_protein_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Protein_Targets\n" >> "$sampl_out"
    printf "No_Protein_Targets" >> TEMP_table_results.txt
else
    count=0
    prot_target=()
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
	    #print the header line to the 'SAMPL' output file
            printf "\t\t$line\n" >> "$sampl_out"
        else
	    #If surface/secretory target is greater than the contamination threshold then add that
            #target to the output array 'prot_target'
            printf "\t\t$line\n" >> "$sampl_out"
            justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
	    justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
		echo "Target $justTarget is a match"
                #printf "$justTarget-($justDepth);" >> TEMP_table_results.txt
                prot_target+=("$justTarget($justDepth)")
            fi
        fi
    done < TEMP_no-Tantigen_Results.txt
    #done < TEMP_protein_Results.txt
    #if the output array 'prot_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [ ${#prot_target[@]} -eq 0 ];
    then
        printf "No_Protein_Targets" >> TEMP_table_results.txt
    else
	printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';' >> TEMP_table_results.txt
    fi
fi
printf "\n" >> "$sampl_out"
printf "\n" >> TEMP_table_results.txt

#cat "$sampl_out" >> "$batch_out"/SAMPL_GAS_"$batch_name"_Typing_Results.txt
#cat TEMP_table_results.txt >> "$batch_out"/TABLE_GAS_"$batch_name"_Typing_Results.txt


###Unload Modules###
module unload perl/5.12.3
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload Python/2.7
module unload samtools/0.1.18
module unload bowtie2/2.1.0
module unload freebayes/0.9.21
module unload prodigal/2.60



###PBP_ID OUTPUT###
#printf "\tPBP_ID Code:\n" >> "$sampl_out"
#count=0
#while read -r line
#do
#    count=$(( $count + 1 ))
#    justPBPs=$(echo "$line" | cut -f2-4)
#    if [[ "$count" -eq 1 ]]
#    then
#        printf "\t\t$justPBPs\n" >> "$sampl_out"
#        #printf "$justPBPs\t" >> TEMP_table_title.txt
#    else
#        printf "\t\t$justPBPs\n" >> "$sampl_out"
#        printf "$justPBPs\t" >> TEMP_table_results.txt
#    fi
#done < Final_pbpID_output.txt
