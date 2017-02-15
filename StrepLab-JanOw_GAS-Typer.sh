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
module load cutadapt/1.8
module load srst2/0.1.7

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
out_nameARG=ARG_"$just_name"
out_nameRES=RES_"$just_name"
out_namePLAS=PLAS_"$just_name"

###Call MLST###
srst2 --samtools_args "\-A" --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_pyogenes.fasta" --mlst_definitions "$allDB_dir/spyogenes.txt" --min_coverage 99.999
###Check and extract new MLST alleles###
MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt "$out_nameMLST"__*.Streptococcus_pyogenes.sorted.bam "$allDB_dir/Streptococcus_pyogenes.fasta"

###Call emm Type###
module unload perl/5.22.1
module load perl/5.16.1-MT
emm_typer.pl -1 "$readPair_1" -2 "$readPair_1" -r "$allDB_dir" -n "$out_nameEMM"
module unload perl/5.16.1-MT
module load perl/5.22.1

###Call GAS Misc Resistance###
GAS_miscRes_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -m GAS_miscR_Gene-DB_Final.fasta -v GAS_vancR_Gene-DB_Final.fasta -n "$just_name"

###Type Surface and Secretory Proteins###
GAS_surface-secretory_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -p GAS_protein_Gene-DB_Final.fasta -n "$out_namePROT"

###Type ARG-ANNOT Resistance Genes###
srst2 --samtools_args '\\-A' --input_pe "$readPair_1" "$readPair_2" --output "out_nameARG" --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db "$allDB_dir/ARGannot_r1.fasta"

###Type ResFinder Resistance Gene###
srst2 --samtools_args '\\-A' --input_pe "$readPair_1" "$readPair_2" --output "out_nameRES" --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db "$allDB_dir/ResFinder.fasta"

###Type PlasmidFinder Resistance Gene###
srst2 --samtools_args '\\-A' --input_pe "$readPair_1" "$readPair_2" --output "out_namePLAS" --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db "$allDB_dir/PlasmidFinder.fasta"

#Add in perl script to find contamination threshold here 
contamination_level=10





###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
tabl_out="TABLE_Isolate_Typing_results.txt"
sampl_out="SAMPLE_Isolate__Typing_Results.txt"
extract_arr=(PARCGBS-1 GYRAGBS-1 23SWT-1 23SWT-3 RPLDGBS-1 RPLDGBS-2 RPLVGBS-1 RPLVGBS-2 RPOBgbs-1 RPOBgbs-2 RPOBgbs-3 RPOBgbs-4)

printf "$just_name\n" >> "$sampl_out"
printf "$just_name\t" >> "$tabl_out"
###EMM TYPE OUTPUT###
printf "\tEmm_Type:\n" >> "$sampl_out"
lineNum=$(cat "$out_nameEMM"__emm-Type__Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no emm type was found
    firstLine=$(head -n1 Serotype_results.txt)
    printf "\t\t$firstLine\n\t\tNo_emm_Type\n" >> "$sampl_out"
    printf "No_emm_Type\t" >> "$tabl_out"
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
            printf "$justTarget;" >> "$tabl_out"
        fi
    done < "$out_nameEMM"__emm-Type__Results.txt
fi
printf "\t" >> "$tabl_out"

###T-ANTIGEN TYPE OUTPUT###
printf "\tT-Antigen:\n" >> "$sampl_out"
lineNum=$(cat TEMP_protein_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no T-antigens were typed
    firstLine=$(head -n1 TEMP_protein_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_T-Antigen_Targets\n" >> "$sampl_out"
    printf "No_T-Antigen_Targets" >> "$tabl_out"
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
		if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
		then
		    echo "Target $justTarget is a match"
		    Tantigen_target+=("$justTarget")
		fi
	    else
		#if target is not a 'T_Antigen' then print it out to the non-T Antigen protein typing file
		printf "$line\n" >> TEMP_no-Tantigen_Results.txt
	    fi
        fi
    done < TEMP_protein_Results.txt
    #if the output array 'Tantigen_target' is not empty, print out the sorted types to the output files
    if [[ ${#Tantigen_target[@]} -eq 0 ]] 
    then
	printf "No_T-Antigen_Targets" >> "$tabl_out"
	printf "\t\tNo_T-Antigen_Targets\n" >> "$sampl_out"
    else
	printf '%s\n' "${Tantigen_target[@]}" | sort | tr '\n' ';'
	printf '%s\n' "${Tantigen_target[@]}" | sort | tr '\n' ';' >> "$tabl_out"
    fi
fi
printf "\t" >> "$tabl_out"

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
        printf "$MLST_tabl\t" >> "$tabl_out"
    fi
done < "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt

###MISC. RESISTANCE###
printf "\tMisc. Resistance:\n" >> "$sampl_out"
lineNum=$(cat TEMP_miscR_Results.txt | wc -l)
misc_contamination_level=7
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no T-antigens were typed
    firstLine=$(head -n1 TEMP_miscR_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Resistance\n" >> "$sampl_out"
    printf "No_Resistance" >> "$tabl_out"
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
	    if [[ $(echo "$justDepth > $misc_contamination_level" | bc) -eq 1 ]]
            then
                #this condition checks if the target is one of the extraction alleles. If so, it will append (extract) to output
                if [[ "$justMatchType" == "imperfect" && " ${extract_arr[@]} " =~ " ${justTarget} " ]]
                then
                    echo "Target $justTarget is a match but needs extraction"
                    #misc_target+=("$justTarget($justDepth|$justMatchType)")
                    #printf "$justTarget;" >> TEMP_table_results.txt
                    misc_target+=("$justTarget(extract)")
                else
                    echo "Target $justTarget is a match"
                    misc_target+=("$justTarget")
                fi
	    fi
        fi
    done < TEMP_miscR_Results.txt
    #if the output array 'misc_target' is not empty, print out the sorted types to the output files
    if [ ${#misc_target[@]} -eq 0 ];
    then
        printf "No_Resistance" >> "$tabl_out"
    else
	printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' >> "$tabl_out"
    fi
fi
printf "\t" >> "$tabl_out"

###Surface / Secretory Protein Output (Not including T-Antigens)###
printf "\tProtein Targets:\n" >> "$sampl_out"
lineNum=$(cat TEMP_no-Tantigen_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no surface/secretory targets were typed
    firstLine=$(head -n1 TEMP_no-Tantigen_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Protein_Targets\n" >> "$sampl_out"
    printf "No_Protein_Targets" >> "$tabl_out"
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
                #prot_target+=("$justTarget($justDepth)")
                prot_target+=("$justTarget")		
            fi
        fi
    done < TEMP_no-Tantigen_Results.txt
    #if the output array 'prot_target' is not empty, print out the sorted types to the output files
    if [ ${#prot_target[@]} -eq 0 ];
    then
        printf "No_Protein_Targets" >> "$tabl_out"
    else
	printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';' >> "$tabl_out"
    fi
fi
printf "\t" >> "$tabl_out"
#printf "\n" >> "$sampl_out"
#printf "\n" >> "$tabl_out"

###ARG-ANNOT and ResFinder Resistance Gene Typing Output###
printf "\tGeneral Resistance Targets:\n\t\tDB_Target\tMatch_Type\tDepth\n" >> "$sampl_out"
genRes_target=()
if [[ -s out_nameARG__fullgenes__ARGannot_r1__results.txt ]]
then
    count=0
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -ne 1 ]]
        then
            isIdentical="identical"
            justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
            if [[ -n "$justDiffs" ]]
            then
                isIdentical="imperfect"
            fi
            justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
            justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
            printf "\t\tARGannot_r1_$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
                echo "Target $justTarget is a match"
                genRes_target+=("ARGannot-r1_$justTarget")
            fi
        fi
    done < out_nameARG__fullgenes__ARGannot_r1__results.txt
fi
if [[ -s out_nameRES__fullgenes__ResFinder__results.txt ]]
then
    count=0
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -ne 1 ]]
        then
            isIdentical="identical"
            justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
            if [[ -n "$justDiffs" ]]
            then
                isIdentical="imperfect"
            fi
            justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
            justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
            printf "\t\tResFinder_$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
                echo "Target $justTarget is a match"
                genRes_target+=("ResFinder_$justTarget")
            fi
        fi
    done < out_nameRES__fullgenes__ResFinder__results.txt
fi
#if the output array 'genRes_target' is not empty, print out the sorted types to the output files
if [ ${#genRes_target[@]} -eq 0 ];
then
    printf "No_Gen_Resistance_Targets" >> "$tabl_out"
    printf "\t\tNo_Gen_Resistance_Targets\n" >> "$sampl_out"
else
    printf '%s\n' "${genRes_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
    printf '%s\n' "${genRes_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
fi
printf "\t" >> "$tabl_out"

###PlasmidFinder Plasmid Typing Output###
printf "\tPlasmid Prediction Targets:\n\t\tTarget\tMatch_Type\tDepth\n" >> "$sampl_out"
if [[ -s out_namePLAS__fullgenes__PlasmidFinder__results.txt ]]
then
    count=0
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -ne 1 ]]
        then
            isIdentical="identical"
            justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
            if [[ -n "$justDiffs" ]]
            then
                isIdentical="imperfect"
            fi
            justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
            justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
            printf "\t\t$justTarget\t$isIdentical\t$justDepth" >> "$sampl_out"
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
                echo "Target $justTarget is a match"
                plas_target+=("$justTarget")
            fi
        fi
    done < out_namePLAS__fullgenes__PlasmidFinder__results.txt
fi
#if the output array 'genRes_target' is not empty, print out the sorted types to the output files
if [ ${#plas_target[@]} -eq 0 ];
then
    printf "No_Plasmid_Targets" >> "$tabl_out"
    printf "\t\tNo_Plasmid_Targets" >> "$sampl_out"
else
    printf '%s\n' "${plas_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
    printf '%s\n' "${plas_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
fi
printf "\n\n" >> "$sampl_out"
printf "\n" >> "$tabl_out"




###Unload Modules###
module unload perl/5.22.1
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload freebayes/0.9.21
module unload prodigal/2.60
module unload cutadapt/1.8
module unload srst2/0.1.7
