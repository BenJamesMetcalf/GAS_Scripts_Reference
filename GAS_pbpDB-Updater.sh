#!/bin/bash -l

. /usr/share/Modules/init/bash
module load Python/2.7
module load ncbi-blast+/2.2.29

#Comment blah...#

while getopts :r:u:t: option
do
    case $option in
        r) refSeq_path=$OPTARG;;
        u) update_PBP=$OPTARG;;
        t) table_out=$OPTARG;;
    esac
done

###Check if the PBP reference paths exist and are correct###
if [[ -z "$refSeq_path" ]]
then
    echo "No PBP reference database path argument given."
    exit 1
else
    refSeq_p2=$(echo $refSeq_path | sed 's/\/$//g')
    refSeq="${refSeq_p2}"
    echo "refSeq: $refSeq"
fi

#Ref_1A=$(ls $refSeq/*1A*.faa)
Ref_2X=$(ls $refSeq/*2X*.faa)

#if [[ -e "$Ref_1A" && -e "$Ref_2X" ]]
if [[ -e "$Ref_2X" ]]
then
    echo "The PBP 2X reference database is: $Ref_2X"
else
    echo "The PBP 2X reference database is not in the correct format or doesn't exist."
    echo "Make sure it's in the $refSeq directory and has the correct naming format."
    echo "Specifically, the PBP sequence file must contain '2X' and have a '.faa' extension."
    exit 1
fi

if [[ -z "$update_PBP" ]]
then
    echo "A file containing the new PBP sequence info hasn't been given."
    exit 1
else
    if [[ -e "$update_PBP" ]]
    then
        update_PBP="${update_PBP}"
    else
        echo "The file containing the new PBP sequence info doesn't exist."
	exit 1
    fi
fi

if [[ -z "$table_out" ]]
then
    echo "The argument giving the filename of the TABLE formatted typing output hasn't been given."
    exit 1
else
    if [[ -e "$table_out" ]]
    then
        table_out="${table_out}"
    else
        echo "The argument giving the filename of the TABLE formatted typing output doesn't exist."
	exit 1
    fi
fi

eval refSeq=$refSeq
eval pbp2X=$pbp2X
eval update_PBP=$update_PBP
eval table_out=$table_out





###Subroutines###
blastTyper () {
    pbpName=$(echo "$2" | sed 's/.*bLactam_\(.*\)-.*\.faa/\1/g')
    SampleName=$(echo "$1" | awk -F"/" '{print $(NF-1)}')
    typedPath=$(dirname "$5")
    pbpAlleleID=""
    #echo "1: $1 | 2: $2 | 3: $3 | 4: $4 | 5: $5 | 6: $6"

    pbpStop="no"
    if [[ ! -s "$1" ]]
    then
        printf "The following extracted PBP sequence file is empty:\n"
        printf "$1\n"
	pbpAlleleID="ERROR"
        pbpStop="yes"
    fi

    #isBlastDB=$(ls "$refSeq"/"$pbpName"_prot_blast_db*)
    isBlastDB=$(ls "$refSeq"/Blast_bLactam_"$pbpName"_prot_DB*)
    if [[ "$pbpStop" == "no" ]]
    then
        if [[ -z "$isBlastDB" ]]
        then
            ### Make blast database
            makeblastdb -in "$2" -dbtype prot -out "$refSeq"/Blast_bLactam_"$pbpName"_prot_DB
            blastp -db "$refSeq"/Blast_bLactam_"$pbpName"_prot_DB -query "$1" -outfmt 6 -out "$SampleName"_blast-out_"$pbpName".txt
        else
            blastp -db "$refSeq"/Blast_bLactam_"$pbpName"_prot_DB -query "$1" -outfmt 6 -out "$SampleName"_blast-out_"$pbpName".txt
        fi

        cat "$SampleName"_blast-out_"$pbpName".txt | sort -r -n -k12,12 -k3,3 -k4,4 | sed -n 1p > prot_best_blast.txt
	printf "\n"
        cat prot_best_blast.txt
        geneLen=$(cat "$1" | grep -v '>' | tr '\n' ' ' | sed 's/ //g' | wc -c)
        blastIden=$(cat prot_best_blast.txt | awk -F"\t" '{print $3}')
        blastLen=$(cat prot_best_blast.txt | awk -F"\t" '{print $4}')
        blastName=$(cat prot_best_blast.txt | awk -F"\t" '{print $2}')

        echo "Best blast identity is: $blastIden"
        echo "Best blast length is: $blastLen"
        echo "PBP gene length is: $geneLen"

        if [[ $(echo "$blastIden == 100" | bc) -eq 1 && $(echo "$blastLen == $geneLen" | bc) -eq 1 ]]
        then
            pbpAlleleID=$(echo "$blastName" | sed 's/^\([0-9]\+\)||.*/\1/g')
	    echo "Found match in database. Sequence match is $pbpAlleleID."
        else
            echo "Didn't find match.  Will add sequence to database."
            ###Modify header (including adding number ID)###
            oldMaxID=$(cat "$2" | grep '>' | sed 's/>\([0-9]\+\)||.*/\1/g' | sort -n | tail -n1)
            pbpAlleleID=$(echo "$oldMaxID + 1" | bc)

            while read line
            do
                if [[ $line =~ ^\>.* ]]
                then
                    echo $line | sed "s/>.*/>$pbpAlleleID||GBS_$pbpName/g" >> Single_newRef_Seq.faa
                else
                    echo $line >> Single_newRef_Seq.faa
                fi
            done < "$1"

	    printf "New PBP Sequence\n"
	    cat Single_newRef_Seq.faa
            ###Output the new allele to Final_newRef_Seq.faa and add it to the database###
            cat Single_newRef_Seq.faa >> "$2"
            echo "Remaking blast database"
	    rm "$refSeq"/Blast_bLactam_"$pbpName"_prot_DB*
	    makeblastdb -in "$2" -dbtype prot -out "$refSeq"/Blast_bLactam_"$pbpName"_prot_DB
            rm Single_newRef_Seq.faa
        fi
    fi  

       ###Add in the code to update the TABLE pbp output file with the pbp IDs for the new alleles###
       while read line
       do
           lineName=$(echo "$line" | awk -F"\t" '{print $1}')
	   samplName=$(echo "$3" | sed 's/PBP_//g')
           if [[ "$samplName" == "$lineName" ]]
           then
	       pbp_out=$(echo "$line" | awk -F"\t" '{print $22}')
               #if [[ "$4" == "1A" ]]
               #then
	       #   pbp1A_out=$(echo $pbp_out | sed 's/NEW:/'$pbpAlleleID':/g')
               #    echo "$line" | awk -v var="$pbp1A_out" -v OFS='\t' '{$11=var; print }' >> "$5"_PRE
               #elif [[ "$4" == "2X" ]]
	       if [[ "$4" == "2X" ]]
               then
		   pbp2X_out=$(echo $pbp_out | sed 's/NEW/'$pbpAlleleID'/g')
		   echo "$line" | awk -v var="$pbp2X_out" -v OFS='\t' '{$22=var; print }' >> "$5"_PRE
               fi
           else
               echo "$line" >> "$5"_PRE
           fi
       done < "$5"
       rm "$5"
       mv "$5"_PRE "$5"

    ###Remove Files###
    rm prot_best_blast.txt
    rm "$SampleName"_blast-out_"$pbpName".txt
}




###Start Doing Stuff###
cp "$table_out" "$table_out"_OLD
while read line
do
    sampleID=$(echo "$line" | awk -F"\t" '{print $1}')
    pbpPath=$(echo "$line" | awk -F"\t" '{print $2}')
    pbpAllele=$(echo "$line" | awk -F"\t" '{print $3}')
    #echo "sample: $sampleID | path: $pbpPath | allele: $pbpAllele"
    #if [[ "$pbpAllele" == "1A" ]]
    #then
    #    blastTyper "$pbpPath" "$Ref_1A" "$sampleID" "$pbpAllele" "$table_out"
    #elif [[ "$pbpAllele" == "2X" ]]
    if [[ "$pbpAllele" == "2X" ]]
    then
        blastTyper "$pbpPath" "$Ref_2X" "$sampleID" "$pbpAllele" "$table_out"
    fi
done < "$update_PBP"


module unload Python/2.7
module unload ncbi-blast+/2.2.29
