#!/bin/bash
set -u

sample_id=${1}
split=${2}
concordant=${3}

# Function to output an error if files do not exist
file_exists () {
    [[ ! -s $1 ]] && \
    echo "
ERROR - CIRCLE_FINDER - $sample_id
===================================
Error $1 does not exist or is empty.
Stopped circle_finder.
No circular DNA was identified.
===================================
" && exit
}

awk '{print $4}' ${split} | sort | uniq -c > ${sample_id}.split.id-freq.txt
#This file "${sample_id}.split.id-freq.txt" will be used for collecting split id that have frequency equal to 4.
awk '$1=="2" {print $2}' ${sample_id}.split.id-freq.txt > ${sample_id}.split.id-freq2.txt
awk '$1=="4" {print $2}' ${sample_id}.split.id-freq.txt > ${sample_id}.split.id-freq4.txt

awk '{print $4}' ${concordant} | sort | uniq -c > ${sample_id}.concordant.id-freq.txt
#The following command will chose (may not be always true) one concordant and 2 split read

awk '$1=="3" {print $2}' ${sample_id}.concordant.id-freq.txt > ${sample_id}.concordant.id-freq3.txt
awk '$1>3 {print $2}' ${sample_id}.concordant.id-freq.txt > ${sample_id}.concordant.id-freqGr3.txt

# check if output files exists and are not empty
file_exists ${sample_id}.concordant.id-freq3.txt
file_exists ${sample_id}.concordant.id-freqGr3.txt

grep -w -Ff ${sample_id}.split.id-freq2.txt ${split} > ${sample_id}.split_freq2.txt
grep -w -Ff ${sample_id}.split.id-freq4.txt ${split} > ${sample_id}.split_freq4.txt

# check if output files exist and are not empty
file_exists ${sample_id}.split_freq2.txt
file_exists ${sample_id}.split_freq4.txt

#Selecting concordant pairs that were 1) mapped uniquely and 2) mapped on more than one loci (file "freqGr3.txt")
grep -w -Ff ${sample_id}.concordant.id-freq3.txt ${concordant} > ${sample_id}.concordant_freq3.txt
grep -w -Ff ${sample_id}.concordant.id-freqGr3.txt ${concordant} > ${sample_id}.concordant_freqGr3.txt

# check if output files exist and are not empty
file_exists ${sample_id}.concordant_freq3.txt
file_exists ${sample_id}.concordant_freqGr3.txt

#Step 7: Putting split read with same id in one line
sed 'N;s/\n/\t/' ${sample_id}.split_freq2.txt > ${sample_id}.split_freq2.oneline.txt
sed 'N;s/\n/\t/' ${sample_id}.split_freq4.txt > ${sample_id}.split_freq4.oneline.txt

# check if output files exist and are not empty
file_exists ${sample_id}.split_freq2.oneline.txt
file_exists ${sample_id}.split_freq4.oneline.txt

#Step 8: Split reads map on same chromosome and map on same strand. Finally extracting id (split read same chromosome, split read same strand), collecting all the split reads that had quality >0
awk '$1==$10 && $7==$16 && $6>0 && $15>0 {print $4} ' ${sample_id}.split_freq2.oneline.txt > \
    ${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

# check if output files exist and are not empty
file_exists ${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

#Step 9: Based on unique id I am extracting one continuously mapped reads and their partner mapped as split read (3 lines for each id) 
grep -w -Ff ${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt ${sample_id}.concordant_freq3.txt > \
    ${sample_id}.concordant_freq3.2SPLIT-1M.txt

# check if output files exist and are not empty
file_exists ${sample_id}.concordant_freq3.2SPLIT-1M.txt

#Step 10: Sorting based on read-id followed by length of mapped reads.
awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", $8)} 1' ${sample_id}.concordant_freq3.2SPLIT-1M.txt | \
    awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", $8)} 1' | \
    awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", $8)} 1' | \
    awk 'BEGIN{FS=OFS=" "} {if (($9=="M" && $NF=="H") || \
        ($9=="M" && $NF=="S"))  {printf ("%s\tfirst\n",$0)} \
        else if (($9=="S" && $NF=="M") || ($9=="H" && $NF=="M")) {printf ("%s\tsecond\n",$0)} \
        else  {printf ("%s\tconfusing\n",$0)}}' | \
    awk 'BEGIN{FS=OFS="\t"} {gsub("\ ", "", $8)} 1' | \
    awk '{printf ("%s\t%d\n",$0,($3-$2)+1)}' | \
    sort -k4,4 -k10,10n | sed 'N;N;s/\n/\t/g' | \
    awk '{if ($5==$15) {print $0}  \
        else if (($5=="1" && $15=="2" && $25=="1") || ($5=="2" && $15=="1" && $25=="2")) \
            {printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20)} \
        else if (($5=="1" && $15=="2" && $25=="2") || ($5=="2" && $15=="1" && $25=="1")) \
        {printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n", $11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)} }' \
    > ${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt

# check if output files exist and are not empty
file_exists ${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt

#Step 11: Unique number of microDNA with number of split reads
awk '$1==$11 && $1==$21 && $7==$17'  ${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt | \
    awk '($7=="+" && $27=="-") || ($7=="-" && $27=="+")' | \
    awk '{if ($17=="+" && $19=="second" && $12<$2 && $22>=$12 && $23<=$3) {printf ("%s\t%d\t%d\n",$1,$12,$3)} \
        else if ($7=="+" && $9=="second" && $2<$12 && $22>=$2 && $23<=$13) {printf ("%s\t%d\t%d\n",$1,$2,$13)} \
        else if ($17=="-" && $19=="second" && $12<$2 && $22>=$12 && $23<=$3) {printf ("%s\t%d\t%d\n",$1,$12,$3)} \
        else if ($7=="-" && $9=="second" && $2<$12 && $22>=$2 && $23<=$13) {printf ("%s\t%d\t%d\n",$1,$2,$13)} }' | \
    sort | uniq -c | awk '{printf ("%s\t%d\t%d\t%d\n",$2,$3,$4,$1)}' > ${sample_id}.microDNA-JT.txt

file_exists ${sample_id}.microDNA-JT.txt