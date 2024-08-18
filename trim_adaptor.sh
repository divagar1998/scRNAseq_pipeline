#! /bin/bash
set -euo pipefail

# trims adaptors from fastq files given a csv file with fastq file paths and adapter sequences
#format of input csv file name_of_group,run1,run2,run3,run4, adapter sequence
csv_file=$1
#directory to download files to
out_dir=$2

ml fastqc
ml cutadapt

while IFS="," read -r col1 col2 col3 col4 col5 col6 
do
    cd $out_dir
    #creates directory with name col1 if it doesn't exist
    mkdir -p $col1
    cd $col1
    runs=([1]=$col2 [2]=$col3 [3]=$col4 [4]=$col5)
    for i in ${runs[@]}
    do
        cutadapt -a $col6 -A $col6 -o ./"${i}_1_cutadapt.fastq.gz" -p ./"${i}_2_cutadapt.fastq.gz" ./"${i}_1.fastq.gz" ./"${i}_2.fastq.gz"
        
        # quality control of fastq file
        fastqc -t 4 ./"${i}_1_cutadapt.fastq.gz"
        fastqc -t 4 ./"${i}_2_cutadapt.fastq.gz"
    done
done < $1