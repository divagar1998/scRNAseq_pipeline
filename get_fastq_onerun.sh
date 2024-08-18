#! /bin/bash
set -euo pipefail

#get fastq file from SRA and zip
#path to csv file with SRA run numbers
#format of input csv file name_of_group,run1,run2,run3...
csv_file=$1
#directory to download files to
out_dir=$2

ml sratoolkit
ml fastqc

while IFS="," read -r col1 col2
do
    cd $out_dir
    #creates directory with name col1 if it doesn't exist
    mkdir -p $col1
    cd $col1
    runs=([1]=$col2)
    for i in ${runs[@]}
    do
        # extract the fastq file
        fasterq-dump --threads 32 $i
        #gunzip the fastq file
        gzip --verbose ./"${i}_1.fastq"
	    gzip --verbose ./"${i}_2.fastq"
        # quality control of fastq file
        fastqc -t 32 ./"${i}_1.fastq.gz"
	    fastqc -t 32 ./"${i}_2.fastq.gz"
    done
done < $1