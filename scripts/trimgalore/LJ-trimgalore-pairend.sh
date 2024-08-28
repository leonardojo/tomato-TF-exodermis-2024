#!/bin/bash

# This script perform trimming in pair-end sequenced reads using Trim-galore!
# author: Leonardo Jo

#Listing read1 (R1) and read2 (R2) files
ls FASTQ/*_R1_001.fastq.gz > FQnames1.txt;
ls FASTQ/*_R2_001.fastq.gz > FQnames2.txt;
FQ1=FQnames1.txt
FQ2=FQnames2.txt

#Storing how many files in FQ1 into NUMFQ
NUMFQ=$(wc -l $FQ1 | awk '{print $1}')

#Create a folder for the trimmed reads
mkdir -p TRIMMED_READS/;

#Loop to trim reads
i=1
while [ "$i" -le "$NUMFQ" ]; do
  FOR=$(awk 'NR == n' n=$i $FQ1 | awk '{print $1}')
  REV=$(awk 'NR == n' n=$i $FQ2 | awk '{print $1}')
  FORname=$(basename $FOR | cut -d '_' -f 1)
  REVname=$(basename $REV | cut -d '_' -f 1)
  echo $FOR
  echo $REV
  echo $FORname
  echo $REVname
  if [ "$FORname" = "$REVname" ]
  then
	command1="trim_galore $FOR $REV -a GATCGGAAGAGCACA --paired  -o TRIMMED_READS/ --basename $FORname"
	echo "$file TrimGalore: trimming reads"  >>  logfile_trimming.txt
	echo $command1  >>  logfile_trimming.txt
  eval $command1
  else
    echo "Something went wrong, R1 is not equal R2, skip this sample" >> logfile_trimming.txt
  fi
  i=$(($i + 1))
done
