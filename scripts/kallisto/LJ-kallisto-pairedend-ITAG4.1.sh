#!/bin/bash

# This script perform kallisto for generating a pseudocount table
# author: Leonardo Jo

#Listing trimmed read1 (R1) and read2 (R2) files
ls TRIMMED_READS/*_val_1.fq.gz > TRIMnames1.txt;
ls TRIMMED_READS/*_val_2.fq.gz > TRIMnames2.txt;
TRIM1=TRIMnames1.txt
TRIM2=TRIMnames2.txt

#Storing how many files in FQ1 into NUMFQ
NUMFQ=$(wc -l $TRIM1 | awk '{print $1}')

#Create a folder for the pseudcount tables
mkdir -p KALLISTO_COUNTS/

#Loop to Run Kallisto to get the number of counts for each sample library
i=1
while [ "$i" -le "$NUMFQ" ]; do
  ## isolate samples names
  FOR=$(awk 'NR == n' n=$i $TRIM1 | awk '{print $1}')
  REV=$(awk 'NR == n' n=$i $TRIM2 | awk '{print $1}')
  FORname=$(basename $FOR | cut -d '_' -f 1)
  REVname=$(basename $REV | cut -d '_' -f 1)

  echo "R1: $FOR" >>  kallisto.logfile.txt
  echo "R2: $REV" >>  kallisto.logfile.txt
  echo "Sample name R1: $FORname" >>  kallisto.logfile.txt
  echo "Sample name R2: $REVname" >>  kallisto.logfile.txt

  if [ $FORname == $REVname ]
  then
  outdir1=KALLISTO_COUNTS/${FORname}
  mkdir -p ${outdir1}
  ## Report the run of kallisto into the logfile
	echo "$FORname Kallisto: pseudocounts to tomato transcriptome ITAG4.1"  >>   kallisto.logfile.txt
	#Pseudo-mapping with Kallisto to generate counts for RNAseq (from Pachter lab; eXpress/Cufflinks are disconstinued and not recommended/supported)
	command1="kallisto quant -t 20 -b 100 -i /home/uu_bio_pep/ljo/genomeindexes/kallistoindex/ITAG4.1_cDNA -o ${outdir1} $FOR $REV >> kallisto.logfile.txt 2>&1"
	eval $command1
  echo $command1 >>   kallisto.logfile.txt
	echo " ------------ " >>   kallisto.logfile.txt

  else
  echo "Something went wrong, R1 is not equal R2, skip this sample"  >> errorLOG.txt
  fi
  i=$(($i + 1))
done
