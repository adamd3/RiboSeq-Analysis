#!/bin/bash

# Set parameters:

# Directories
datadir="/mnt/sdb1/RibosomeProfiling/RAWDATA"
databasedir="/mnt/sdb1/RibosomeProfiling/DATABASES"
stardbdir="/mnt/sdb1/RibosomeProfiling/STAR_DATABASES"
scriptsdir="/home/adinan/Annotated_scripts/"
plotsdir="/home/adinan/PlotScripts_AD/"
templatedir="/mnt/sdb1/RibosomeProfiling/TEMPLATES"

# Database names
databases1="rRNA:rRNA/rRNA"
databases2="mRNA:mRNA/mRNA ncRNA2:ncRNA_other/ncRNA_other gDNA:genome/genome"

# colours used for R plots:
dbdatacol1="rRNA:616 vRNA:552 mRNA:494 ncRNA_other:78"
#blue, red, light green, brown
dbdatacol2="gDNA:142"
#yellow

# Adaptor sequence to be trimmed from reads + min length of trimmed reads
adaptor="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCA"
trimlen=25
export adaptor
export trimlen

#-----------------------------------------------------------------------

# Place all gzipped Fastq files in the current working directory.
#    (files should have the extension .fq.gz)

# gunzip files in parallel
parallel 'gunzip {}' ::: *.fq.gz

#Check number of cycles (most abundant read length and number of reads of this
#  length in the first 100 reads in each file):
for library in $(awk '{print $1}' libraries.txt)
do
   head -400 $library.fq | \
    awk '{if (NR%4==2) print length($1)}' | sort -n | uniq -c | sort -nr | \
    head -1 | awk '{printf "%s: %s x %s nt\n","'"$library"'",$1,$2}'
done

#-------------------------------------------------------------------------------

# trim reads in parallel
parallel 'fastx_clipper -Q33 -l "$trimlen" -a "$adaptor" -c -n \
    -v -i {} > {.}.trimmed.fq 2>> {.}.log.txt' ::: *.fq

#The read names line in trimmed fastq files (e.g.
#  "@M00964:54:000000000-A3D68:1:1101:16462:1518 1:N:0:1") should contain two
#  fields. So if it is >2, just save the first 2 fields, and if it is only 1,
#  then add the dummy field "1:N:0:1" to each read name line.

for library in $(awk '{print $1}' libraries.txt)
do
  test=$(head -1 $library.trimmed.fq | awk '{print NF}') #check no. fields in 1st line
  if [ $test != 2 ]; then
    awk '{if (NR%4==1) {print $1,"1:N:0:1"} else {print $0}}' \
      $library.trimmed.fq | \
    awk '{if (NR%4==1) {print $1,$2} else {print $0}}' > $library.temp1
    mv $library.temp1 $library.trimmed.fq
  fi
done


#-----------------------------------------------------------------------

#Check that all bowtie databases are present for mapping

for host in $(awk '{print $3}' libraries.txt | sort | uniq)
do
  for dbline in $(echo $databases1 $databases2)
  do
    dbbowtie=$(echo $dbline | awk -F: '{print $2}')
    dbdir=$(ls $databasedir/$host/$dbbowtie.*ebwt | wc -l | awk '{print $1}')
    if [ -z $dbdir ]; then
      echo "Can't find bowtie database $databasedir/$host/$dbbowtie"
    fi
  done
done

for virus in $(awk '{print $2}' libraries.txt | sort | uniq)
do
  test=$(echo $virus | awk '{if ($1=="NA") {print 0} else {print 1}}')
  if [ -n $test ]; then
    dbdir=$(ls $databasedir/$virus/$virus.*ebwt | wc -l | awk '{print $1}')
    if [ -z $dbdir ]; then
      echo "Can't find bowtie database $databasedir/$virus/$virus"
    fi
  fi
done


#-----------------------------------------------------------------------


# map to ribosomal RNA

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "rRNA" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonrRNA.fq \
       $databasedir/$hostname/rRNA/rRNA \
       -q ./$library.trimmed.fq \
        > $library.rRNA.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done

#-----------------------------------------------------------------------

# map remaining reads to vRNA

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "vRNA" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonvRNA.fq \
       $databasedir/$virus/$virus \
       -q $library.nonrRNA.fq \
        > $library.vRNA.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done


#-----------------------------------------------------------------------

# map remaining reads to mRNA
for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "mRNA" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonmRNA.fq \
       $databasedir/$hostname/mRNA/mRNA \
       -q $library.nonvRNA.fq \
        > $library.mRNA.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done

#-----------------------------------------------------------------

# Map to other ncRNAs

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "ncRNA_other" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonncRNA.fq \
       $databasedir/$hostname/ncRNA_other/ncRNA_other \
       -q $library.nonmRNA.fq \
        > ./$library.ncRNA_other.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done

#-----------------------------------------------------------------------

# map to the host genome with STAR

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    STAR --runMode alignReads \
    --runThreadN 8  --outFileNamePrefix $library. --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 8 --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 2  \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outMultimapperOrder Random --outSAMmultNmax 1 --genomeLoad LoadAndKeep \
    --readFilesIn $library.nonncRNA.fq \
    --genomeDir $stardbdir/$hostname
done

