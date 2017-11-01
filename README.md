# RiboSeq-Analysis
#!/bin/bash

#-----------------------------------------------------------------------

# Analysis of Ribosome Profiling and RNA sequence data for
# virus infections of host cells

# Outline:
#  - Trim reads.
#  - Map sequentially to rRNA; vRNA; mRNA; ncRNA; host genome.
#  - Plot the read mapping statistics.
#  - Calculate the length distribution of reads + plot
#  - Calculate framing of reads on host mRNA transcripts + plot
#  - Assess contamination of samples by ribonucleoproteins (RNPs)

#-----------------------------------------------------------------------

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
dbdatacol1="rRNA:616 vRNA:552 mRNA:494 ncRNA_other:78" #blue, red, light green, purple, brown
dbdatacol2="gDNA:142" #yellow

# Adaptor sequence to be trimmed from reads + min length of trimmed reads
adaptor="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCA"
trimlen=25
