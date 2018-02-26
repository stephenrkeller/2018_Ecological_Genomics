#!/bin/bash

########################################################################
# Script for calling variants using bcftools mpileup and call commands
# Feb 2018
# S.R. Keller
########################################################################


############# Pre-processing steps ##############
# Should be sorted, fixmated, and dups removed

# samtools view -@ 4 -bS file.sam >file.bam

# samtools fixmate file.bam file_fixmate.bam

# Sorting
# samtools sort -@ 8 file_fixmate.bam -o file_fixmate.sorted.bam 

# Using version 0.19 (bug with newer versions on this step)
# samtools_019 rmdup file_fixmate.sorted.bam file_fixmate.sorted.rmdup.bam

# samtools index file_fixmate.sorted.rmdup.bam

#################################################

# Change dir to location of your bam file. 
cd /data/otau/cleanreads

/data/programs/bcftools/bin/bcftools mpileup \
  -Ou -f /data/otau/reference/OTAU.fna merged_fixmate.sorted.rmdup.bam \
  -C 50 --min-MAQ 20 --min-BAQ 20 --threads 4 --skip-indels | \
  bcftools call -Ov -mv >OTAU_2018.vcf
  
  
