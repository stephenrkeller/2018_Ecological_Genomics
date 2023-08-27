#!/bin/bash

########################################################################
# Script for calling variants using reads2snps commands
# Feb 2018
# S.R. Keller
########################################################################


/reads2snps/    \
  -bam merged_fixmate.sorted.rmdup.bam  \
  -bamref /data/otau/reference/OTAU.fna \
  -out OTAU_2018_reads2snps \ 
  -min 5 \
  -th1 0.9 \
  -nbth 10 \