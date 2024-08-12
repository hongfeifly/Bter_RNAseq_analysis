#!/bin/sh
########################################################################
##
## Author: Joe Colgan (joscolgan), Hongfei Xu (hongfeifly)
## Name: make_star_index.sh
##
## Purpose:
## This script generates a genome index using STAR, a common short read
## aligner for RNA-seq datasets. This script takes three input arguments
## from the command line to generate the genome index. The arguments are:
## 1) A reference genome assembly in FASTA format.
## 2) A corresponding GFF for the reference genome assembly.  
## 3) The overhang, a value of the maximum read length minus 1.
## The script outputs a folder containing genome index.
##
########################################################################

## Create STAR indices
STAR --runThreadN 10
--runMode genomeGenerate
--genomeDir ./index_dir \
--genomeFastaFiles ./reference/genome.fa \
--sjdbGTFfile ./reference/genome.gtf \
--sjdbOverhang 149
