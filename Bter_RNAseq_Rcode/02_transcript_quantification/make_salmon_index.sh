#!/bin/sh
########################################################################
##
## Author: Joe Colgan (joscolgan), Hongfei Xu (hongfeifly)
## Name: make_star_index.sh
##
## Purpose:
## This script generates a genome index using salmon, a common short read
## aligner for RNA-seq datasets. This script takes one input argument
## from the command line to generate the genome index. The arguments is:
## A reference genome assembly in FASTA format.
## The script outputs a folder containing the genome index.
##
########################################################################

## Take inputs from the command line:
reference_dir=./reference_dir
# dirctory storing reference genome assembly
index_dir=./index_dir
# dirctory storing transcriptomic index file generated for running Salmon

## Create salmon indices:
salmon index 
-t $reference/transcriptome.fa
-i $reference 
--decoys $index_dir/decoys.txt 
-k 31
