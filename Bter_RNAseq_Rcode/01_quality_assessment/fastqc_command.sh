#!/bin/sh
########################################################################
##
## Author: Joe Colgan (joscolgan), Hongfei Xu (hongfeifly)
## Name: run_fastqc.sh
##
## Purpose:
## This script takes pairs of compressed fastq files and performs 
## a basic quality assessment analysis using fastqc. 
## The script produces a html file as output summarises the results
## of the assessment.
##
########################################################################

## Take inputs from the command line:
input_dir=./input_dir
# dirctory storing compressed fastq files

## Create a results output directory:
mkdir output_dir

## Run fastqc for all compressed files in input directory and output
## to newly created results directory:
fastqc -t 10 input_dir/*.gz -o output_dir
