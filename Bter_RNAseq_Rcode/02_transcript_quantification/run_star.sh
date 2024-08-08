#!/bin/sh

for name in ./data_dir/B*R1*.gz;
do

name_2="$(echo "$name" | cut -d '/' -f 3 | cut -d '.' -f 1)"
echo "$name_2"

STAR --runThreadN 4 \
--genomeDir ./index_dir/ \
--readFilesCommand gunzip -c \
--readFilesIn "$name" ./data_dir/"$name_2".R2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix ./output_dir/"$name_2"