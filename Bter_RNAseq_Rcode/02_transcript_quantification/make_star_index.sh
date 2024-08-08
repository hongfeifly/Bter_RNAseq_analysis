#!/bin/sh

## Create STAR indices
STAR --runThreadN 10
--runMode genomeGenerate
--genomeDir ./index_dir \
--genomeFastaFiles ./genome_dir/Bombus_terrestris.Bter_1.0.dna.toplevel.fa \
--sjdbGTFfile ./annotation_dir/Bombus_terrestris.Bter_1.0.55.gtf \
--sjdbOverhang 149