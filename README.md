Context & Citation:
===
Localised tissue-specific gene expression and gene duplications are important sources of morph differences in a social bumblebee
===
This repository contains scripts related to the transcriptomic analysis of tissues (brains, fat bodies, ovaries, and reproductive tissues (spermatheca, median oviduct, vagina)) collected from worker and queen bumblebees. Findings of the analysis are reported in the following manuscript:

Hongfei Xu, ..., Thomas J. Colgan. Localised tissue-specific gene expression and gene duplications are important sources of morph differences in a social bumblebee.

The following directory contains scripts for:

* **The quality assessment of RNA-seq FASTQ sequences.**
    * Performed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
* **The alignment of filtered reads.**
    * Performed using [STAR](https://github.com/alexdobin/STAR).
* **Data exploration, including differential expression analysis.**
    * Performed using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
* **Tissue specificity analysis.**
    * Performed using [tispec](https://rdrr.io/github/roonysgalbi/tispec).
* **dN/dS estimation.**
    * Performed using [orthologr](https://github.com/drostlab/orthologr?tab=readme-ov-file).
* **Weighted gene co-expression network analysis.**
    * Performed using [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html).
* **Gene Ontology term enrichment analysis.**
    * Performed using [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html).

Scripts here, especially those for transcript quantification, differential expression, and Gene Ontology enrichment analysis used modifications of scripts published along with the following manuscripts:

* [Colgan et al. (2019), Molecular ecology](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15047) - [Github repository](https://github.com/wurmlab/Bter_neonicotinoid_exposure_experiment)
* [Zhuang et al. (2023), Nature Communications](https://www.nature.com/articles/s41467-023-41198-6) - [Github repository](https://github.com/Joscolgan/bombus_mated_worker_analysis?tab=readme-ov-file)
