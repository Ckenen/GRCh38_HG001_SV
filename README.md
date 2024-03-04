# GRCh38_HG001_SV

Generation of the genotyped and phased structure variants for HG001.

In brief, the datasets (PacBio CCS: PRJNA540705) of the Whole-Genome-Sequencing of HG001 were mapped to GRCh38 by using Minimap2. Then, the alignments were phased according to the GIAB phased SNPs (https://github.com/Ckenen/GRCh38_HG001_SNP_Indel) by using whatshap (haplotag). Next, the SVs were called, genotyped, and phased by using sniffles2 (see `1_Snakemake.smk` and `jupyter.ipynb`).