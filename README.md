# as_analysis
This repository contains a set of scripts for processing RNA-seq data. The analysis pipeline starts from aligned .bam files, from which gene annotation files (.gtf) are generated. These GTF files are then parsed to construct and populate structured MySQL databases for downstream analyses.

The dataset includes multiple Arabidopsis thaliana genotypes:
Col-0, tfiis, upf1, tfiis-upf1, upf3, and tfiis-upf3,
each subjected to the following experimental treatments:
NT (untreated), 1h (1 hour heat stress), and 1d (1 day heat stress),
with three biological replicates per condition.
