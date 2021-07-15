coverageMaster User Guide
================

## Introduction
coverageMaster is a copy number variant (CNV) calling algorithm based on depth-of-coverage maps to detect CNVs of size>50pb at nucleotide level in whole exome sequencing and whole genome sequencing data. coverageMaster analyzes the sequencing coverage in a multidimensional Wavelet compressed (nucleotide-like) space and the CNVs are inferred with Hidden Markov Models at the nucleotide-like level in the multiple Wavelet spaces. Through « zooming » across dimensions, this approach enables the punctual analysis of regions with altered depth and the visual inspection of the outcome at nucleotide level to further reduce the false positive rate. coverageMaster provides the graphical representations of the predicted CNVs for all the genes of interest, and optionally, a wig formatted file compatible with UCSC Genome Browser for detailed coverage visualization of the target regions.

## Method Overview

## Input Requirements
coverageMaster requires input depth-of-coverage maps to be mapped and processed by external tools. The coverage files provided as input are generated using samtools depth. For each coverage file, a statistics file should be generated using samtools flagstat.

coverageMaster utilizes a reference file with the average coverage and standard deviation of 15-20 samples processed with the same technology. This reference file isgenerated as following:
* create the coverage and statistics file for each sample with samtools. These pairs of files should be in the same folder.
* run python create_average_cov.py path_to_input_folder path_to_output_folder
* for i in $(ls path_to_output_folder); do cat $i >> tmp_ref
* run python coverage_process_total_ref tmp_ref > path_to_final_output



## Output
