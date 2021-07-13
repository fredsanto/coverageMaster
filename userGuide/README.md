coverageMaster User Guide
================

## Introduction
coverageMaster is a copy number variant (CNV) calling algorithm based on depth-of-coverage maps to detect CNVs of size>50pb at nucleotide level in Whole Exome Sequencing and Whole Genome Sequencing data. coverageMaster analyzes the sequencing coverage in a multidimensional Wavelet compressed (nucleotide-like) space and the CNVs are inferred with Hidden Markov Models at the nucleotide-like level in the multiple Wavelet spaces. Through « zooming » across dimensions, this approach enables the punctual analysis of regions with altered depth and the visual inspection of the outcome at nucleotide level to further reduce the false positive rate. coverageMaster provides the graphical representations of the predicted CNVs for all the genes of interest, and optionally, a wig formatted file compatible with UCSC Genome Browser for detailed coverage visualization of the target regions.

## Method Overview

## Input Requirements

## Output
