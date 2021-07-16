# CoverageMaster

CoverageMaster is a copy number variant (CNV) calling algorithm based on depth-of-coverage maps to detect CNVs of any size at nucleotide level in Whole Exome Sequencing and Whole Genome Sequencing data. CoverageMaster analyzes the sequencing coverage in a multidimensional Wavelet compressed (nucleotide-like) space and the CNVs are inferred with Hidden Markov Models at the nucleotide-like level in the multiple Wavelet spaces. Through « zooming » across dimensions, this approach enables the punctual analysis of regions with altered depth and the visual inspection of the outcome at nucleotide level to further reduce the false positive rate. 
CoverageMaster provides the graphical representations of the predicted CNVs for all the genes of interest, and optionally, a wig formatted file compatible with UCSC Genome Browser for detailed coverage visualization of the target regions. 

Dependencies
------------
The following libraries need to be previously installed:

numpy (>1.16.2)

sympy (>1.0)

PyWavelets (>3.5)

matplotlib (>2.2.3)

scipy (>1.2.1)


Quick Start
-----------

1) COV file creation

BED = location of relevant genomic regions (e.g. REFSEQ or Probes or Gene Panels) in .bed format, one region per line

> samtools depth -a -b \<BED\> \<BAM file\> \> \<samplename\>.cov

2) report.txt file creation (stats file)
  
> samtools flagstat \<BAM\> > \<samplename\>.report.txt
  
3) Reference COV creation

3) DEMO running

