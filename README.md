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
===========

running DEMO (python > 3.5) - single control (-s)
-------------------------------------------------
```
> cd <coverageMaster_dir>\DEMO    
> gene=PGM1 && co=control.PGM1.cov && python ~/coverageMaster/coverageMaster.py test.PGM1.cov test.PGM1.report.txt $gene -s $co  -r ref.PGM1 -o test.PGM1
```

How to use it
=============

__COV file creation__

_BED = location of relevant genomic regions (e.g. REFSEQ or Probes or Gene Panels) in .bed format, one region per line_  

`> samtools depth -a -b <BED> <BAM> > <samplename>.cov`

__report.txt file creation__
  
`> samtools flagstat <BAM> > <samplename>.report.txt`
  
__Reference COV creation__

  1. copy or link your COV files in \<COV_folder\>

  2. create your reference per chromosome (scripts in TOOLS)

  `> python create_reference_cov.py <COV_folder>`

  3. concatenate all chromosomes

  `> for i in {1..22} X Y; do cat total*chr$i.res >> total_ref; done`

  4. calculate mean and std 

  `> python create_total_ref.py total_ref > total_ref_m_std`



Tips
----

*    To use more controls, put  all controls.cov and the related controls.report.txt in the same folder. Create a txt file with the absolute location of the .cov files (e.g. `ls -1 COV_folder > controls`) and use `-c controls` instead of -s  
*   to inspect more genes, create a txt file with genes separated by one space or one per line and give the filename as input (e.g `gene=<genelist> && ...` )  
*  to inspect a region just replace gene with chromosomal position chr:start-end (e.g `gene=chr1:123456-234567 && ...`). Zooming is not active for chromosomal position  
*   to inspect multiple regions, create a bed file with one chromosomal position per line and use -b option.  
