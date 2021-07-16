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

running DEMO (python > 3.5) - single control (-s)

> cd \<coverageMaster_dir\>

> gene=PGM1 && co=control.PGM1.cov && python ~/coverageMaster/coverageMaster.py test.PGM1.cov test.PGM1.report.txt $gene -s $co  -r ref.PGM1 -o test.PGM1


Normal Procedure
----------------

1) COV file creation

BED = location of relevant genomic regions (e.g. REFSEQ or Probes or Gene Panels) in .bed format, one region per line

> samtools depth -a -b \<BED\> \<BAM file\> \> \<samplename\>.cov

2) report.txt file creation (stats file)
  
> samtools flagstat \<BAM\> > \<samplename\>.report.txt
  
3) Reference COV creation

  copy or link your COV files in \<COV_folder\>

  create your reference per chromosome

  > python create_reference_cov.py <COV_folder> 

  concatenate all chromosomes

  > for i in {1..22} X Y; do cat total*chr$i.res >> total_ref ;done

  calculate mean and std 

  > python create_total_ref.py total_ref > total_ref_m_std

3) 


Tips
----

i)    To use more controls, put in the same folder all controls.cov and the related controls.report.txt. Create a txt file with the absolute location of the .cov and use the -c option instead of -s

ii)   to inspect more genes, create a txt file with genes separated by one space or one per line and give the filename as input in the gene position

iii)  to inspect a region just replace $gene with chromosomal position chr:start-end. Zooming is not active for chromosomal position

iv)   to inspect multiple region, create a bed file with one chromosomal position per line and use the option -b 
