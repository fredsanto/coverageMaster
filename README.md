coverageMaster User Guide
================

## Introduction
coverageMaster is a copy number variant (CNV) calling algorithm based on depth-of-coverage maps to detect CNVs of any size at nucleotide level in Whole Exome Sequencing and Whole Genome Sequencing data. CoverageMaster analyzes the sequencing coverage in a multidimensional Wavelet compressed (nucleotide-like) space and the CNVs are inferred with Hidden Markov Models at the nucleotide-like level. Through « zooming » across dimensions, the punctual analysis of regions with altered depth and the visual inspection of the outcome at nucleotide level to further reduce the false positive rate are enabled.   
coverageMaster provides the graphical representations of the predicted CNVs for all the genes of interest, and optionally, a wig formatted file compatible with UCSC Genome Browser for detailed coverage visualization of the target regions. 

## Dependencies
* Python in version 3.6.8

The following libraries need to be previously installed:  
* numpy (version 1.16.2)  
* pandas (version 1.1.5)
* sympy (version 1.0)  
* PyWavelets (version 3.5)  
* matplotlib (version 2.2.3)  
* scipy (version 1.2.1)
* more_itertools (version 8.8.0)

## Quick Start

__Running the DEMO (python > 3.5) with one control (-s)__

```
> cd <coverageMaster_dir>\DEMO    
> gene=PGM1 && co=control.PGM1.cov && python ~/coverageMaster/coverageMaster.py test.PGM1.cov test.PGM1.report.txt $gene -s $co -r ref.PGM1 -o test.PGM1
```

## Input Requirements 
coverageMaster requires input depth-of-coverage maps to be mapped and processed by external tools.

### COV files
The coverage files provided as input are generated using samtools depth. 

_BED = location of relevant genomic regions (e.g. REFSEQ or Probes or Gene Panels) in .bed format, one region per line_  

`> samtools depth -a -b <BED> <BAM> > <samplename>.cov`

### Statistics files
For each coverage file, a statistics file should be generated.

`> samtools flagstat <BAM> > <samplename>.report.txt`
  
### Reference COV file
coverageMaster utilizes a reference file with the average coverage and standard deviation of a 15-20 samples set. These samples should be processed with the same technology. The reference file is generated as following:

  1. copy or link your COV files in \<COV_folder\>

  2. create your reference per chromosome after creating your OUTP_folder (scripts in TOOLS)

  `> python create_reference_cov.py <COV_folder> <OUTP_folder>`

  3. concatenate all chromosomes

  `> cd <OUTP_folder && for i in {1..22} X Y; do cat total*chr$i.res >> total_ref; done`

  4. calculate mean and std 

  `> python create_total_ref.py total_ref > total_ref_m_std`

## Standard Usage 
`python coverageMaster.py <samplename>.cov <samplename>.report.txt GENE/GENELIST/REGION/REGIONLIST [-s <control.cov>]|[-c <control_list>]  -r total_ref_m_std -o <output_prefix>`

## Tips
*  To compare with more controls, put  all controls.cov and the related controls.report.txt in the same folder. Create a .txt file with the absolute location of the .cov files (e.g. `ls -1 COV_folder > controls`) and use `-c controls` instead of -s.  
*  To inspect more genes, create a .txt file `genelist` with gene names separated by one space and give the filename as input instead of the gene name.  
*  To inspect a chromosomal region, replace `gene` with chromosomal position chr:start-end (e.g `gene=chr1:123456-234567 && ...`). Zooming is not active for chromosomal position.  
*  To inspect multiple regions, create a bed file `positions.bed` with one chromosomal position per line and use `-b <positions>.bed`. Zooming is not active for this option. 
## More Tips
*  The standard windows size is +/-10 exons from the gene of interest. To extend (or reduce) this window, set `-e <number of exons outside the gene>`.
*  To get quickly the coverage plot of a given gene, add `-f` to the command line. It will produce the coverage plot even if there is no CNV call.
*  To add the mode of inheritance of all genes from the CGD database (v7.2021) to the coverage plot, use `-g <coverageMaster folder>/REF/CGD_07_21.inh`.
*  To speed up the process if you are interested in large CNVs only, set `-m 5`
## Output files
* __.CMcalls__
    * Position of the detected CNV as: chrom-gene-start-end and inheritance as identified in the Genomic Clinical Database (https://research.nhgri.nih.gov/CGD/).
* __.CMpositives.pdf__
    * Graphical representations of the predicted CNVs for all the genes of interest. In the first plot, the gene's structure in the exonic space is shown. In the second plot, the coverage profiles in the exonic space of the test sample, control and reference are represented, along with the HMM prediction and call.
* __.CMreport__
    * List of all the genes for which CNVs have been detected.
* __.CM.log__
    * .log file of the run.
## Docker Image
 To avoid the burden of installing all python library dependencies, a Docker image of coverageMaster can be created:
 ```
 > cd coverageMasterfolder (i.e. where coverageMaster.py is)
 > docker build -t coveragemaster .
 ```
 coveragemaster image should be now be listed when typing
 `> docker image ls`
 
 For example, to execute the demo:
 ```
 > cd DEMO
 > gene=PGM1 && co=control.PGM1.cov && docker run --rm -v `pwd`/:/data coveragemaster /data/test.PGM1.cov /data/test.PGM1.report.txt $gene -s /data/$co -r /data/ref.PGM1 -o /data/test.PGM1
 ```
