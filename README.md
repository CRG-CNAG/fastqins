# FASTQINS
FASTQINS is a Python pipeline to map transponson insertions from Tn-seq data. 

## Requirements
  [Fastuniq](https://sourceforge.net/projects/fastuniq/) <br /> 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)<br />
  [Samtools](http://www.htslib.org/)<br/>
  [Bedtools](https://bedtools.readthedocs.io/en/latest/)

## Installation & Help

Download this repository and run:

`python3 setup.py install`

You may require to call it using sudo. Once installed, `fastqins` should be available anywhere in your terminal. 

## Example

Requirements to run an experiment are: 

  -i [fastq files with transposon mapped, if no -i2 is passed, single-end mapping by default] <br />
  -t [IR transposon sequence, expected to be found contiguous genome sequence] <br />
  -g [genome sequence, fasta or genbank format]  <br />
  -o [output directory to locate the results]

As example, we included a pair of files that you can use to test the pipeline functioning as:

`fastqins -i ./test/test_read2.fastq.gz -i2 ./test/test_read2.fastq.gz -t TACGGACTTTATC -g ./test/NC_000912.fna -o test -v -r 0`

To see additional arguments:
`fastqins --help`

### 2020 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
