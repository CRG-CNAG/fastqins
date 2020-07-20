---

<p align="center">
  <img src=".logo/fastqins.png"/>
</p>

FASTQINS is a Python pipeline to map transponson insertions from Tn-seq data. 

## Requirements
Specific libraries are required by FASTQINS. We provide a [requirements](./requirements.txt) file to install everything at once. To do so, you will need first to have [pip](https://pip.pypa.io/en/stable/installing/) installed and then run:

```bash
pip3 --version                      # Check if installed
sudo apt-get install python3-pip    # if you need to install pip, you can check installation with the previous command
pip3 install -r requirements.txt

```


In addition, we have as dependencies standard tools commonly used in high-throughput sequencing analysis:

  [Fastuniq](https://sourceforge.net/projects/fastuniq/) <br /> 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)<br />
  [Samtools](http://www.htslib.org/)<br/>
  [Bedtools](https://bedtools.readthedocs.io/en/latest/)

Path to this dependencies can still be defined by the user editing the file [dependencies](./bin/dependencies.py).

## Installation & Help

Download this repository and run:

```bash
python3 setup.py install
```

You may require to call it using sudo. Once installed, `fastqins` should be available anywhere in your terminal.

In the case you need to install the package in a specific directory of your system, you can call the argument *--install-lib* followed by a directory path:

```bash
python3 setup.py install --install-lib /custom/path/
```

## Example

Requirements to run an experiment are: 

  -i [fastq files with transposon mapped, if no -i2 is passed, single-end mapping by default] <br />
  -t [IR transposon sequence, expected to be found contiguous genome sequence] <br />
  -g [genome sequence, fasta or genbank format]  <br />
  -o [output directory to locate the results]

As example, we included a pair of files that you can use to test the pipeline functioning as:

```bash
fastqins -i ./test/test_read2.fastq.gz -i2 ./test/test_read1.fastq.gz -t TACGGACTTTATC -g ./test/NC_000912.fna -o test -v -r 0
```

To see additional arguments:
```bash
fastqins --help
```

##  Output Information:

The following files are generated as default output:
- \*_fw.qins - read counts of insertions mapping to forward strand \[[example](./test/output_test/test_read2_fw.qins)\]
- \*_rv.qins - read counts of insertions mapping to reverse strand \[[example](./test/output_test/test_read2_rv.qins)\]
- \*.qins - read counts of insertions mapping to both strands \[[example](./test/output_test/test_read2.qins)\]
- \*.bam - file generated with the aligned reads
- \*.log - log file with general features of the process run \[[example](./test/output_test/test_read2.log)\]

## Contact

This project has been fully developed at [Centre for Genomic Regulation](http://www.crg.eu/) at the group of [Design of Biological Systems](http://www.crg.eu/en/luis_serrano).

If you experience any problem at any step involving the program, you can use the 'Issues' page of this repository or contact:

[Miravet-Verde, Samuel](mailto:samuel.miravet@crg.eu)         
[Lluch-Senar, Maria](mailto:maria.lluch@crg.eu)           
[Serrano, Luis](mailto:luis.serrano@crg.eu)

## License

FASTQINS is under a common GNU GENERAL PUBLIC LICENSE. Plese, check [LICENSE](./LICENSE) for further information.

###### [2020] - Centre de Regulació Genòmica (CRG) - All Rights Reserved*

