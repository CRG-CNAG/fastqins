#!/usr/bin/env python3

# 2020 - Centre de Regulacio Genomica (CRG) - All Rights Reserved


# --------------------------------------------------
# environment
# --------------------------------------------------
import sys, os
import argparse
from fastqins.fastqins import fastqins as fq

def run_fastqins():
    fq(tn_reads=options.tn_reads, paired_reads=options.paired_reads,
       tn_seq=options.tn_seq,
       genome=options.genome,
       output_folder=options.output_folder,
       separate_orientations=options.separate_orientations,
       rm_pcr=options.rm_pcr, barcode_size=options.barcode_size,
       mismatches=options.mismatches, extension=options.extension,
       threads=options.threads, align_qual=options.align_qual,
       ins_calling=options.ins_calling, zeroes=options.zeroes, keep_multiple=options.keep_multiple,
       rm_inter=options.rm_inter, verbose=options.verbose)

    print('Moving files...')
    basename = options.tn_reads.split('/')[-1].split('.')[0]
    intermediate_dir = '{}/{}_intermediate_files/'.format(options.output_folder, basename)
    cmd = 'mv {}*.qins {}*.bam {}*.log {}/'.format(intermediate_dir, intermediate_dir, intermediate_dir, options.output_folder)
    print(cmd)
    os.system(cmd)

    if options.rm_inter:
        cmd = 'rm -fR {}'.format(intermediate_dir)
        print(cmd)
        os.system(cmd)

    print('Finished! You will find the following files in your directory:\n\t*_fw.qins - read counts of insertions mapping to forward strand\n\t*_rv.qins - read counts of insertions mapping to reverse strand\    ↪\n\t*.qins - read counts of insertions mapping to both strands\n\t*.bam - file generated with the aligned reads. Allows visual revision and additional processing\n\t*.log - log file with general features of     ↪\the process run\n\n')


# --------------------------------------------------
# argument parser
# --------------------------------------------------

parser = argparse.ArgumentParser(description = "FastQins extracts insertion positions from a fastq file. Please cite ''")
parser.add_argument('-i', '--reads_with_tn',
                    dest="tn_reads",
                    required=True,
                    type=str,
                    help="Reads with the transposon IR.")
parser.add_argument('-i2', '--paired_reads',
                    dest="paired_reads",
                    default='0',
                    type=str,
                    help="Paired reads (if available) for --reads_with_tn")
parser.add_argument('-t', '--tn_seq',
                    dest="tn_seq",
                    default="TACGGACTTTATC",
                    type=str,
                    help="Inverted Repeat Transposon sequence to trim. Default= TACGGACTTTATC")
parser.add_argument('-g', '--genome',
                    dest="genome",
                    required=True,
                    type=str,
                    help="Genome of reference to map insertions. Default= /users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk")
parser.add_argument('-o', '--output_folder',
                    dest="output_folder",
                    default="./",
                    type=str,
                    help="Output directory to write the files.")
parser.add_argument('-rm_pcr', '--rm_pcr',
                    dest="rm_pcr",
                    default=1,
                    type=int,
                    help="Either to remove or not the pcr dup. Default=True. Not recommended in highly selected passages.")
parser.add_argument('-so', '--separate_orientations',
                    dest="separate_orientations",
                    default=1,
                    type=int,
                    help="Return insertions separated by strand, default=True")
parser.add_argument('-b', '--barcode_size',
                    dest="barcode_size",
                    default=0,
                    type=int,
                    help="Extract barcodes or not. [NOT TESTED IN EVERY BARCODE CONDITION]")
parser.add_argument('-m', '--mismatches',
                    dest="mismatches",
                    default=0,
                    type=int,
                    help="Accepted mismacthes during bowtie alignment: 0 or 1, default=0")
parser.add_argument('-e', '--extension',
                    dest="extension",
                    default='',
                    type=str,
                    help="Extension to add to the name of the file. Default='' (nothing).")
parser.add_argument('-r', '--rm_inter',
                    dest="rm_inter",
                    default=0,
                    type=int,
                    help="Remove intermediate files. By default intermediate files are removed.")
parser.add_argument('-p', '--threads',
                    dest="threads",
                    default=1,
                    type=int,
                    help="Number of threads to use in bowtie2 alignment")
parser.add_argument('-q', '--align_qual',
                    dest="align_qual",
                    default=10,
                    type=int,
                    help="Minimum alignment quality")
parser.add_argument('-ic', '--ins_calling',
                    dest="ins_calling",
                    default='bed',
                    type=str,
                    help="Mode to map insertions:\n\t-awk [problems with Ins/Del]\n\t-pyt [faster for small datasets]\n\t-bed [default, more efficient for general cases]")
parser.add_argument('-z', "--zeroes",
                    dest="zeroes",
                    action="store_true",
                    help="Output reports positions with no insertions as 0 reads. Default: only inserted positions")
parser.add_argument('-mult', "--keep_multiple",
                    dest="keep_multiple",
                    action="store_true",
                    help="Only recommended for PE sequencing. If selected, keeps reads mapping multiple times.")
parser.add_argument('-v', "--verbose",
                    dest="verbose",
                    action="store_true",
                    help="increase output verbosity")
options = parser.parse_args()

if __name__=='__main__':
    run_fastqins()

# 2020 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
