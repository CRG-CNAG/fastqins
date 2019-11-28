#!/usr/bin/env python3

####################################
# 2019 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
#
# Writen by Miravet-Verde, Samuel
# Last updated = 11/26/2019 (MM/DD/YYYY)
#
# Command line example:
# python /software/ls/fastqins/fastqins_web.py -i1 /users/lserrano/smiravet/ins_results/web_tests/test_read1.fastq.gz -i2 /users/lserrano/smiravet/ins_results/web_tests/test_read2.fastq.gz -t TACGGACTTTATC -g /users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk -o /users/lserrano/smiravet/ins_results/web_tests/your_test/ -p 1 -v -r 0
####################################


# --------------------------------------------------
# environment
# --------------------------------------------------
import re            # Only required by python mapper, can be removed without that part
import glob
import json
import os.path
import sys, os
import argparse
import subprocess
import numpy as np
from ruffus import *
from Bio import SeqIO
from collections import Counter

# --------------------------------------------------
# dependencies
# --------------------------------------------------
softdir  = '/software/ls/'
sys.path.insert(0, softdir+"fastqins")
fastuniq = softdir+'FastUniq/source/fastuniq'
bowtie2  = softdir+'bowtie2-2.2.9/bowtie2'
samtools = softdir+'/samtools-1.4/samtools'
bedtools = softdir+'bedtools2/bedtools'

fastuniq = 'fastuniq'
bowtie2  = 'bowtie2'
samtools = 'samtools'
bedtools = 'bedtools'

# --------------------------------------------------
# general functions
# --------------------------------------------------

def empty_transform(i, o):
    with open(o, 'w') as fo:
        fo.write(i)

def revcomp(seq):
    """
    Given a sequence string <seq>,
    returns the reverse complementary sequence
    """
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([comp.get(i, 'N') for i in seq])[::-1]

def return_original_read_length(i):
    """ Given an uncompressed fastq, returns the read length """
    cmd = "head -2 {} | tail -1".format(i)
    return len(subprocess.check_output(cmd, shell=True).strip())

def return_genome_lentgh(genome):
    if any([ext in genome for ext in ['gb', 'gbk', 'genbank']]):
        fmt = 'genbank'
    else:
        fmt = 'fasta'
    handle = open(genome, 'rU')
    for record in SeqIO.parse(handle, fmt):
        return len(record.seq)
    handle.close()

def process_cigar(cigar):
    if 'D' not in cigar and 'I' not in cigar:
        return 0
    else:
        dels = sum(np.array([int(i) for i in re.findall(r"(\d+)D", cigar)]))
        inss = sum(np.array([int(i) for i in re.findall(r"(\d+)I", cigar)]))
        return dels-inss

def load_qins(i):
    """ """
    rs = {}
    with open(i) as fi:
        for line in fi:
            line = line.strip().split()
            rs['{}\t{}'.format(line[0], line[1])]=int(line[2])
    return rs

# --------------------------------------------------
# pipeline functions
# --------------------------------------------------
def validate_genome(genome, outgenome):
    if genome!=outgenome:
        handle = open(genome, 'rU')
        for record in SeqIO.parse(handle, 'genbank'):
            handleout = open(outgenome, 'w')
            handleout.write(">{}\n{}".format(record.id, record.seq))
            handleout.close()
        handle.close()

def create_log_file(output_file, kwargs):
    """ Creates log file """
    pipelineDoc = """
#########################################
FastQins: Pipeline to analyze Tn-seq data.
#########################################

Python >=2.7 (>3 included) |
Author: Miravet-Verde, Samuel |
Last updated: 2019.11.21 |
Affiliation: Center for Genomic Regulation, Luis Serrano's lab |
email: samuelmiver@gmail.com; samuel.miravet@crg.eu |\n\n"""
    separator = '\n-----------------\n\n'
    with open(output_file, 'w') as logfile:
        logfile.write(pipelineDoc)
        logfile.write('Execution info\n-----------------\nversion:FASTQINS v1.0'+separator+'GENERAL INFORMATION:\n\nTn reads:{0}\nPaired reads:{1}\nTnSeq:{2}\nGenome:{3}\nOutDir:{4}\nPCR dup removal:{5}\nMismatches:{6}\nMapping threads:{7}\nMinimum quality alignment:{8}\nInsertion caller:{9}\nZeroes:{10}\nKeep reads multiple mapping:{11}'.format(kwargs['tn_reads'], kwargs['paired_reads'], kwargs['tn_seq'], kwargs['genome'], kwargs['output_folder'], kwargs['rm_pcr'], kwargs['mismatches'], kwargs['threads'], kwargs['align_qual'], kwargs['ins_calling'], kwargs['zeroes'], kwargs['keep_multiple'])+separator+'BOWTIE2 ALIGN. INFORMATION:\n\n')

def decompress(i, o):
    """ Gunzip input <i> to output <o> """
    cmd = "gunzip -c {} > {}".format(i, o)
    os.system(cmd)

def create_lis_file(inputs, output):
    """ Create the lis.ls file required to run FastUniq """
    with open(output, 'w') as fo:
        fo.write('{}\n{}'.format(inputs[0], inputs[1]))

def rm_pcr_duplicates(lisfile, outputs):
    """ Fastqins [https://sourceforge.net/projects/fastuniq/] parser to run fastuniq from a list of paired-seq list of files """
    cmd = '{} -i {} -t q -o {} -p {}'.format(fastuniq, lisfile, outputs[0], outputs[1])
    # -i lisfile > input files in list format from txt
    # -t q > output in two different files
    # -o outputs[0] > first output
    # -p outputs[1] > second output
    os.system(cmd)

def tn_trimming(reads, output, tn_seq):
    """ Trims the transposon from reads """
    cmd = "awk 'NR%4==2{s=match($0,/"+tn_seq+"/)+RLENGTH} NR%4~/[02]/{$0=substr($0,s)} 1' "+reads+" > "+output
    os.system(cmd)

def barcode_calling(reads, output, barcode_size):
    """
    Keeps record of the barcodes from reads,
    expected read structure: barcode-DNAg
    enerate a fasta with the identifiers and the associated barcode
    Trim the barcode in the read1 file (when line does not start with @IDE or +
    """
    ide = subprocess.run(['head', '-1', reads],stdout=subprocess.PIPE).stdout.decode('utf-8').split(':')[0]
    cmd =  "sed -n '1~4s/^"+ide+"/>"+ide[1:]+":/p;2~4p' "+reads+" | awk '{if($0~/^>/) {print $0} else {print substr($0,0,"+str(barcode_size)+")}}' > "+output
    os.system(cmd)

def genome_indexing(genome, output):
    """ Bowtie2 genome indexing """
    cmd = '{}-build --quiet -f {} {}'.format(bowtie2, genome, output)
    os.system(cmd)

def mapping(inputs, output, other_reads, kwargs):
    """ Bowtie2 mapping """
    # -p : threads
    # -x : genome
    # -1 : read1
    # -2 : read2
    # -S : sam output
    if other_reads=='0':
        cmd = '{} -p {} -x {} -N {} -U {} -S {} 2>> {}'.format(bowtie2, kwargs['threads'], inputs[1].replace('.1.bt2', ''), kwargs['mismatches'], inputs[0], output, inputs[2]) # SE mapping
    else:
        cmd = '{} -p {} -x {} -N {} -1 {} -2 {} -S {} 2>> {}'.format(bowtie2, kwargs['threads'], inputs[1].replace('.1.bt2', ''), kwargs['mismatches'], other_reads, inputs[0], output, inputs[2]) # PE mapping
    os.system(cmd)

def len_mult_selection(i, o, original_reads, keep_multiple):
    """ Filters out reads longer than they should and multiple mapped insertions if required """
    read_length = str(return_original_read_length(original_reads))
    if keep_multiple:
        keep_multiple = ' > {}'.format(o)
    else:
        keep_multiple = "| grep -v 'XS:i:' > {}".format(o)
    cmd = "awk 'length($10) <"+read_length+"' "+i+" | grep -v '^@'"+keep_multiple     # Read length filter and multiple mapped filter, transposon in the shorter ones
    os.system(cmd)

def recover_header(i, o):
    cmd = '{} view -H {} > {}'.format(samtools, i, o)
    os.system(cmd)

def sam_fixing2bam(i, o, align_qual):
    """ Reheads, sorts and transform to bam a sam filtered file """
    if align_qual:
        align_qual = '-q {}'.format(align_qual)
    cmd = 'cat {1} {2} | {0} view {3} -b | {0} sort > {4} '.format(samtools, i[0], i[1], align_qual, o)
    os.system(cmd)

def awk_insertion_calling(i, o, genome, orientation, zeroes):
    """ Insertion calling using awk. IT CANNOT PROCESS CIGARS """
    genome_size = str(return_genome_lentgh(genome))
    if orientation=='fw':
        cmd  = '{} view -F 0x10 {} | '.format(samtools, i)
        cmd += "awk '{print $4, $3}' "
    else:
        cmd  = '{} view -f 0x10 {} | '.format(samtools, i)
        cmd += "awk '{if($4+length($10)-1<="+genome_size+") {print $4+length($10)-1, $3} else if ($4+length($10)-1>"+genome_size+") {print $4+length($10)-1-"+genome_size+", $3}}' "
    cmd += " | sort | uniq -c | awk '{print $3, $2, $1}' > "+o
    os.system(cmd)

def bed_insertion_calling(i, o, genome, orientation, zeroes):
    """ Wrapper of bedtools to extract oriented insertions """
    if orientation=='fw':
        strand = '+'
    else:
        strand = '-'

    if zeroes:
        z = '-d'
        correction = 0
    else:
        z = '-dz'
        correction = 1
    cmd  = "{} genomecov {} -5 -strand {} -ibam {} -g {}".format(bedtools, '-dz', strand, i, genome)
    cmd += " | awk '{print $1, $2+"+str(correction)+", $3}' > "+o
    os.system(cmd)


def pyt_insertion_calling(i, o, genome, orientation, zeroes):
    """ Python insertion caller, also corrects CIGARS """

    intero = o.replace('.qins', '.sam')
    if orientation=='fw':
        cmd  = '{} view -F 0x10 {} > {} '.format(samtools, i, intero)
    else:
        cmd  = '{} view -f 0x10 {} > {} '.format(samtools, i, intero)
    os.system(cmd)

    genome_size = return_genome_lentgh(genome)
    positions = []
    with open(intero) as fi:
        for line in fi:
            if line[0]!='@':
                line = line.strip().split()
                gid  = line[2]         # genome id
                pos  = int(line[3])    # position
                l    = len(line[9])    # len of the read
                if orientation=='rv':
                    pos += l-1+process_cigar(line[5])
                    if pos>genome_size:
                        pos-=genome_size
                positions.append('{} {}'.format(gid, pos))
    with open(o, 'w') as fo:
        for k, v in Counter(positions).items():
            fo.write('{} {}\n'.format(k, v))

def count_insertions(i, o):
    """ Merge two orientations files """
    i1 = load_qins(i[0])
    i2 = load_qins(i[1])
    i1i2 = dict(Counter(i1)+Counter(i2))
    with open(o, 'w') as fo:
        for k, v in i1i2.items():
            fo.write('{} {}\n'.format(k, v))

# --------------------------------------------------
# pipeline control
# --------------------------------------------------

def fastqins(tn_reads, genome,
             paired_reads='0', output_folder='./', tn_seq='TACGGACTTTATC',
             separate_orientations=1, rm_pcr=1, barcode_size=0, mismatches=0, extension='',
             threads=1, align_qual=10,
             ins_calling='bed', zeroes=False, keep_multiple=False,
             rm_inter=0, verbose=False):
    """ Main function to run the fastqins pipeline. """

    tn_seq = tn_seq.upper()

    # Some useful variables
    if output_folder[-1]!='/':
        output_folder+='/'
    if mismatches>1:
        mismatches=1
    basename = tn_reads.split('/')[-1].split('.')[0]
    intermediate_dir = '{}{}_intermediate_files/'.format(output_folder, basename)
    kwargs           = {'tn_reads':tn_reads, 'paired_reads':paired_reads,
                        'tn_seq':tn_seq,
                        'genome':genome,
                        'output_folder':output_folder,
                        'rm_pcr':rm_pcr, 'mismatches':mismatches,
                        'threads':threads, 'align_qual':align_qual,
                        'ins_calling':ins_calling, 'zeroes':zeroes, 'keep_multiple': keep_multiple}

    ###
    # Initialize the pipeline
    fastqins_pipeline = Pipeline(name='fastqins')

    # make directory for all the files
    regexInputFiles = r'^(.+/)*(?P<SAMPLENAME>.+)\.fastq(\.gz)?$'                                    # Take as SAMPLENAME the string from the last '/' to the first '.' in filename
    mkd = fastqins_pipeline.mkdir(input    = tn_reads,
                                  filter   = formatter(regexInputFiles),
                                  output   = output_folder+"{SAMPLENAME[0]}_intermediate_files/")
    # Create log file:
    fastqins_pipeline.originate(task_func  = create_log_file,
                                output     = '{}{}.log'.format(intermediate_dir, basename),
                                extras     = [kwargs]).follows(mkd)


    # If genome in genbank format, create a fasta file
    genome_id = genome.split('/')[-1].split('.')[0]
    genome_ex = genome.split('/')[-1].split('.')[-1]
    # Determine the file type:
    if genome_ex in ['gb', 'gbk', 'genbank']:
        outgenome = '{}{}.fna'.format(intermediate_dir, genome_id)
    elif genome_ex in ['fa', 'fna', 'fasta', 'fast']:
        outgenome = genome
    else:
        print("GENOME FORMAT NOT SUPPORTED")
        raise SystemExit
    fastqins_pipeline.transform(task_func  = validate_genome,
                                input      = genome,
                                filter     = formatter(),
                                output     = outgenome).follows(mkd)

    # validate the genome and create fasta if required
    tn_reads_zip_check = tn_reads.endswith('.gz')
    tn_basename        = tn_reads.split('/')[-1].split('.')[0]
    fastqins_pipeline.transform(name       = 'decompress_tn_reads',
                                task_func  = decompress,
                                input      = tn_reads,
                                filter     = formatter(),
                                output     = intermediate_dir+tn_basename+'.fastq').active_if(tn_reads_zip_check).follows(mkd)

    paired = paired_reads!='0'
    paired_reads_zip_check = paired_reads.endswith('.gz')
    fastqins_pipeline.transform(name       = 'decompress_paired_reads',
                                task_func  = decompress,
                                input      = paired_reads,
                                filter     = formatter(),
                                output     = intermediate_dir+paired_reads.split('/')[-1].split('.')[0]+'.fastq').active_if(paired_reads_zip_check).follows(mkd)

    # Control input for removing pcr duplicates
    inputlis = []
    if tn_reads_zip_check:
        inputlis.append(fastqins_pipeline['decompress_tn_reads'])
    else:
        inputlis.append(tn_reads)
    if paired_reads_zip_check:
        inputlis.append(fastqins_pipeline['decompress_paired_reads'])
    else:
        inputlis.append(paired_reads)
    fastqins_pipeline.merge(task_func  = create_lis_file,
                            input      = inputlis,
                            output     = '{}lis.ls'.format(intermediate_dir)).active_if(paired, [rm_pcr])
    fastqins_pipeline.split(task_func  = rm_pcr_duplicates,
                            input      = create_lis_file,
                            output     = [intermediate_dir+tn_basename+'.filt',
                                          intermediate_dir+paired_reads.split('/')[-1].split('.')[0]+'.filt']).active_if(paired, [rm_pcr])

    # Tn calling
    if rm_pcr and paired:
        inputtnc = fastqins_pipeline['rm_pcr_duplicates'].parsed_args['output'][0]
    else:
        if tn_reads_zip_check:
            inputtnc = fastqins_pipeline['decompress_tn_reads']
        else:
            inputtnc = tn_reads
    fastqins_pipeline.transform(task_func  = tn_trimming,
                                input      = inputtnc,
                                filter     = formatter(),
                                extras     = [tn_seq],
                                output     = intermediate_dir+tn_basename+'.tn').follows(mkd).follows(fastqins_pipeline['rm_pcr_duplicates'])

    # Barcode calling
    inputbc = paired_reads
    if paired:
        if rm_pcr:
            inputbc = fastqins_pipeline['rm_pcr_duplicates'].parsed_args['output'][1]
        else:
            if tn_reads_zip_check:
                inputbc = fastqins_pipeline['decompress_paired_reads'].parsed_args['output']
    fastqins_pipeline.transform(task_func  = barcode_calling,
                                input      = inputbc,
                                filter     = formatter(),
                                output     = intermediate_dir+tn_basename+'_barcodes.fa',
                                extras     = [barcode_size]).follows(mkd).active_if(paired).active_if(int(barcode_size)>0)

    # Genome indexing
    fastqins_pipeline.subdivide(task_func  = genome_indexing,
                                input      = validate_genome,
                                filter     = formatter(),
                                output     = intermediate_dir+'{basename[0]}').active_if(os.path.isfile('{}{}.1.bt2'.format(intermediate_dir, genome_id))==False)

    # Mapping
    fastqins_pipeline.merge(task_func      = mapping,
                            input          = [fastqins_pipeline['tn_trimming'], '{}{}.1.bt2'.format(intermediate_dir, genome_id), create_log_file],
                            output         = intermediate_dir+tn_basename+'.sam',
                            extras         = [inputbc, kwargs]).follows(genome_indexing).follows(fastqins_pipeline['decompress_paired_reads'])

    # Sam filtering by quality and length of reads
    if tn_reads_zip_check:
        original_reads = fastqins_pipeline['decompress_tn_reads'].parsed_args['output']
    else:
        original_reads = tn_reads

    fastqins_pipeline.transform(task_func  = len_mult_selection,
                                input      = mapping,
                                filter     = suffix('.sam'),
                                output     = "_filt.sam",
                                extras     = [original_reads, keep_multiple])
    fastqins_pipeline.transform(task_func  = recover_header,
                                input      = mapping,
                                filter     = suffix('.sam'),
                                output     = "_head.sam")    # Required to properly work with bam files


    ###
    # Insertion calling
    # 3 different modes, run 2 times one per strand, differeneces between them, abstractly written
    outfmt = '.bam'
    if ins_calling=='awk':
        f = awk_insertion_calling
    elif ins_calling=='pyt':
        f = pyt_insertion_calling
    else:
        f = bed_insertion_calling    # Default --> bed
    fastqins_pipeline.merge(task_func      = sam_fixing2bam,
                            input          = [recover_header, len_mult_selection],
                            output         = intermediate_dir+tn_basename+outfmt,
                            extras         = [align_qual])
    for orientation in ['fw', 'rv']:
        fastqins_pipeline.transform(name       = '{}_{}_insertion_calling'.format(ins_calling, orientation),
                                    task_func  = f,
                                    input      = sam_fixing2bam,
                                    filter     = suffix(outfmt),
                                    output     = '_{}.qins'.format(orientation),
                                    extras     = [genome, orientation, zeroes])

    fastqins_pipeline.merge(task_func      = count_insertions,
                            input          = [fastqins_pipeline[ins_calling+'_fw_insertion_calling'], fastqins_pipeline[ins_calling+'_rv_insertion_calling']],
                            output         = intermediate_dir+tn_basename+'.qins') # Count insertions

    # Printout pipeline
    fastqins_pipeline.printout_graph(stream = 'flowchart.svg',
                                     pipeline_name   = 'fastqins | {} mode'.format(ins_calling),
                                     output_format   = 'svg',
                                     no_key_legend   = False,
                                     draw_vertically = True)

#    fastqins_pipeline.transform(task_func  = empty_transform,
#                                input      = fastqins_pipeline['merge_and_count'],
#                                filter     = suffix('.qins'),
#                                output     = ".qins2")


#    fastqins_pipeline.run(multiprocess = 4)
    fastqins_pipeline.run()


if __name__=='__main__':

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
                        default="/users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk",
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
                        help="Mode to map insertions:\n\t-awk [problems with Ins/Del]\n\t-pyt [faster for small datasets]\n\t-bed [default, more efficient for general cases")
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


    fastqins(tn_reads=options.tn_reads  , paired_reads=options.paired_reads,
             tn_seq=options.tn_seq,
             genome=options.genome,
             output_folder=options.output_folder,
             separate_orientations=options.separate_orientations,
             rm_pcr=options.rm_pcr, barcode_size=options.barcode_size,
             mismatches=options.mismatches, extension=options.extension,
             threads=options.threads, align_qual=options.align_qual,
             ins_calling=options.ins_calling, zeroes=options.zeroes, keep_multiple=options.keep_multiple,
             rm_inter=options.rm_inter, verbose=options.verbose)

####
# FUNCTIONS
####


sys.exit()

def project_identifier(fil, outputfolder):
    """
    Required to discriminate between @
    in the identifier and in the quality string
    """
    out = outputfolder+'id.txt'
    cmd='head -1 '+fil+' > '+out
    os.system(cmd)
    with open(out, 'rU') as fi:
        for line in fi:
            return line.split(':')[0]

def extract_insertions(out, read_length, genome_size, separate_orientations):
    """ core function to call insertions """

    # Filter by read length: only take reads that are shorter than usual: they have the transposon trimmed!
    # The position is equal to the first base mapped before the transposon starts
    # 'XS:' --> remove reads with multiple mapping

    # Forward strand
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_forward_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{print $4, $3}' > "+out+"_forward.ins"
    os.system(cmd)
    # Reverse strand
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_reverse_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{if($4+length($10)-1<="+genome_size+") {print $4+length($10)-1, $3} else if ($4+length($10)-1>"+genome_size+") {print $4+length($10)-1-"+genome_size+", $3}}' > "+out+"_reverse.ins"
    os.system(cmd)

    # Merge and extract info:
    cmd = "cat "+out+"_reverse.ins "+out+"_forward.ins | sort | uniq -c | awk '{print $3, $2, $1}' > "+out+".qins"
    os.system(cmd)

    if separate_orientations:
        # Merge and extract info:
        cmd = "cat "+out+"_reverse.ins | sort | uniq -c | awk '{print $3, $2, $1}' > "+out+"_reverse.qins"
        os.system(cmd)
        cmd = "cat "+out+"_forward.ins | sort | uniq -c | awk '{print $3, $2, $1}' > "+out+"_forward.qins"
        os.system(cmd)


def extract_identifier_position(out, read_length, genome_size):
    # Extracting barcode information
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_forward_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{print $1, $4}' > "+out+"_forward.ides"
    os.system(cmd)
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_reverse_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{if($4+length($10)-1<="+genome_size+") {print $1, $4+length($10)-1} else if ($4+length($10)-1>"+genome_size+") {print $1, $4+length($10)-1-"+genome_size+"}}' > "+out+"_reverse.ides"
    os.system(cmd)
    cmd = "cat "+out+"_reverse.ides "+out+"_forward.ides > "+out+".ides"
    os.system(cmd)

def finish_log_file(fastq_dir, fastq_filt, tn_seq, out, genome_size):

    # Compute stats:
    cmd = "awk '{s++}END{print s/4}' "+fastq_dir
    nr_reads = float(subprocess.check_output(cmd, shell=True).strip())
    if fastq_filt:
        cmd = "awk '{s++}END{print s/4}' "+fastq_filt
        nr_remain = float(subprocess.check_output(cmd, shell=True).strip())
        filt_line = str(int(nr_remain))+" ("+str(round(100.0*nr_remain/nr_reads, 2))+"%) passed the PCR duplicates filter; of these:"

        cmd = "grep "+tn_seq+" "+fastq_filt+" | wc -l"
        nr_IR = float(subprocess.check_output(cmd, shell=True).strip())
    else:
        nr_remain = nr_reads
        filt_line = "PCR duplicate removal not selected. "+str(int(nr_reads))+" considered; of these:"

        cmd = "grep "+tn_seq+" "+fastq_dir+" | wc -l"
        nr_IR = float(subprocess.check_output(cmd, shell=True).strip())

    cmd = "cat "+out+"_reverse.ins "+out+"_forward.ins | wc -l"
    nr_Tn = float(subprocess.check_output(cmd, shell=True).strip())
    cmd = "cat "+out+".qins | wc -l"
    nr_pos = int(subprocess.check_output(cmd, shell=True).strip())
    if nr_pos!=0:
        cmd = "awk '{sum+=$3} END{print sum}' "+out+".qins"
        nr_post_reads = int(subprocess.check_output(cmd, shell=True).strip())
    else:
        nr_post_reads = 0
    # Write log
    separator = '\n-----------------\n\n'
    with open(out+'.log', 'a') as logfile:
        logfile.write(separator+'FASTQINS INFORMATION:\n\n')
        text = "{0} reads provided; of these:\n  ".format(int(nr_reads))+filt_line
        text+= "\n    {0} ({1}%) presented the IR:{2} sequence; of these:\n      ".format(int(nr_IR), round(100.0*nr_IR/nr_remain, 2), tn_seq)
        text+= "{0} having IR ({1}%) were mapped unambiguously; in total:\n        ".format(int(nr_Tn), round(100.0*nr_Tn/nr_IR, 2))
        text+= "{0} insertion positions can be extracted from {1} reads (this number should be the same than previous line)\n".format(nr_pos, nr_post_reads)
        text+= "{0}% of the fastq was informative.".format(round(100.0*nr_Tn/nr_reads, 2))
        text+= "{0}% coverage (number of insertions per genome base).".format(round(100*int(nr_pos)/int(genome_size)))
        logfile.write(text)


def fastqins(read1        , read2   , Tn_seq      , genome    ,
             output_folder, separate_orientations , pcr_dup   , barcode_size, mismatches,
             extension    , rm_inter, verbose    ):

    check_everything([read1, read2])

    Tn_seq = Tn_seq.upper()

    if verbose:
        print('read1:', read1)
        print('read2:', read2)
        print('TnSeq:', Tn_seq)
        rcTn_seq = revcomp(Tn_seq)
        print('RC-TnSeq:', rcTn_seq)
        print('genome:', genome)
        print('outFolder:', output_folder)

    # Create the directory and the environment
    if output_folder[-1]!='/':
        output_folder+='/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_folder += extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # If genome in genbank format, create a fasta file
    genome, genome_size = testfasta_and_length(genome, output_folder)

    # Create log file:
    pipeline_name = 'pipeline_fastqins'
    pipelineDoc = pipeline_name + """
    Pipeline to analyze Tn-seq data.
    Python >=2.7 script. |
    Author: Miravet-Verde, Samuel |
    Last updated: 2018.12.10 |
    Affiliation: Center for Genomic Regulation, Luis Serrano's lab |
    email: samuelmiver@gmail.com; samuel.miravet@crg.eu |\n
    """
    separator = '\n-----------------\n\n'
    with open(output_folder+read1.split('/')[-1].split('.')[0]+'.log', 'w') as logfile:
        logfile.write(pipelineDoc)
        logfile.write('Execution info\n-----------------\nversion:FASTQINS v1.0'+separator+'GENERAL INFORMATION:\n\nR1:{0}\nR2:{1}\nTnSeq:{2}\nGenome:{3}\nGenome length:{4}\nOutDir:{5}\nPCR dup removal:{6}\nMismatches:{7}'.format(read1, read2, Tn_seq, genome, genome_size, output_folder, str(pcr_dup), str(mismatches))+separator+'BOWTIE2 ALIGN. INFORMATION:\n\n')

    midfile1 = output_folder+extension+read1.split('/')[-1]
    midfile2 = output_folder+extension+read2.split('/')[-1]

    # Decompress the file if compressed
    update_progress(options.analysisId, 'Uncompressing files', 4)
    marker = False
    if read1.endswith('gz'):
        cmd = 'gunzip -c '+read1+' > '+midfile1.replace('.gz', '')
        if verbose:
            print('decompressing read1 ... ')
        _ = os.system(cmd)
        midfile1 = midfile1.replace('.gz', '')
        marker = True

    if read2.endswith('gz'):
        cmd = 'gunzip -c '+read2+' > '+midfile2.replace('.gz', '')
        if verbose:
            print('decompressing read2 ... ')
        _ = os.system(cmd)
        midfile2 = midfile2.replace('.gz', '')
        read2_backup = str(midfile2)
    else:
        read2_backup = str(read2)

    # Get the read length
    cmd = 'head -2 '+read2_backup+' | tail -1'
    read_length = str(len(subprocess.check_output(cmd, shell=True))-1)

    # Remove duplicates
    update_progress(options.analysisId, 'Filtering PCR duplicates', 5)
    read2_filt = False
    if pcr_dup:
        fo = open(output_folder+'lis.ls', 'w')
        if marker:
            fo.write(midfile1+'\n'+midfile2)
        else:
            fo.write(read1+'\n'+read2)
        fo.close()
        if verbose:
            print('Removing PCR duplicates...')
        cmd = fastuniq+' -i '+output_folder+'lis.ls -t q -o '+midfile1+'.filt -p '+midfile2+'.filt'
        _ = os.system(cmd)
        if verbose:
            print('PCR Duplicates removed.')
        midfile1 += '.filt'
        midfile2 += '.filt'
        read2_filt = str(midfile2)

    # Transposon calling and trimming
    update_progress(options.analysisId, 'Transposon calling and trimming', 5)
    if marker or pcr_dup:
        cmd = "awk 'NR%4==2{s=match($0,/"+Tn_seq+"/)+RLENGTH} NR%4~/[02]/{$0=substr($0,s)} 1' "+midfile2+" > "+midfile2+".tn"
    else:
        cmd = "awk 'NR%4==2{s=match($0,/"+Tn_seq+"/)+RLENGTH} NR%4~/[02]/{$0=substr($0,s)} 1' "+read2+" > "+midfile2+".tn"
    if verbose:
        print('Calling the transposon and trimming it ... ')
    os.system(cmd)
    midfile2 += '.tn'

    # Barcode processement if present
    if barcode_size:
        if verbose:
            print('Calling the barcode with size '+str(barcode_size)+' and trimming it ...')
        # Generate a fasta with the identifiers and the associated barcode
        # Trim the barcode in the read1 file (when line does not start with @IDE or +
        if marker or pcr_dup:
            ide = project_identifier(midfile1, output_folder)
            cmd1 = "sed -n '1~4s/^"+str(ide)+"/>/p;2~4p' "+midfile1+" | awk '{if($0~/^>/) {print $0} else {print substr($0,0,"+str(barcode_size)+")}}' > "+midfile1+"_barcodes.fa"
            cmd2 = "awk '{if($0 ~/^"+str(ide)+"/ || $0 ~/+/) {print $0} else {print substr($0,"+str(int(barcode_size)+1)+")}}' "+midfile1+" > "+midfile1+".bc"
        else:
            ide = project_identifier(read1, output_folder)
            cmd1 = "sed -n '1~4s/^"+str(ide)+"/>/p;2~4p' "+read1+" | awk '{if($0~/^>/) {print $0} else {print substr($0,0,"+str(barcode_size)+")}}' > "+midfile1+"_barcodes.fa"
            cmd2 = "awk '{if($0 ~/^"+str(ide)+"/ || $0 ~/+/) {print $0} else {print substr($0,"+str(int(barcode_size)+1)+")}}' "+read1+" > "+midfile1+".bc"
        os.system(cmd1)
        os.system(cmd2)
        midfile1 += '.bc'

    # Index the genome
    update_progress(options.analysisId, 'Building genome index', 6)
    if verbose:
        print('Indexing the genome '+genome+' ...')
    cmd = bowtie2+'-build -f '+genome+' '+output_folder+'genome_bowtie'
    os.system(cmd)

    # Map the reads
    update_progress(options.analysisId, 'Mapping reads', 7)
    out = output_folder+midfile1.split('/')[-1].split('.')[0]
    # -p : threads
    # -x : genome
    # -1 : read1
    # -2 : read2
    # -S : sam output
    if marker or pcr_dup:
        cmd = bowtie2+' -p 12 -x '+output_folder+'genome_bowtie -N '+str(mismatches)+' -1 '+midfile1+' -2 '+midfile2+' -S '+out+'.sam 2>>'+out+'.log'
    else:
        cmd = bowtie2+' -p 12 -x '+output_folder+'genome_bowtie -N '+str(mismatches)+' -1 '+read1+' -2 '+midfile2+' -S '+out+'.sam 2>>'+out+'.log'
    if verbose:
        print('Mapping ...')
    os.system(cmd)

    # Transform to bam
    # -b : bam output
    # -S : sam input
    # -o : output name
    cmd = samtools+' view -b -S -o '+out+'.bam '+out+'.sam'
    if verbose:
        print('Creating bam ...')
    os.system(cmd)

    # Filter paired end
    # -F 0x04 : filtering unmapped reads
    # -f 0x02 : required to map PE
    # -q 30 : minimum alignment quality
    # cmd = samtools+' view -F 0x04 -f 0x02 -q 30 -b '+out+'.bam > '+out+'.paired_mapped'
    update_progress(options.analysisId, 'Filtering alignments', 8)
    cmd = samtools+' view -q 10 -F 0x10 '+out+'.bam > '+out+'_forward_filt_paired.sam'   # Analog to -F 16, 0x10 is 16 in decimal value
    os.system(cmd)
    cmd = samtools+' view -q 10 -f 0x10 '+out+'.bam > '+out+'_reverse_filt_paired.sam'
    os.system(cmd)
    if verbose:
        print('Filtering unmmaped and separating forward/reverse...')

    # Take unpaired but mapped reads
    # -F 0x04 : filtering unmapped reads
    # -F 0x02 : filtering PE reads
    # cmd = samtools+' view -F 0x04 -F 0x02 -b '+out+'.bam > '+out+'.unpaired_mapped'
    # print 'Extracting mapped that are non-paired...'
    # os.system(cmd)
    # Back to sam the 2 files
    # cmd = samtools+' view -h -o '+out+'_filt_paired.sam '+out+'.paired_mapped'
    # print 'Back to sam the paired mapped reads ...'
    # os.system(cmd)
    # cmd = samtools+' view -h -o '+out+'_filt_single.sam '+out+'.unpaired_mapped'
    # print 'Back to sam the singletons mapped ...'
    # os.system(cmd)

    # Extract insertions and relation with identifiers
    if verbose:
        print('Selecting reads with the transposon from paired reads N<'+str(read_length)+'...')
    update_progress(options.analysisId, 'Extracting transposition events', 9)
    extract_insertions(out, read_length, genome_size, separate_orientations=separate_orientations)
    if barcode_size:
        if verbose:
            print('Extracting barcode information...')
        extract_identifier_position(out, read_length, genome_size)
        cmd = 'cp '+output_folder+'*.ides '+output_folder.replace(extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/', '')
        os.system(cmd)

    # Finish the log file
    update_progress(options.analysisId, 'Writing files', 10)
    finish_log_file(read2_backup, read2_filt, Tn_seq, out, genome_size)
    # Put important files in its place
    cmd = 'cp '+output_folder+'*.qins '+output_folder.replace(extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/', '')
    os.system(cmd)
    cmd = 'cp '+output_folder+'*.log '+output_folder.replace(extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/', '')
    os.system(cmd)

    if rm_inter and '_intermediate_files' in output_folder:
         cmd = 'rm '+output_folder+'*'
         os.system(cmd)
         cmd = 'rmdir '+output_folder
         os.system(cmd)

         if verbose:
             print('Intermediate files removed.\n')

    # Say goodbye
    if verbose:
        print('Insertions extracted to '+out+".qins file\n\nEnjoy your day and remember:\n\nSTAND AND BE TRUE :D!\n")


# FOR WEB SERVER
def communicate_with_server(analysisId, status, hint='', port=50001):

    URL = "http://dbspipes.crg.es/api/v1"
    HEADERS = {'Content-type': 'application/json'}

    # sending get request and saving the response as response object
    URL = URL + '/analyses/' + str(analysisId)

    PARAMS = json.dumps({'status': status, 'hint': hint})

    try:
        print("Send server request: " + URL + " status: " + str(status) + " hint: " + hint)
        r = requests.put(url = URL, data=PARAMS, headers=HEADERS)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:  # This is the correct syntax
        print(e)
        sys.exit(1)

def update_progress(analysisId, job, n):

    URL = "http://dbspipes.crg.es/api/v1"
    HEADERS = {'Content-type': 'application/json'}

    # sending get request and saving the response as response object
    URL = URL + '/analyses/' + str(analysisId) + '/progress'

    PARAMS = json.dumps({"progress": {'job_name': job, 'task_n': n}})

    try:
        print("Send server request: " + URL)
        r = requests.post(url = URL, data=PARAMS, headers=HEADERS)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:  # This is the correct syntax
        print(e)
        sys.exit(1)

def write_jobid_files(OUTFILES, analysisId, bc):
    # Check that final output files exist
    # task_path is the last task's output directory
    if bc!=0:
        n = 6
        indexes, extensions = range(n), ['.qins']*3+['.log', '.ides']
    else:
        n = 5
        indexes, extensions = range(n), ['.qins']*3+['.log']

    print(OUTFILES)
    if len(OUTFILES)!=n:
        status = -1
        hint = 'No output files generated'
    else:
        for index, extension in zip(indexes, extensions):
            # We just check that the coverage qins file exists and is not empty
            if OUTFILES[index].endswith(extension) and os.path.getsize(OUTFILES[index]) > 0:
                status = 2
                hint=''
            else:
                status = -1
                hint='Empty output files'
    communicate_with_server(analysisId, status, hint)

# --------------------------------------------------------------

#####
# PARSER
#####
# --------------------------------------------------------------
####
# MAIN RUN
####

# If in cluster
if options.run_on_cluster:
    update_progress(options.analysisId, 'Checking files and directory structure', 1)
    check_everything([options.read1, options.read2])
    # Check directory structure
    if options.output_folder[-1]!='/':
        new_output = options.output_folder+'/'
    if not os.path.exists(new_output):
        os.makedirs(new_output)

    # Write bash script
    bash_script = '#!/bin/sh\n\n'
    bash_script += 'python /software/ls/fastqins/fastqins_web.py -i1 '+str(options.read1)+' -i2 '+str(options.read2)
    bash_script += ' -t '+str(options.Tn_seq)+' -g '+str(options.genome)+' -o '+str(new_output)+' -p '+str(options.pcr_dup)+' -v'
    bash_script += ' -id '+str(options.analysisId)+' -r '+str(options.rm_inter)
    fname = new_output+"qsub_fastqins_"+str(options.analysisId)+".sh"
    with open(fname, "w") as fo:
        fo.write(bash_script)

    # Submit to queue
    update_progress(options.analysisId, 'Submitting to queue', 2)
    if ' -w ' not in bash_script:
        if options.pcr_dup!=0:
            cmd = 'qsub -q long-sl7 -l virtual_free=12G,h_rt=48:00:00 -pe smp 8 -N qsub_fastqins_'+str(options.analysisId)+' -o '+new_output+' -e '+new_output+' '+fname
        else:
            cmd = 'qsub -q long-sl7 -N qsub_fastqins_'+str(options.analysisId)+' -o '+new_output+' -e '+new_output+' '+fname
        # cmd = 'qsub -q long-sl7 -l virtual_free=24G,h_rt=12:00:00 -pe smp 8 -N qsub_fastqins_'+str(options.analysisId)+' -o '+new_output+' -e '+new_output+' '+fname
        os.system(cmd)
    else:
        print("ERROR: double web call")
        raise SystemExit
else:
    # This try-except is required to communicate with server if killed by keyboard.
    try:
        # The job entered the queue at this step, send status
        if options.analysisId:
            status = 1
            communicate_with_server(options.analysisId, status)

        # Run
        update_progress(options.analysisId, 'Running FastQins', 3)
        fastqins(read1=options.read1  , read2=options.read2,
                 Tn_seq=options.Tn_seq,
                 genome=options.genome,
                 output_folder=options.output_folder,
                 separate_orientations=options.separate_orientations,
                 pcr_dup=options.pcr_dup, barcode_size=options.barcode_size,
                 mismatches=options.mismatches, extension=options.extension,
                 rm_inter=options.rm_inter, verbose=options.verbose)
        # Check output and send messages if in cluster
        if options.analysisId:
            handle = options.output_folder
            if handle[-1]!='/':
                handle+='/'
            OUTFILES =  glob.glob(handle+'*.qins')
            OUTFILES += glob.glob(handle+'*.log')
            if options.barcode_size:
                OUTFILES += glob.glob(handle+'*.ides')
            print(OUTFILES)
            write_jobid_files(OUTFILES, analysisId=options.analysisId, bc=options.barcode_size)
    except:
        if options.analysisId:
            communicate_with_server(options.analysisId, -1, 'No connection to queue')

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

