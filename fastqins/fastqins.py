#!/usr/bin/env python3

####################################
# 2020 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
#
# Writen by Miravet-Verde, Samuel
# Last updated = 04/28/2020 (MM/DD/YYYY)
#
# Command line example:
# python /software/ls/fastqins/fastqins_web.py -i1 /users/lserrano/smiravet/ins_results/web_tests/test_read1.fastq.gz -i2 /users/lserrano/smiravet/ins_results/web_tests/test_read2.fastq.gz -t TACGGACTTTATC -g /users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk -o /users/lserrano/smiravet/ins_results/web_tests/your_test/ -p 1 -v -r 0
####################################


# --------------------------------------------------
# environment
# --------------------------------------------------
import re            # Only required by python mapper, can be removed without that part
import glob
import os.path
import sys, os
import argparse
import subprocess
import numpy as np
from ruffus import *
from Bio import SeqIO
from collections import Counter

# --------------------------------------------------
# dependencies #
# Edit this section to full paths and reinstall if
# programs not accessible in every part of the system
# --------------------------------------------------

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
            rs['{}\t{}'.format(line[0], line[1])]=float(line[2])
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
Python >=2.7 (>=3 included) |
Author: Miravet-Verde, Samuel |
Last updated: 2020.07.20 |
Affiliation: Center for Genomic Regulation, Luis Serrano's lab |
email: samuelmiver@gmail.com; samuel.miravet@crg.eu |\n\n
"""
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
    # TODO Control memory
    memory = '-m 5G'
    #
    cmd = 'cat {1} {2} | {0} view {3} -b | {0} sort {4} > {5} '.format(samtools, i[0], i[1], align_qual, memory, o)
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
            fo.write('{}\t{}\n'.format(k, v))

def finish_log_file(i, o, tn_seq, genome, flag_pcr=False):
    """
    Adds quality of the process information to the log file
    i = [fastq_dir, fastq_filt, fastq_tn, forward ins, reverse ins, ins, log file]"""
    genome_length = return_genome_lentgh(genome)
    # Compute stats:
    cmd = "awk '{s++}END{print s/4}' "+i[0]
    nr_reads = float(subprocess.check_output(cmd, shell=True).strip())

    cmd = "awk '{s++}END{print s/4}' "+i[1]
    nr_remain = float(subprocess.check_output(cmd, shell=True).strip())
    if flag_pcr:
        filt_line = str(int(nr_remain))+" ("+str(round(100.0*nr_remain/nr_reads, 2))+"%) passed the PCR duplicates filter; of these:"
    else:
        filt_line = "PCR duplicate removal not selected. "+str(int(nr_reads))+" considered; of these:"

    cmd = "cat "+i[2]+" | wc -l"
    nr_IR = float(subprocess.check_output(cmd, shell=True).strip())

    cmd = "cat "+i[3]+" "+i[4]+" | wc -l"
    nr_Tn = float(subprocess.check_output(cmd, shell=True).strip())
    cmd = "cat "+i[5]+ "| wc -l"
    nr_pos = int(subprocess.check_output(cmd, shell=True).strip())
    if nr_pos!=0:
        cmd = "awk '{sum+=$3} END{print sum}' "+i[4]
        nr_post_reads = int(subprocess.check_output(cmd, shell=True).strip())
    else:
        nr_post_reads = 0
    # Write log
    separator = '\n-----------------\n\n'

    with open(i[-1], 'a') as logfile:
        logfile.write(separator+'FASTQINS INFORMATION:\n\n')
        text = "{0} reads provided; of these:\n  ".format(int(nr_reads))+filt_line
        text+= "\n    {0} ({1}%) presented the IR:{2} sequence; of these:\n      ".format(int(nr_IR), round(100.0*nr_IR/nr_remain, 2), tn_seq)
        text+= "{0} having IR ({1}%) were mapped unambiguously; in total:\n        ".format(int(nr_Tn), round(100.0*nr_Tn/nr_IR, 2))
        text+= "{0} insertion positions can be extracted from {1} reads\n".format(nr_pos, nr_post_reads)
        text+= "{0}% of the sample was informative: ".format(round(100.0*nr_Tn/nr_reads, 2))
        text+= "{0}% coverage (number of insertions per genome base).".format(round(100*int(nr_pos)/int(genome_length)))
        logfile.write(text)
        logfile.write(separator+'OUTPUT INFORMATION:\n\n')
        logfile.write('The following files are generated as default output:\n')
        text  =  "\t- *_fw.qins - read counts of insertions mapping to forward strand\n"
        text +=  "\t- *_rv.qins - read counts of insertions mapping to reverse strand\n"
        text +=  "\t- *.qins - read counts of insertions mapping to both strands\n"
        text +=  "\t- *.bam - file generated with the aligned reads. Allows visual revision and additional processing\n"
        text +=  "\t- *.log - log file with general features of the process run\n"
        logfile.write(text)



# --------------------------------------------------
# pipeline control
# --------------------------------------------------

def fastqins(tn_reads, genome,
             paired_reads='0', output_folder='./', tn_seq='TACGGACTTTATC',
             separate_orientations=1, rm_pcr=1, barcode_size=0, mismatches=0, extension='',
             threads=1, align_qual=10,
             ins_calling='bed', zeroes=False, keep_multiple=False,
             rm_inter=0, verbose=False, printout_pipeline=True):
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
                            output     = '{}lis.ls'.format(intermediate_dir)).active_if(rm_pcr).active_if(paired)
    fastqins_pipeline.split(task_func  = rm_pcr_duplicates,
                            input      = create_lis_file,
                            output     = [intermediate_dir+tn_basename+'.filt',
                                          intermediate_dir+paired_reads.split('/')[-1].split('.')[0]+'.filt']).active_if(rm_pcr).active_if(paired)

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

    fastqins_pipeline.merge(task_func  = finish_log_file,
                            input      = [inputtnc,
                                          fastqins_pipeline['tn_trimming'],
                                          fastqins_pipeline['len_mult_selection'],
                                          fastqins_pipeline[ins_calling+'_fw_insertion_calling'],
                                          fastqins_pipeline[ins_calling+'_rv_insertion_calling'],
                                          fastqins_pipeline['count_insertions'],
                                          fastqins_pipeline['create_log_file']],
                            output     = '{}{}.tmplog'.format(intermediate_dir, basename),
                            extras     = [tn_seq, genome, rm_pcr]).follows(fastqins_pipeline['count_insertions'])

    # Printout pipeline
    if printout_pipeline:
        fastqins_pipeline.printout_graph(stream = 'flowchart.svg',
                                         pipeline_name   = 'fastqins | {} mode'.format(ins_calling),
                                         output_format   = 'svg',
                                         no_key_legend   = False,
                                         draw_vertically = True)

#    fastqins_pipeline.transform(task_func  = empty_transform,
#                                input      = fastqins_pipeline['merge_and_count'],
#                                filter     = suffix('.qins'),
#                                output     = ".qins2")


    fastqins_pipeline.run(multiprocess = 4)
#    fastqins_pipeline.run()


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

    print('Moving files...')
    basename = options.tn_reads.split('/')[-1].split('.')[0]
    intermediate_dir = '{}/{}_intermediate_files/'.format(options.output_folder, basename)
    cmd = 'mv {}*.qins {}*.bam {}*.log {}/'.format(intermediate_dir, intermediate_dir, intermediate_dir, options.output_folder)
    print(cmd)
    os.system(cmd)

    if options.rm_inter:
        print('CLEANING...................')
        cmd = 'rm -fR {}'.format(intermediate_dir)
        print(cmd)
        os.system(cmd)
