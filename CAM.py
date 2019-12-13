#!/usr/bin/python3

import os
import sys

current_path = os.path.realpath(__file__)
cam_directory = os.path.dirname(current_path)
current_path = os.path.dirname(current_path) + '/PRAGUI'

sys.path.append(current_path)
import rnaseq_pip_util as pragui

current_path = current_path + '/cell_bio_util'
sys.path.append(current_path)

import cell_bio_util as util


PROG_NAME = 'CAM'
DESCRIPTION = 'CRISPR Analysis Module'


# Function to convert sam to bam using samtools:
def convert_sam_to_bam(sam,is_single_end):
  if is_single_end:
    bam_file = sam.strip('.sam') + '.bam'
    cmdArgs = ['samtools','view','-bh',sam,'-o',bam_file]
    util.call(cmdArgs)
    #file = bam_file
    os.remove(sam)
    return(bam_file)


# Function to run bowtie
def run_aligner(trimmed_fq,fastq_dirs,aligner='bowtie2',reference_fasta=None,genome_index=None,num_cpu=util.MAX_CORES, is_single_end=True,pair_tags=['r_1','r_2'],aligner_args=None,convert_to_bam=True):
  # Generate genome indexes if not provided
  if aligner == 'bowtie':
    genome_index = os.path.dirname(reference_fasta) + '/bt-genome/'
    index_builder = 'bowtie-build'
  elif aligner == 'bowtie2':
    genome_index = os.path.dirname(reference_fasta) + '/bt2-genome/'
    index_builder = 'bowtie2-build'
  if aligner in [ 'bowtie', 'bowtie2']:
    if genome_index is None:
      util.warn('Folder where %s indices are located hasn\'t been specified. Program will default to %s...' % (aligner,genome_index))
      base = os.path.basename(reference_fasta).split('.')[:-1]
      base = '.'.join(base)
      base = genome_index + base
      if not os.path.exists(genome_index):
        os.mkdir(genome_index)
        util.info('Bowtie2 indices not found. Generating indices...')
        cmdArgs = [index_builder,reference_fasta,base]
        util.call(cmdArgs)
      genome_index = base
    
    # Alignment
    util.info('Aligning reads using %s...' % aligner)
    
    def format_aligner_input(trimmed_fq,aligner,aligner_args,is_single_end):
      if aligner == 'bowtie':
        ext = 'bt'
      elif aligner == 'bowtie2':
        ext = 'bt2'
      k = 0 
      if is_single_end:
        file_list = []
        for f in trimmed_fq:
          cmdArgs = [aligner] + aligner_args
          fo = os.path.basename(f)
          fo = fastq_dirs[k]+ '/' + fo
          sam = fo + '.%s.sam' % ext
          log = fo + '.%s.log' % ext
          return([sam,log])
   
    if aligner == 'bowtie':
      if aligner_args is None:
        aligner_args = ['-v', '0', '-m', '1', '--strata'] # allow no mismatches and report reads that align only once
      sam , log = format_aligner_input(trimmed_fq=trimmed_fq,aligner=aligner,aligner_args=aligner_args,is_single_end=is_single_end)
      file = sam
      if pragui.exists_skip(sam):
        cmdArgs = aligner_args + ['-p',str(num_cpu), genome_index,f, '-S','--sam-nohead', '--no-unal',sam]
        util.call(cmdArgs,stderr=log)
          
    if aligner == 'bowtie2':
      if aligner_args is None:
        aligner_args = ['-N','0','--no-1mm-upfront','-L','25', '--no-unal', '--no-hd'] # set seed to read length, allow no mismatches and no pre-alignment before multiseed heuristic
      sam , log = format_aligner_input(trimmed_fq=trimmed_fq,aligner=aligner,aligner_args=aligner_args,is_single_end=is_single_end)
      file = sam
      if pragui.exists_skip(sam):
        print([aligner,aligner_args])
        cmdArgs = aligner_args + ['-p',str(num_cpu),'-x', genome_index,'-U', f, '-S', sam]
        util.call(cmdArgs,stderr=log)
            
    # Convert sam to bam
    if convert_to_bam is True:
      file = convert_sam_to_bam(sam=sam,is_single_end=is_single_end)
    
    file_list.append(file)
    return(file_list)


# Function to run sam_parser_to_guide_counts.sh and to convert sam files to bam
def sam_parser_parallel(file_list, num_cpu=util.MAX_CORES, remove_sam = True):
  
  util.info('Parsing sam files to get guide counts...')
  
  def sam_parser(sam_file,remove_sam = True):
    ext = '.sam'
    if '.bam' in sam_file:
      ext = '.bam'
    counts_file = sam_file.strip(ext) + '_lib_guidecounts.txt'
    counts_log = sam_file.strip(ext) + '_lib_guidecounts.log'
    # temp = sam_file.strip(ext) + '_temp.sam'
    # Remove header, non-aligned reads and multi-mapped reads from sam/bam file:
    # cmdArgs = ["samtools","view","-q","2","-F","4",sam_file,'-o',temp] # Keep only mapped reads (-F 4). (-q 2) should remove multimapped reads but doesn't work for bt2 given that it doesn't provide the correct flag.  
    #util.call(cmdArgs)
    sam_parser_to_guide_counts = os.path.dirname(os.path.realpath(__file__)) + '/sam_parser_to_guide_counts.sh'
    util.info(sam_parser_to_guide_counts)
    # cmdArgs = [sam_parser_to_guide_counts,temp,counts_file]
    cmdArgs = [sam_parser_to_guide_counts,sam,counts_file]
    util.call(cmdArgs,stderr=counts_log)
    # os.remove(temp)
    if remove_sam is True and ext is '.sam':
      os.remove(sam_file)
    return(counts_file)
    
  common_args=[remove_sam]
  counts_file_list = util.parallel_split_job(sam_parser,file_list,common_args,num_cpu)
  return(counts_file_list)

  
# Wrapper function
def CAM(samples_csv, reference_fasta=None, trim_galore=None, skipfastqc=False, fastqc_args=None, is_single_end=True, pair_tags=['r_1','r_2'], aligner='bowtie2', genome_index=None, aligner_args=None, sam_output='convert_to_bam', multiqc=True, num_cpu=util.MAX_CORES):
  
  convert_to_bam = False
  remove_sam = True
  
  if sam_output not in ['sam','convert_to_bam','delete']:
    util.critical('sam_output flag has been misassigned. Please assign one of the following option: sam, convert_to_bam or delete. For help please type python3 CAM.py --help')
  elif sam_output == 'convert_to_bam':
    convert_to_bam = True
  elif sam_output == 'sam':
    remove_sam = False
  
  if isinstance(pair_tags, str):
    pair_tags = pair_tags.split(',')
  
  header, csv = pragui.parse_csv(samples_csv)
  
  # Fastq file trimming using trimgalore
  # Defaults: --length 12 -e 0.2 -a GTTTAAGAGCTA ( only first 12 bp of adapter GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA)
  if trim_galore is None:
    trim_galore='--length 11 -e 0.2 -a GTTTAAGAGCTA --clip_R1 1' # Allow up to two differences in adapter
  trimmed_fq, fastq_dirs = pragui.trim_bam(samples_csv=samples_csv, csv=csv, trim_galore=trim_galore, skipfastqc=skipfastqc, fastqc_args=fastqc_args, 
                                    is_single_end=is_single_end, pair_tags=pair_tags)

  # Alignment using bowtie2
  # Defaults: --no-sq -5 1 -N 1
  file_list = run_aligner(trimmed_fq=trimmed_fq,fastq_dirs=fastq_dirs,aligner=aligner,reference_fasta=reference_fasta,genome_index=genome_index,num_cpu=num_cpu, is_single_end=is_single_end,pair_tags=pair_tags,aligner_args=aligner_args,convert_to_bam=convert_to_bam)


  # Bam files processing to create input for MAGeCK
  # Run sam_parser_to_guide_counts.sh
  counts_file_list = sam_parser_parallel(file_list=file_list, num_cpu=num_cpu)
  
  # Run Multiqc for quality control 
  pragui.run_multiqc(multiqc=multiqc)



if __name__ == '__main__':

  from argparse import ArgumentParser

  epilog = 'For further help on running this program please email paulafp@mrc-lmb.cam.ac.uk.\n\n'
  epilog += 'Generic example use: python3 CAM.py samples.csv reference.fa \n\n'
  epilog += 'Example use for bassik library: python3 CAM.py -aligner_args "-5 1" samples.csv reference.fa \n\n'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('samples_csv', metavar='SAMPLES_CSV',
                         help='File path of a tab-separated file containing the samples names, the file path for read1, the file path for read2, the experimental condition (e.g. Mutant or Wild-type) and any other information to be used as contrasts for differential expression calling. For single-ended experiments, please fill read2 slot with NA.')

  arg_parse.add_argument('reference_fasta', metavar='REFERENCE_FASTA',
                         help='File path of guide RNAs\' reference sequence FASTA file (for use by genome aligner)')
                         
  arg_parse.add_argument('-trim_galore', # metavar='TRIM_GALORE_OPTIONS',
                         default=None,
                         help='options to be provided to trimgalore. They should be provided under double quotes. If not provided, trimgalore will run with developer\'s default options.')

  arg_parse.add_argument('-fastqc_args', metavar='FASTQC',
                         default=None,
                         help='options to be provided to fastqc. They should be provided under double quotes. If not provided, fastqc will run with developer\'s default options.')

  arg_parse.add_argument('-skipfastqc', default=False, action='store_true',
                         help='Option to skip fastqc step. If this option is set, the option -fastqc_args will be ignored.')

  arg_parse.add_argument('-al', metavar='ALIGNER_NAME', default='bowtie2',
                         help='Name of the program to perform the genome alignment/mapping: Default: bowtie2')
                      
  arg_parse.add_argument('-aligner_index', metavar='BOWTIE_REFERENCE_SEQUENCE_INDEX', default=None,
                         help='Path to directory where bowtie2 indices are stored.')

  arg_parse.add_argument('-aligner_args', default=None,
                         help='Options to be provided to bowtie2. They should be provided under double quotes. If not provided, bowtie2 will be using the following options: -N 1')
  
  arg_parse.add_argument('-sam_output', default='convert_to_bam',
                         help='Specify what to do with the sam file. Options are: sam (keep sam file),convert_to_bam (convert sam file to bam format), delete (delete sam file - best option to save disk space). Default is set to convert_to_bam')

  arg_parse.add_argument('-cpu', metavar='NUM_CORES', default=util.MAX_CORES, type=int,
                         help='Number of parallel CPU cores to use. Default: All available (%d)' % util.MAX_CORES)

  arg_parse.add_argument('-pe', nargs=2, metavar='PAIRED_READ_TAGS', default=['r_1','r_2'],
                        help='The subtrings/tags which are the only differences between paired FASTQ file paths. Default: r_1 r_2')

  arg_parse.add_argument('-se', default=True, action='store_true',
                         help='Input reads are single-end data, otherwise defaults to paired-end.')
  
  arg_parse.add_argument('-disable_multiqc', default=False, action='store_true',
                         help='Specify whether to disable multiqc run. Defaults to False.')

  args = vars(arg_parse.parse_args())

  samples_csv   = args['samples_csv']
  reference_fasta  = args['reference_fasta']
  trim_galore   = args['trim_galore']
  skipfastqc    = args['skipfastqc']
  fastqc_args   = args['fastqc_args']
  aligner       = args['al']
  genome_index     = args['aligner_index']
  aligner_args      = args['aligner_args']
  sam_output    = args['sam_output']
  num_cpu       = args['cpu'] or None # May not be zero
  pair_tags     = args['pe']
  is_single_end = args['se']
  multiqc       = not args['disable_multiqc']
  
  CAM(samples_csv=samples_csv, reference_fasta=reference_fasta, trim_galore=trim_galore, skipfastqc=skipfastqc, fastqc_args=fastqc_args, is_single_end=is_single_end, pair_tags=pair_tags, aligner=aligner, genome_index=genome_index, aligner_args=aligner_args, sam_output=sam_output, multiqc=multiqc, num_cpu=num_cpu)

