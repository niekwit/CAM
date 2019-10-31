#!/usr/bin/python3

import os

current_path = os.path.realpath(__file__)
cam_directory = os.path.dirname(current_path)
current_path = os.path.dirname(current_path) + '/PRAGUI'

sys.path.append(current_path)
import PRAGUI as pragui

PROG_NAME = 'CAM'
DESCRIPTION = 'CRISPR Analysis Module'

# Function to run bowtie
def run_aligner(trimmed_fq,fastq_dirs,aligner='bowtie2',genome_fasta=None,bt2_index=None,num_cpu=pragui.util.MAX_CORES, is_single_end=True,mapq=20,pair_tags=['r_1','r_2'],bt2_args=None):
  if aligner == "bowtie":
    if bt2_index is None:
      bt2_index = os.path.dirname(genome_fasta) + '/bt2-genome/'
      bt2_base = os.path.basename(genome_fasta).split(,'.')[:-1]
      bt2_base = '.'.join(bt2_base)
      pragui.util.info('Bowtie2 indices not found. Generating indices...')
      os.mkdir(bt2_index)
      cmdArgs = ['bowtie2-build',genome_fasta,bt2_base]
      pragui.util.call(cmdArgs)
    
    util.info('Aligning reads using bowtie2...')
    if bt2_args is None:
      bt2_args = ['--no-sq','-5','1','-N','1'] # defaults set by Niek
    cmdArgs = [aligner] + bt2_args
    
    if is_single_end is True:
      cmdArgs = cmdArgs + ['-U'] + trimmed_fq +
                [""]
      
      

# Fastq file trimming using trimgalore
# Defaults: --length 12 -a GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA

# Alignment using bowtie2
# Defaults: --no-sq -5 1 -N 1


# Bam files processing to create input for MAGeCK
# Run sam_parser_to_guide_counts.sh


