#!/usr/bin/python3

PROG_NAME = 'CAM'
DESCRIPTION = 'CRISPR Analysis Module'

# Fastq file trimming using trimgalore
# Defaults: --length 12 -a GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA

# Alignment using bowtie2
# Defaults: --no-sq -5 1 -N 1


# Bam files processing to create input for MAGeCK
## Remove multi-mapped reads: sed '/XS:/d'
## Keep guide name only: cut -f3
## Sort guides: sort
## Count guides and report: uniq -c


