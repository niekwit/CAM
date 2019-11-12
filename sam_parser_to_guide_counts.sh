#!/bin/bash

## Script to parse bam into computing counts for guide RNAs in a CRISPR experiment

# This script needs to arguments: samfile ($1) and output file name ($2)

# Remove multi-mapped reads: sed '/XS:/d'
# Keep guide name only: cut -f3
# Sort guides: sort
# Count guides and report: uniq -c

cut -f3 $1 | sort | uniq -c > $2

