#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 10:47:22 2020

@author: Niek Wit (MRC LMB)

This function returns a list of guides names extracted from the FASTA CRISPR library file
"""

import pandas as pd

def guide_names(library):
    library = pd.read_csv(library + '.fasta', names = ['guide'])
    #creates new dataframe with only guide names:
    library = library[library['guide'].str.contains('>')]
    #removes '<' character from each row:
    library = library['guide'].str.strip('>')
    library = library.tolist()
    return library