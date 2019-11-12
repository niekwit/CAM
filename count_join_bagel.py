#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2019
@author: Niek Wit (MRC-LMB)

This script will join all your individual Bowtie alignment files (.txt) into one file that is suitable for analysis with BAGEL

Intructions:
1. Create a new folder
2. Copy in your library file that only contains gene/sgRNA info (e.g. A1BG_sgA1BG_1) in one column
3. Copy in your Bowtie output .txt files (can be any name) that need to be joined for analysis with BAGEL
    - Important: do not include any other .txt files
4. Copy in script file and run script from terminal
5. The output file (counts-aggregated.tsv) can be used with BAGEL
"""

import glob
import pandas as pd
import csv

#Generates reference list of all library sgRNAs
sgrnas_list00 = list(csv.reader(open("bassik-guides-sorted.csv"))) #Change file name to your library file

sgrnas_list0 = []

for x in sgrnas_list00: #Flattens the list
    for y in x:
        sgrnas_list0.append(y)

#Generates sgRNA and gene columns for final output
sgRNA_output = []
gene_output = []

for n in sgrnas_list0:
    s,g = n.split("_", 1)
    sgRNA_output.append(g)
    gene_output.append(s)

#Generates Pandas data frame from sgRNA list
d0 = {'SEQID':pd.Series(sgRNA_output),'GENE':pd.Series(gene_output),'sgRNA2':pd.Series(sgrnas_list0)}
dfjoin1 = pd.DataFrame(d0)

#Generates a list of all count .txt files
file_list = glob.glob('*.txt')

#Counts number of .txt files in script folder
txtnumber = len(file_list)

#Generates list of lists for join function output
cycle = 1
master_count_list0 = []
while cycle <= txtnumber:
    master_count_list0.append("count_list"+ str(cycle))
    cycle +=1
master_count_list1 = []
for i in master_count_list0:
    master_count_list1.append([i])
    
cycle = 1
master_sgrna_list0 = []
while cycle <= txtnumber:
    master_sgrna_list0.append("sgrna_list"+ str(cycle))
    cycle +=1
master_sgrna_list1 = []
for i in master_sgrna_list0:
    master_sgrna_list1.append([i])

cycle = 1
master_joined_count_list0 = []
while cycle <= txtnumber:
    master_joined_count_list0.append("joined_count_list"+ str(cycle))
    cycle +=1
master_joined_count_list1 = []
for i in master_joined_count_list0:
    master_joined_count_list1.append([i])

#Generates data Pandas frame from count list
counter = 0
while counter < txtnumber:
    #Opens count files and extract counts and sgRNA names
    file = list(csv.reader(open(file_list [counter])))
    
    for x in file:
        a = str(x)
        if a.count(' ') > 1:
            z,b,c = a.split()
            bint = int(b)
            cmod = c.replace("']","")
            master_count_list1 [counter].append(bint)
            master_sgrna_list1 [counter].append(cmod)
        else:
            b,c = a.split()
            bint = b.replace("['","")
            bint = int(bint)
            cmod = c.replace("']","")
            master_count_list1 [counter].append(bint)
            master_sgrna_list1 [counter].append(cmod)
        
    d1 = {'sgRNA2':pd.Series(master_sgrna_list1 [counter]),
        file_list [counter]:pd.Series(master_count_list1 [counter])}
    df1 = pd.DataFrame(d1)

    #Performs left join to merge data sets:
    dfjoin1 = pd.merge(dfjoin1, df1, on='sgRNA2', how='left')
    dfjoin1 = dfjoin1.fillna(0) #Replaces nan with zero
       
    counter +=1

#Deletes sgRNA2 column from dataframe (only needed for joining, not for BAGEL)
dfjoin2 = dfjoin1.drop(columns='sgRNA2')
 
#Writes all data to a single .tsv file, ready for BAGEL
dfjoin2.to_csv('counts-aggregated.tsv', sep='\t',index=False)