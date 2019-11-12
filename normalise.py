#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 14:51:11 2019

@author: nwit

Normalises sgRNA count from join file output for direct visual comparisons.

"""
import pandas as pd
import glob

file_list = glob.glob('*.tsv')

#Copies MAGeCK file into dataframe
df = pd.read_csv(file_list[0],sep='\t')

#Creates new dataframe with just guide counts
df1 = df.select_dtypes(include=['integer'])

header_list = []
for col in df1.columns:
    header_list.append(col)
    
header_count = len(header_list)

column_sum = df1.sum()
column_sum2 = list(column_sum)
total_sum = column_sum.sum()

exp_fract = 1/header_count

counter = 0
while counter < header_count:
        for x in column_sum2:
            y = x / total_sum #fraction of total reads
            z = exp_fract / y #normalisation factor
            df[header_list[counter]] = df[header_list[counter]].multiply(z) #normalisation of guide counts
            df[header_list[counter]] = df[header_list[counter]].astype(int)
            counter += 1

df.to_csv('normalised-counts-aggregated.tsv', sep='\t',index=False)