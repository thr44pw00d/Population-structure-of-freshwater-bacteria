#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 12:17:11 2021
@author: Chelsea Ramsin
modified by Matthias Hoetzinger on 15.7.2022
"""
import pandas as pd
import os

#Read csv with full dataset including accessions and coordinates
all_organisms = '/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/geo_coord.csv'
df_all_organisms = pd.read_csv(all_organisms, sep=',', index_col=0)

#Define directory in which each organisms result directory lies
dir = '/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/test/'

#Merge pogenom accessions with coordinates by iterating over each organism POGENOM result directory 
for proj in os.listdir(dir):
	df_in = pd.read_csv(dir+proj+'/'+proj+'_cov30_bdth50_samples.txt', sep = ',', index_col=0)
	#Merge accessions from pogenom with the full dataset to get pogenom accessions with coordinates
	merge = pd.merge(df_all_organisms, df_in, how="right", left_index=True, right_index=True)
	merge.to_csv(dir+proj+'/'+proj+'_cov30_bdth50_accessions_coordinates.csv', encoding='utf-8')
