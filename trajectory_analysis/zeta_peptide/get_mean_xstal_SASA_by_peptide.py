# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:42:11 2021

@author: nissley
"""

import os, sys
import numpy as np

### PURPOSE: Compute the mean SASA of various peptides within 2QCU xstal simulations
###          These values will be used to compute Chi location for these sites

f_file_list = open('xstal_sasa_file_list.txt')
file_list = f_file_list.readlines()
f_file_list.close()

# peptides with 1 < |log2(R/N)| < 6 at t1 = 1 min, t2 = 5 min, t3 = 120 min
# all of these peptides are more exposed in R than N
t1_sig_peptides1 = ['203', '333-354']
t2_sig_peptides1 = ['293', '351', '487', '387']
t3_sig_peptides1 = ['333-354', '125-136', '313', '284-302', '437', '422', '351', '244-254', '293']

all_sig_peptides1 = t1_sig_peptides1 + t2_sig_peptides1 + t3_sig_peptides1

unique_sig_peptides1 = []
for peptide in all_sig_peptides1:
	if peptide not in unique_sig_peptides1:
		unique_sig_peptides1.append(peptide)

peptide_totals = {}
counts = 0

for peptide in unique_sig_peptides1:
	peptide_totals[peptide] = 0.0

for x in file_list:
	data = np.loadtxt('SASA/'+x.strip())
	for peptide in unique_sig_peptides1:
		if '-' in peptide:
			r1, r2 = peptide.split('-')[0:2]
			r1, r2 = int(r1), int(r2)
			resids = np.arange(r1, r2+1, 1, dtype=int)
			peptide_totals[peptide] += np.sum(data[:, resids])
		else:
			peptide_totals[peptide] += np.sum(data[:, int(peptide)])
			
	counts += len(data)
	
	print ('Done extracting SASA data from', x.strip())
	
ofile = open('xstal_sig_peptide1_means.dat', 'w')
for peptide in all_sig_peptides1:
	ofile.write(peptide + ' ' + '%.9f' %(peptide_totals[peptide]/float(counts))+'\n')
	print (peptide + ' ' + '%.9f' %(peptide_totals[peptide]/float(counts)))
ofile.close()