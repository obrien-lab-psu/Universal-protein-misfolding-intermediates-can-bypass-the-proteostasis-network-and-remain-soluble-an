# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 20:33:17 2021

@author: nissley
"""
import sys
import numpy as np

data = np.load('msm_data.npz', allow_pickle=True)

f_ref = open('xstal_sig_peptide1_means.dat')
ref = f_ref.readlines()
f_ref.close()

xstal_ref = {}
for line in ref:
	line = line.strip().split()
	xstal_ref[line[0]] = float(line[1])

state_trajs = data['meta_dtrajs']
states = [0, 1, 2, 3, 4, 5, 6, 7, 8]

trajs = ['1_i45_post_trans_2qcu_sasa.txt.gz', '2_i42_post_trans_2qcu_sasa.txt.gz', '3_i48_post_trans_2qcu_sasa.txt.gz',
		 '4_i43_post_trans_2qcu_sasa.txt.gz', '5_i48_post_trans_2qcu_sasa.txt.gz', '6_i46_post_trans_2qcu_sasa.txt.gz',
		 '7_i43_post_trans_2qcu_sasa.txt.gz', '8_i43_post_trans_2qcu_sasa.txt.gz', '9_i44_post_trans_2qcu_sasa.txt.gz',
		 '10_i53_post_trans_2qcu_sasa.txt.gz', '11_i42_post_trans_2qcu_sasa.txt.gz', '12_i48_post_trans_2qcu_sasa.txt.gz',
		 '13_i42_post_trans_2qcu_sasa.txt.gz', '14_i56_post_trans_2qcu_sasa.txt.gz', '15_i41_post_trans_2qcu_sasa.txt.gz',
		 '16_i43_post_trans_2qcu_sasa.txt.gz', '17_i40_post_trans_2qcu_sasa.txt.gz', '18_i61_post_trans_2qcu_sasa.txt.gz', 
		 '19_i47_post_trans_2qcu_sasa.txt.gz', '20_i44_post_trans_2qcu_sasa.txt.gz', '21_i42_post_trans_2qcu_sasa.txt.gz',
		 '22_i44_post_trans_2qcu_sasa.txt.gz', '23_i43_post_trans_2qcu_sasa.txt.gz', '24_i44_post_trans_2qcu_sasa.txt.gz',
		 '25_i44_post_trans_2qcu_sasa.txt.gz', '26_i45_post_trans_2qcu_sasa.txt.gz', '27_i43_post_trans_2qcu_sasa.txt.gz',
		 '28_i42_post_trans_2qcu_sasa.txt.gz', '29_i44_post_trans_2qcu_sasa.txt.gz', '30_i47_post_trans_2qcu_sasa.txt.gz',
		 '31_i47_post_trans_2qcu_sasa.txt.gz', '32_i43_post_trans_2qcu_sasa.txt.gz', '33_i48_post_trans_2qcu_sasa.txt.gz',
		 '34_i51_post_trans_2qcu_sasa.txt.gz', '35_i49_post_trans_2qcu_sasa.txt.gz', '36_i45_post_trans_2qcu_sasa.txt.gz',
		 '37_i48_post_trans_2qcu_sasa.txt.gz', '38_i42_post_trans_2qcu_sasa.txt.gz', '39_i43_post_trans_2qcu_sasa.txt.gz',
		 '40_i43_post_trans_2qcu_sasa.txt.gz', '41_i49_post_trans_2qcu_sasa.txt.gz', '42_i44_post_trans_2qcu_sasa.txt.gz',
		 '43_i42_post_trans_2qcu_sasa.txt.gz', '44_i45_post_trans_2qcu_sasa.txt.gz', '45_i42_post_trans_2qcu_sasa.txt.gz',
		 '46_i43_post_trans_2qcu_sasa.txt.gz', '47_i43_post_trans_2qcu_sasa.txt.gz', '48_i45_post_trans_2qcu_sasa.txt.gz',
		 '49_i44_post_trans_2qcu_sasa.txt.gz', '50_i49_post_trans_2qcu_sasa.txt.gz']

# peptides with 1 < |log2(R/N)| < 6 at t1 = 1 min, t2 = 5 min, t3 = 120 min
# all of these peptides are more exposed in R than N
t1_sig_peptides1 = ['203', '333-354']
t2_sig_peptides1 = ['293', '351', '487']
t3_sig_peptides1 = ['333-354', '125-136', '313', '284-302', '437', '422', '351', '244-254']
all_sig_peptides1 = t1_sig_peptides1 + t2_sig_peptides1 + t3_sig_peptides1

overlapping_sig_peptides1 = ['333-354', '351', '293']

unique_sig_peptides1 = []
for peptide in all_sig_peptides1:
	if peptide not in unique_sig_peptides1:
		unique_sig_peptides1.append(peptide)

# keys are '0:203' or '5:125-136', signifying 'state:peptide' pairs
# values are total SASA for peptide in state
totals = {}

# same keys as totals
# values are number of frames counted for totals (should be the same for all peptides with the same state)
counts = {}

# populate totals and counts
for state in states:
	for peptide in overlapping_sig_peptides1:
		totals[str(state)+':'+peptide] = 0.0
		counts[str(state)+':'+peptide] = 0

# get mean SASA of each site classified by state
ofile = open('markov_state_zeta_peptide_by_frame.txt', 'w')
for i in range (0, len(trajs)):
	sasa = np.loadtxt('SASA/'+trajs[i])[-1333:, :]
	for peptide in overlapping_sig_peptides1:
		for j in range (0, 1333):
			if '-' in peptide:
				r1, r2 = peptide.split('-')[0:2]
				r1, r2 = int(r1), int(r2)
				resids = np.arange(r1, r2+1, 1, dtype=int)
				val = ((np.sum(sasa[j, resids])/xstal_ref[peptide])-1.0)*100.
				#if val >= 10.0:
					#totals[str(state_trajs[i][j])+':'+peptide] += 1.0
				totals[str(state_trajs[i][j])+':'+peptide] += val
				#ofile.write(str(state_trajs[i][j])+':'+peptide+','+'%.9f' %(val)+'\n')
			else:
				resid = int(peptide)
				val = ((sasa[j, resid]/xstal_ref[peptide])-1.0)*100.
				
				#if val >= 10.0:
				#	totals[str(state_trajs[i][j])+':'+peptide] += 1.0
				totals[str(state_trajs[i][j])+':'+peptide] += val
				#print (sasa[j, resid], xstal_ref[peptide], ((sasa[j, resid]/xstal_ref[peptide])-1.0)*100., state_trajs[i][j])
			ofile.write(trajs[i]+','+str(j)+','+str(state_trajs[i][j])+':'+peptide+','+'%.9f' %(val)+'\n')
			counts[str(state_trajs[i][j])+':'+peptide] += 1
	print ('Done parsing data for', trajs[i])

ofile = open('2qcu_state-wise_sig_peptide1_SASA.dat', 'w')
ofile.write('#state peptide total counts total/counts\n')
for state in states:
	for peptide in overlapping_sig_peptides1:
		print (state, peptide, totals[str(state)+':'+peptide], counts[str(state)+':'+peptide], totals[str(state)+':'+peptide] / float(counts[str(state)+':'+peptide]))
		ofile.write(str(state)+' '+peptide+' '+'%.9f' %totals[str(state)+':'+peptide]+' '+'%.9f' %counts[str(state)+':'+peptide]+' '+'%.9f' %(totals[str(state)+':'+peptide]/float(counts[str(state)+':'+peptide]))+'\n')
ofile.close()