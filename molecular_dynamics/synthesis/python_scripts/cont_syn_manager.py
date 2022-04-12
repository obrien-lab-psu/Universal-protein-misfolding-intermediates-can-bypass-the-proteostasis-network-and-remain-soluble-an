#!/bin/usr/python

import os
import sys
import subprocess
from subprocess import call
import numpy as np
from random import randint
from random import uniform

######################
# PRODUCTION VERSION #
######################

# SUPERCEDED VERSION 9 ON 4/03/18

# SAMPLE CNTRL FILE FOR RUNNING WITH VARIABLE CODON TRANSLATION TIMES
#######################################################################################################
# prot         = 1jys                                                                                 #
# charmm_exec  = /storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm #
# timestp      = 0.015                                                                                #
# fbsolu       = 0.050                                                                                #
# batht        = 310                                                                                  #
# elong_script = charmm_scripts/rnc_syn_elong_mini_flex.inp                                           #
# dynam_script = charmm_scripts/rnc_syn_dyna_flex.inp                                                 #
# top          = flexi_ribo/1jys/inpfiles/1jys_rnc_complex.top                                        #
# param        = flexi_ribo/1jys/inpfiles/1jys_rnc_complex.prm                                        #
# init_psf     = flexi_ribo/1jys/inpfiles/1jys_r1_rnc_complex.psf                                     #
# init_cor     = flexi_ribo/1jys/inpfiles/1jys_r1_rnc_complex.cor                                     #
# init_res     = 1                                                                                    #
# max_res      = 232                                                                                  #
# ribostruc    = inpfiles/ribo_no_l22.cor                                                             #
# ribopsf      = inpfiles/ribo_no_l22.psf                                                             #
# uniform_ta   = 0                                                                                    #
# mrna_seq     = flexi_ribo/1jys/inpfiles/fake_sequence.txt                                           #
# trans_times  = inpfiles/fluit_trans_times_mean_840000.txt                                           #
# log_file_dir = flexi_ribo/1jys/outfiles/                                                            #
# output_dir   = flexi_ribo/1jys/output/                                                              #
#######################################################################################################

# SAMPLE CNTRL FILE FOR RESTARTING A RUN
#######################################################################################################
# prot         = 1jys                                                                                 #
# charmm_exec  = /storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm #
# timestp      = 0.015                                                                                #
# fbsolu       = 0.050                                                                                #
# batht        = 310                                                                                  #
# elong_script = charmm_scripts/rnc_syn_elong_mini_flex.inp                                           #
# dynam_script = charmm_scripts/rnc_syn_dyna_flex.inp                                                 #
# top          = flexi_ribo/1jys/inpfiles/1jys_rnc_complex.top                                        #
# param        = flexi_ribo/1jys/inpfiles/1jys_rnc_complex.prm                                        #
# init_psf     = flexi_ribo/1jys/inpfiles/1_r50_steps100000_1jys.psf                                  #
# init_cor     = flexi_ribo/1jys/inpfiles/1_r50_steps100000_1jys.cor                                  #
# init_res     = 50                                                                                   #
# max_res      = 232                                                                                  #
# ribostruc    = inpfiles/ribo_no_l22.cor                                                             #
# ribopsf      = inpfiles/ribo_no_l22.psf                                                             #
# uniform_ta   = 0                                                                                    #
# mrna_seq     = flexi_ribo/1jys/inpfiles/fake_sequence.txt                                           #
# trans_times  = inpfiles/fluit_trans_times_mean_840000.txt                                           #
# log_file_dir = flexi_ribo/1jys/outfiles/                                                            #
# output_dir   = flexi_ribo/1jys/output/                                                              #
# restart      = 1                                                                                    #
#######################################################################################################

      #
########
######### Purpose
########
      #

# This python script makes a series of calls to CHARMM scripts that
# amount to the commands necessary to elongate, minimize, and run
# dynamics on a ribosome nascent chain complex. A single instance
# of this script corresponds to a single statistically independent
# trajectory.

# This script requires two command line arguments. The first is the 
# control file being used for this trajectory. The second is an integer
# index used to label the trajectory.

      #
########
######### Check input command for proper number and type of arguments
########
      #

# check number of arguments
if len(sys.argv) != 3:

	print 'USAGE: python python_scripts/cont_syn_manager.py [control file] [trajectory index]'
	print 'All arguments are required.'
	sys.exit()

      #
########
######### Function definition(s)
########
      #

# function to generate random seeds
def random_with_N_digits(n):
	range_start = 10**(n-1)
	range_end = (10**n)-1
	return randint(range_start, range_end)

# function to sample a first passage time from an
# exponential distribution with the requested mean
# value; units are 0.015 ps integration time steps
def sample_fpt_dist(mfptaa):
	r1 = uniform(0,1)
	return int(-np.log(r1)*np.float64(mfptaa))

      #
########
######### DEFAULT INPUTS - some may be superseded by control file inputs
########
      #

# pdbid of the protein
prot         = '2amj'

# trajectory index, given by user separate from the control file;
# any number of statistically independent trajectories can be run 
# based on the same control file as long as the value of traj is 
# altered for each run. Unique log files will be written for all trajectories
traj         =  sys.argv[2]

# integration time step in units of ps for CHARMM dynamics
timestp      = 0.015 

# friction coefficient for Langevin dynamics, units of ps**-1
fbsolu       = 0.05

# temperature at which to hold the system bath, units of degrees Kelvin
batht        = 310

# absolute path to charmm executable
charmm_exec  = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm'

# path to CHARMM script for elongation and minimization
elong_script =  'charmm_scripts/rnc_syn_elong_mini_rigid_unif.inp'

# path to CHARMM script for dynamics
dynam_script =  'charmm_scripts/rnc_syn_dyna_rigid_unif.inp'

# path to topology file for ribosome nascent chain complex
top          = '2amj/inpfiles/2amj_rnc_complex.top'

# path to parameter file for ribosome nascent chain complex
param        = '2amj/inpfiles/2amj_rnc_complex.prm'

# path to protein structure file for initial RNC structure
init_psf     = '2amj/inpfiles/19_r40_steps6882_ta840000_2amj.psf'

# path to coordinate file for initial RNC structure
init_cor     = '2amj/inpfiles/19_r40_steps6882_ta840000_2amj.cor'

# initial value of total; CHECK THIS IS CORRECT OFF BY ONE IS SUPER EASY HERE
# the value of init_res MUST be the integer length of the first RNC at which we
# want to run dynamics NOT elongation/minimization
init_res     =  41

# number of residues in the full-length protein
max_res      =  193

# path to ribosome-only protein structure file
ribopsf      = 'inpfiles/minimal_ribo_struc.psf'

# path to ribosome-only coordinate file
ribostruc    = 'inpfiles/minimal_ribo_struc.cor'

# path to file containing codon sequence for this protein
mrna_seq = 'None'

# path to file containing codon:translation time pairs
# default file path is 'None' so we can check, after reading the cntrl file,
# if it has been updated and see if the uniform_ta flag is also invoked
trans_times = 'None'

# flag for whether or not we want to use a constant translation time
# if 1, use uniform translation time. 0 by default
uniform_ta = '0'

# value for the mean translation time to use when uniform_ta = 1
# since uniform_ta = 0 by default, uniform_mfpt = 'None' by default
uniform_mfpt = 'None'

# directory in which a log file named traj+'_log_file.txt' will be placed
log_file_dir = 'None'

# flag for whether or not this is a restart; 0 by default; if restart == 1
# then certain error checks (e.g., log file doesn't already exist) will be
# passed so that the synth.log file can be appended to
restart      = '0'

      # 
########
######### Parse control file
########
      #

# check to see if the control file exists
if os.path.exists(sys.argv[1]) == False:

	# if it doesn't, print this to screen; we can't append it
	# to the log file because without the cntrl file to specify
	# a log file, we might overwrite something. Do not take that chance.
	print 'Requested control file does not exist.'
	sys.exit()

# open the control file, read it into memory, and then close it
f_cntrl = open(sys.argv[1])
cntrl = f_cntrl.readlines()
f_cntrl.close()

# loop over the lines in the cntrl file
for line in cntrl:

	# if this line is commented out
	if line[0] == '#':

		# skip it
		continue

	# split the line on the equals sign
	temp = line.split('=')

	# make sure this line of the control file is formatted correctly
	if len(temp) != 2:

		print 'The following line in '+sys.argv[1]+' does not obey the expected'
		print 'control file format of var_name = var_value'
		print line
		sys.exit()

	# save the left part as the variable/flag name
	var_name  = temp[0].strip()

	# save the right part as the variable/flag value
	var_value = temp[1].strip()

	#print var_name, type(var_name), var_value, type(var_value)

	# now, look for each of the key word variables
	# and update their values as necessary
	if var_name == 'prot':

		prot = var_value

	elif var_name == 'charmm_exec':

		charmm_exec = var_value

	elif var_name == 'timestp':

		timestp = var_value	

	elif var_name == 'fbsolu':

		fbsolu = var_value

	elif var_name == 'batht':

		batht = var_value

	elif var_name == 'elong_script':

		elong_script = var_value

	elif var_name == 'dynam_script':

		dynam_script = var_value

	elif var_name == 'top':

		top = var_value

	elif var_name == 'param':

		param = var_value

	elif var_name == 'init_psf':

		init_psf = var_value

	elif var_name == 'init_cor':

		init_cor = var_value

	elif var_name == 'init_res':

		init_res = int(var_value)

	elif var_name == 'max_res':

		max_res = int(var_value)

	elif var_name == 'ribopsf':

		ribopsf = var_value

	elif var_name == 'ribostruc':

		ribostruc = var_value

	elif var_name == 'restart':

		restart = var_value

		# make sure restart = {'1', '0'}
		if restart != '1' and restart != '0':

			print 'A restart value of '+restart+' is not allowed; choose from 1 or 0.'
			sys.exit()

	# for mrna_seq, read the sequence specified into memory
	elif var_name == 'mrna_seq':

		# make sure the file exists
		if os.path.exists(var_value) == True:

			f_mrna_seq = open(var_value)
			mrna_seq = f_mrna_seq.read().split()[0]
			f_mrna_seq.close()

		else:

			print 'mRNA sequence '+var_value+' does not exist.'
			sys.exit()	

	# for trans_times, read the codon:translation time pairs into memory
	elif var_name == 'trans_times':

		# make sure the file exists
		if os.path.exists(var_value) == True:

			f_trans_times = open(var_value)
			trans_times = f_trans_times.readlines()
			f_trans_times.close()

		else:

			print 'Codon:translation time pair input file '+var_value+' does not exist.'
			sys.exit()

	elif var_name == 'uniform_ta':

		# make sure the user knows uniform_ta is a binary flag
		if (var_value == '1') or (var_value == '0'):

			uniform_ta = var_value

		else:

			print 'uniform_ta is a binary flag, i.e. it must be 0 or 1. A value of '+var_value+' is not allowed'
			sys.exit()

	elif var_name == 'uniform_mfpt':

		uniform_mfpt = var_value

	elif var_name == 'log_file_dir':

		log_file_dir = var_value

		# make sure log_file_dir exists
		if os.path.exists(log_file_dir) == False:

			print 'Requested log file directory '+log_file_dir+' does not exist.'
			sys.exit()

		# generate the log file path as well
		# as paths for dynamics and elong/mini information
		if log_file_dir[-1] == '/':

			log_file = log_file_dir+traj+'_synth.log'

		elif log_file_dir[-1] != '/':

			log_file = log_file_dir+'/'+traj+'_synth.log'

		else:

		        print 'The log file path could not be generated.'
		        sys.exit()

		dyna_out = log_file.split('synth.log')[0]+'dyna.out'
		mini_out = log_file.split('synth.log')[0]+'mini.out'

	elif var_name == 'output_dir':

		output_dir = var_value

		# make sure output_dir ends in '/'
		if output_dir[-1] != '/':

			output_dir = output_dir + '/'

		# if this path doesn't exist
		if os.path.exists(output_dir) == False:

			# then we can't very well write data there, can we?
			print 'The requested output directory '+output_dir+' does not exist.'
			sys.exit()

	else:

		print 'Variable name '+var_name+' is not recognized as a variable by '+sys.argv[0]
		sys.exit()

# check to see if this log_file already exists
if (os.path.exists(log_file) == True) and (restart == '0'):

	# and if it does, print an error message and die
 	print 'The log file you attempted to use for this run already exists; as a'
        print 'precaution, you must manually delete or rename that log file so you'
        print 'do not accidentally delete something important!'
        sys.exit()

# make sure the cntrl file included an updated log_file_dir
if log_file_dir == 'None':

	print 'No log file directory specified in control file.'
	sys.exit()

# print a message to the log file saying this is a restart if it is a restart
if restart == '1':

	with open(log_file, "a") as f_log:
	        f_log.write('\n'+'############################'+'\n')
        	f_log.write('# RESTARTING AT L = '+str(init_res).ljust(7)+'#'+'\n')
	        f_log.write('############################'+'\n\n')

# add a descriptive header to the log file
with open(log_file, "a") as f_log:
	f_log.write('############################'+'\n')
	f_log.write('# CONTROL FILE INFORMATION #'+'\n')
	f_log.write('############################'+'\n\n')

# and then write each line to the log file
for line in cntrl:

        with open(log_file, "a") as f_log:
                f_log.write(line)

with open(log_file, "a") as f_log:
        f_log.write('\n')

      #
########
######### Final quality control for uniform translation rate use case
########
      #

# if uniform_ta = 1, then we must have no mRNA sequence, no codon-specific translation times,
# and we must have a value for uniform_mfpt other than 'None'
if uniform_ta == '1':

            #
        ###### run some final quality controls, mainly just
        ###### to see that correct parameter combination is present in cntrl file
            #

	if uniform_mfpt == 'None':

		with open(log_file, "a") as f_log:
			f_log.write('When running with uniform_ta = 1, you must specify the value of'+'\n')
			f_log.write('uniform_mfpt, the mean first passage time that will be used to '+'\n')
			f_log.write('generate single-molecule first passage times in units of integration'+'\n')
			f_log.write('time steps.'+'\n')

		sys.exit()

	if mrna_seq != 'None':

		with open(log_file, "a") as f_log:
			f_log.write('When running with uniform_ta = 1, you cannot specify a path to an mRNA sequence'+'\n')
			f_log.write('in your control file.'+'\n')

		sys.exit()

	if trans_times != 'None':

		with open(log_file, "a") as f_log:
			f_log.write('When running with uniform_ta = 1, you cannot specify a path to a set of codon:translation'+'\n')
			f_log.write('time pairs in your control file.'+'\n')

		sys.exit()

      #
########
######### Final quality control for variable translation rate use case; if checks passed, parse mrna sequence and trans times into dictionaries
########
      #

# if uniform_ta = 0 is used, then we MUST have a mRNA sequence, a set
# of codon translation times, and uniform_mfpt must = 'None'
if uniform_ta == '0':

	if uniform_mfpt != 'None':

		with open(log_file, "a") as f_log:
			f_log.write('When uniform_ta = 0, uniform_mfpt cannot appear in your cntrl file.'+'\n')

		sys.exit()
		
	if mrna_seq == 'None':

		with open(log_file, "a") as f_log:
			f_log.write('When uniform_ta = 0, you must specify the variable mrna_seq in'+'\n')
			f_log.write('your control file, which gives a path to a file containing the'+'\n')
			f_log.write('translated sequence of the protein of interest'+'\n')

		sys.exit()

	if trans_times == 'None':

		with open(log_file, "a") as f_log:
			f_log.write('When uniform_ta = 0, you must specify the variable trans_times in'+'\n')
			f_log.write('your control file, which gives a path to a file containing codon:'+'\n')
			f_log.write('translation time pairs.'+'\n')

		sys.exit()

            #
        ###### If we pass all of the above quality controls, echo the trans_time sequence
        ###### to the log file and set up dictionaries for codon:translation time mapping
            #

	# make a dictionary to map codons to translation times
	map_codon_to_mfpt = {}

	# add a descriptive header to the log file
	with open(log_file, "a") as f_log:
        	f_log.write('############################'+'\n')
	        f_log.write('# CODON : TRANS TIME PAIRS #'+'\n')
        	f_log.write('############################'+'\n\n')

	# populate the map_codon_to_mfpt dictionary and write the
	# codon/value pairs out to te log file
	for line in trans_times:

        	temp = line.split()
	        map_codon_to_mfpt[temp[0]] = temp[1]

		with open(log_file, "a") as f_log:
			f_log.write(temp[0]+'\t'+temp[1]+'\n')

	with open(log_file, "a") as f_log:
                f_log.write('\n')

	# check to see if the number of nucleotides read in is evenly divisible by 3
	if len(mrna_seq) % 3 != 0:

		with open(log_file, "a") as f_log:
			f_log.write('Number of nucleotides in sequence, '+str(len(mrna_seq))+', is not evenly divisible by 3'+'\n')
			
		sys.exit()

	# check to make sure the number of codons in the mRNA is equal to max_res+1; 
	# that's one codon for each amino acid plus the stop codon
	if (len(mrna_seq)//3) != (int(max_res)+1):

		with open(log_file, "a") as f_log:
			f_log.write('mRNA sequence length of '+str(len(mrna_seq)//3)+' codons does not match the max_res+1 value of '+str(max_res+1)+'\n')
			f_log.write('You must have one codon for each amino acid as well as a stop codon for a total of max_res+1 codons.'+'\n')

		sys.exit()

	# make a dictionary to map codon position, i, to codon type
	map_resid_to_codon = {}

	resid = 0

	# add a descriptive header to the log file
	with open(log_file, "a") as f_log:
        	f_log.write('############################'+'\n')
	        f_log.write('# mRNA SEQUENCE TRANSLATED #'+'\n')
        	f_log.write('############################'+'\n\n')

	# populate map_resid_to_codon with codon positon : codon type pairs
	for nc in range (1, len(mrna_seq), 3):
        	codon = mrna_seq[nc-1]+mrna_seq[nc]+mrna_seq[nc+1]
	        map_resid_to_codon[resid] = codon
        	resid += 1

		with open(log_file, "a") as f_log:
        	        f_log.write(codon)

	with open(log_file, "a") as f_log:
		f_log.write('\n\n')

      #
########
######### Simulate translation of protein via iterative dynamics and elongation/minimization
########
      #

# add a descriptive header to the log file
with open(log_file, "a") as f_log:
	f_log.write('############################'+'\n')
	f_log.write('#   BEGIN SYNTHESIS LOOP   #'+'\n')
	f_log.write('############################'+'\n')

# loop over the desired set of nascent chain lengths
for total in range (init_res, max_res+1):

            #
	###### get the single-molecule first passage time for this nascent chain length
        ###### how this number is generated depends on whether uniform_ta = 1 or 0
            #

	# newline for easy reading
	with open(log_file, "a") as f_log:
		f_log.write('\n')

	# if we are using variable codon translation rates
	if uniform_ta == '0':
	
		# then first identify the codon being decoded at this nascent chain length
		codon = map_resid_to_codon[total]

		# then determine the mean dwell time for this codon type
		codon_specific_mean_fpt = map_codon_to_mfpt[codon]

		# and use the sample_fpt_dist function to get a single-molecule dwell time
		# from an exponential distribution with mean = mean_fpt
		steps = str(sample_fpt_dist(codon_specific_mean_fpt))

		# write the mean fpt and single-molecular fpt used to the log file
		with open(log_file, "a") as f_log:
			f_log.write('Nascent chain length  : '+str(total)+'\n')
			f_log.write('Selected dwell time   : '+steps+'\n')
			f_log.write('Codon in the A-site   : '+codon+'\n')
			f_log.write('Mean dwell time       : '+codon_specific_mean_fpt+'\n')

	# elif we are using a constant mean first passage time
	elif uniform_ta == '1':

		# then simply generate a random selection from an exponential
		# distribution with mean of uniform_mfpt
		steps = str(sample_fpt_dist(int(uniform_mfpt)))

                # write the mean fpt and single-molecular fpt used to the log file
                with open(log_file, "a") as f_log:
                        f_log.write('Selected dwell time   : '+steps+'\n')
			f_log.write('Mean dwell time       : '+uniform_mfpt+'\n')

	# otherwise
	else:

		with open(log_file, "a") as f_log:
			f_log.write('Error in uniform_ta value, value of '+uniform_ta+' is not allowed.'+'\n')

		sys.exit()

	    #
	###### get the .psf/.cor to read in for dynamics; if this is the first nascent chain length simulated,
	###### use the starting structure from cntrl file; otherwise, use previous initstruc .psf/.cor from elong/mini CHARMM script
            #

	# if we are just starting this simulation
	if total == init_res:

		# use the initial structures indicated by the cntrl file
		psf = init_psf
		cor = init_cor

	# otherwise,
	else:

		# use the psf/cor pair printed by the previous subprocess.call() for elongation/minimization
		psf = os.popen('ls -1v '+output_dir+traj+'_r'+str(total)+'_initstruc_'+prot+'.psf').read().split()[-1]
		cor = psf.split('.psf')[0]+'.cor'

	    #
	###### piece together the dynamics command; note that the format is the same for
	###### rigid or flexible ribosomal protein L22 (CHECK THIS)
            #

       	dynam_cmd = charmm_exec+' < '+dynam_script+' prot='+prot+' traj='+traj+' total='+str(total)+\
		     ' outdir='+output_dir+' timestp='+timestp+' fbsolu='+fbsolu+' batht='+batht+' rand='\
               	     +str(random_with_N_digits(9))+' steps='+steps+' top='+top+' param='\
                     +param+' psf='+psf+' cor='+cor+' >> '+dyna_out

        # save the command in a log file
        with open(log_file, "a") as f_log:
                f_log.write('Dynamics command      : '+dynam_cmd+'\n')
        
	# run the dynamics command
	subprocess.call(dynam_cmd, shell=True)

        # print something to line for testing purposes
	print 'Trajectory '+traj+' nascent chain length '+str(total)+' dynamics done.'

	# since we are running dynamics first in this loop, we need to be sure to check
	# if we are in the total = max_res iteration of the loop; the CHARMM script
	# will try to elongate to max_res + 1 if we let it run for that value of total...
	if total == max_res:

		# ...so bust the loop and be done
		break

            #
       	###### generate a command for elongation and minimization; we need to save the final
        ###### structure so it can be read in by the next round of dynamics
       	    #

	psf = os.popen('ls -1v '+output_dir+traj+'_r'+str(total)+'_steps*psf').read().split()[-1]
	cor = psf.split('.psf')[0]+'.cor'

	# generate the CHARMM command
	elong_cmd = charmm_exec+' < '+elong_script+' prot='+prot+' traj='+traj+' total='+str(total+1)+' rand='\
               	     +str(random_with_N_digits(9))+' outdir='+output_dir+' fbsolu='+fbsolu+' top='+top+' param='\
                     +param+' psf='+psf+' cor='+cor+' ribopsf='+ribopsf+' ribostruc='\
	             +ribostruc+' >> '+mini_out

	# save the command in a log file
	with open(log_file, "a") as f_log:
		f_log.write('Elongate/mini command : '+elong_cmd+'\n')

	# run the CHARMM command for elongation of the nascent chain and minimization of the resultant ribosome
	# nascent chain complex
	subprocess.call(elong_cmd, shell=True)

# if we make it this far, write a normal termination indicator to the log file
with open(log_file, "a") as f_log:
	f_log.write('\n'+'NORMAL TERMINATION'+'\n')
