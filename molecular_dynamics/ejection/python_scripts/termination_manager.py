#!/bin/usr/python

import os
import sys
import subprocess
import numpy as np
from random import randint

# SAMPLE CONTROL FILE FOR TRANSLATION TERMINATION PROTOCOL
#######################################################################################################
# prot         = 1cli                                                                                 #
# charmm_exec  = /storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm #
# timestp      = 0.015                                                                                #
# fbsolu       = 0.050                                                                                #
# batht        = 310                                                                                  #
# termi_script = charmm_scripts/rnc_term_protocol.inp                                                 #
# top          = multidom/1cli/inpfiles/1cli_rnc_complex.top                                          #
# param        = multidom/1cli/inpfiles/1cli_rnc_complex.prm                                          #
# max_res      = 345                                                                                  #
# ribostruc    = inpfiles/ribo_no_l22.cor                                                             #
# ribopsf      = inpfiles/ribo_no_l22.psf                                                             #
# log_file_dir = multidom/1cli/outfiles/                                                              #
# output_dir   = multidom/1cli/output/                                                                #
#######################################################################################################

      #
########
######### Purpose
########
      #

# This script picks up where cont_syn_mamager.py leaves off. The final .cor and .psf
# written by the continuous synthesis protocol are fed into a CHARMM script that carries
# out iterative rounds of dynamics until the C-terminal residue has exited to ribosome 
# tunnel. Each frame of dynamics is written to a different DCD files, which are then 
# merged to facilitate future Q analysis. The .cor and .psf files output by this script
# are the proper input for simulations of proteins in bulk to assess post-translational
# kinetic trapping.

      #
########
######### Check input command for proper number and type of arguments
########
      #

# check number of arguments
if len(sys.argv) != 3:

        print 'USAGE: python python_scripts/trans_term_protocol.py [control file] [trajectory index]'
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

      #
########
######### Initialize pointers
########
      #

# pdbid of protein
prot         = '2amj'

# trajectory index
traj         = sys.argv[2]

# path to charmm executable
charmm_exec  = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm'

# integration time step in units of ps
timestep     = '0.015'

# frictional coefficient, units of ps**-1
fbsolu       = '0.050'

# temperature for dynamics
batht        = '310'

# script for translation termination simulations
termi_script = 'charmm_scripts/rnc_term_protocol.inp'

# topology file
top          = '2amj/inpfiles/2amj_rnc_complex.top'

# parameter file
param        = '2amj/inpfiles/2amj_rnc_complex.prm'

# final nascent chain length 
max_res      = '193'

# ribosome only coordinates
ribostruc    = 'inpfiles/ribo_no_l22.cor'

# ribosome only protein structure file
ribopsf      = 'inpfiles/ribo_no_l22.psf'

# directory for log files
log_file_dir = 'None'

# directory for CHARMM and VMD output data
output_dir   = 'None'

      # 
########
######### Parse control file
########
      #

# check to make sure the requested cntrl file exists
if os.path.exists(sys.argv[1]) == False:

        print 'Requested control file does not exist.'
        sys.exit()

# load the cntrl file
f_cntrl = open(sys.argv[1])
cntrl = f_cntrl.readlines()
f_cntrl.close()

# load cntrl file data into memory, superceding hard-coded parameters as necessary
for line in cntrl:

        # if this line is commented out
        if line[0] == '#':

                # skip it
                continue

        temp = line.split('=')

        # make sure this line of the control file is formatted correctly
        if len(temp) != 2:

                print 'The following line in '+sys.argv[1]+' does not obey the expected'
                print 'control file format of var_name = var_value'
                print line
                sys.exit()

        var_name = temp[0].strip()
        var_val  = temp[1].strip()

        if var_name == 'prot':

                prot = var_val

        elif var_name == 'charmm_exec':

                charmm_exec = var_val

        elif var_name == 'timestp':

                timestp = var_val

        elif var_name == 'fbsolu':

                fbsolu = var_val

        elif var_name == 'batht':

                batht = var_val

        elif var_name == 'termi_script':

                termi_script = var_val

        elif var_name == 'top':

                top = var_val

        elif var_name == 'param':

                param = var_val

        elif var_name == 'max_res':

                max_res = var_val

        elif var_name == 'ribostruc':
        
                ribostruc = var_val

        elif var_name == 'ribopsf':

                ribopsf = var_val

        elif var_name == 'log_file_dir':

                log_file_dir = var_val

                # make sure log_file_dir exists
                if os.path.exists(log_file_dir) == False:

                        print 'Requested log file directory '+log_file_dir+' does not exist.'
                        sys.exit()

                # generate the log file path as well
                # as paths for dynamics and elong/mini information
                if log_file_dir[-1] == '/':

                        log_file = log_file_dir+traj+'_termi.log'

                elif log_file_dir[-1] != '/':

                        log_file = log_file_dir+'/'+traj+'_termi.log'

                else:

                        print 'The log file path could not be generated.'
                        sys.exit()

                # use the same dynamics output file
                dyna_out = log_file.split('termi.log')[0]+'termi.out'

        elif var_name == 'output_dir':

                output_dir = var_val

                # make sure output_dir ends in '/'
                if output_dir[-1] != '/':

                        output_dir = output_dir + '/'

                # if this path doesn't exist
                if os.path.exists(output_dir) == False:

                        # then we can't very well write data there, can we?
                        print 'The requested output directory '+output_dir+' does not exist.'
                        sys.exit()

        else:

                print 'The string '+var_name+' is not recognized as a command by '+sys.argv[0]
                sys.exit()

# make sure the cntrl file included an updated log_file_dir
if log_file_dir == 'None':

        print 'No log file directory specified in control file.'
        sys.exit()

# make sure the cntrl file included an updated output_dir
if output_dir == 'None':

        print 'No output directory specified in control file.'
        sys.exit()

# check to see if this log_file already exists
if (os.path.exists(log_file) == True):

        # and if it does, print an error message and die
        print 'The log file you attempted to use for this run already exists; as a'
        print 'precaution, you must manually delete or rename that log file so you'
        print 'do not accidentally delete something important!'
        sys.exit()

# add a descriptive header to the log file
with open(log_file, "a") as f_log:
        f_log.write('############################'+'\n')
        f_log.write('# CONTROL FILE INFORMATION #'+'\n')
        f_log.write('############################'+'\n\n')

# and then write each line from the cntrl file to the log file
for line in cntrl:

        with open(log_file, "a") as f_log:
                f_log.write(line)

with open(log_file, "a") as f_log:
        f_log.write('\n')

# add a descriptive header to the log file
with open(log_file, "a") as f_log:
        f_log.write('############################'+'\n')
        f_log.write('# BEGIN TERMINATION SCRIPT #'+'\n')
        f_log.write('############################'+'\n')
        f_log.write('\n')

      #
########
######### Generate and run translation termination command for this trajectory
########
      #

# get the final protein structure file from rnc_syn_dyna_flex.inp
psf = os.popen('ls '+output_dir+traj+'_r'+max_res+'_steps*psf').read().split()[0]

# generate cor file path based on the psf path
cor = psf.split('.psf')[0]+'.cor'

# piece together the translation termination protocol command
termi_cmd = charmm_exec + ' < '+termi_script+' prot='+prot+' traj='+traj+' top='+top+\
            ' param='+param+' psf='+psf+' cor='+cor+' ribopsf='+ribopsf+' ribostruc='+\
            ribostruc+' total='+max_res+' fbsolu='+fbsolu+' outdir='+output_dir+' rand='+\
            str(random_with_N_digits(9))+' batht='+batht+' timestp='+timestp+' >> '+dyna_out

# append some important information about the job to the log file
with open(log_file, "a") as f_log:
        f_log.write('Termination Command   : '+termi_cmd+'\n\n')

# execute the termination command as a child process
subprocess.call(termi_cmd, shell=True)

      #
########
######### Make a tcl script to merge individual dcd files to facilitate Q analysis
########
      #

# we can use the same psf loaded into termi_script
# but we need to generate a list of the dcd files from the termination protocol
dcd_list = os.popen('ls -1v '+output_dir+traj+'_r'+max_res+'_trans_term*dcd').read().split()

# generate a path to a tcl script that will be used to merge the individual dcd files
tcl_script = 'tcl_scripts/'+prot+'_'+traj+'_merge_termi_dcd.tcl'

# check to see if the tcl script requested already exists
if os.path.exists(tcl_script) == True:

        with open(log_file, "a") as f_log:
                f_log.write('The tcl script '+tcl_script+' already exists. Delete it manually to proceed.')

        sys.exit()

with open(log_file, "a") as f_log:
        f_log.write('############################'+'\n')
        f_log.write('# TCL SCRIPT TO MERGE DCDS #'+'\n')
        f_log.write('############################'+'\n')
        f_log.write('\n')

# append commands to the tcl script for this trajectory
count = 1

for dcd_file in dcd_list:

        tcl_line = 'set dcd'+str(count)+' '+dcd_file

        if count == 1:

                tcl_line = 'mol load psf '+psf+' dcd '+dcd_file+'\n'

        else:

                tcl_line = 'mol addfile '+dcd_file+'\n'

        # write tcl script
        with open (tcl_script, "a") as f_tcl:
                f_tcl.write(tcl_line)

        # and record it in the log file, too
        with open (log_file, "a") as f_log:
                f_log.write(tcl_line)

        count += 1

# finish the tcl script by adding the commands to write the merged dcd and to exit VMD when done
with open (tcl_script, "a") as f_tcl:
        f_tcl.write('animate write dcd '+output_dir+traj+'_termi_merged_'+prot+'.dcd'+'\n')
        f_tcl.write('exit')

# and save these lines to the log file, too
with open (log_file, "a") as f_log:
        f_log.write('animate write dcd '+output_dir+traj+'_termi_merged_'+prot+'.dcd'+'\n')
        f_log.write('exit'+'\n')

vmd_cmd = 'vmd -dispdev text -e '+tcl_script

# save the VMD command to the log file
with open (log_file, "a") as f_log:
	f_log.write('\n')
        f_log.write('VMD command           : '+vmd_cmd+'\n')

# submit the command as a child process
subprocess.call(vmd_cmd, shell=True)

       #
#########
########## Delete the extraneous *term1*dcd, *term2*dcd, etc. files since we have merged their contents
#########
       #

# make a list of the dcd files to delete
dcds_to_rm = os.popen('ls -1v '+output_dir+traj+'_*_trans_term*dcd').read().split()

# make sure that we have found the expected number of dcd files to delete;
# should be equal to the value of count
if len(dcds_to_rm) != (count-1):

        with open (log_file, "a") as f_log:
                f_log.write('\n'+str(len(dcds_to_rm))+' DCD files were selected for deletion but '+str(count+1)+'\n')
                f_log.write('DCD files were created during termination protocol.'+'\n')

        sys.exit()

# generate command to remove dcd files
rm_dcd_cmd = 'rm'
for dcd in dcds_to_rm:

        rm_dcd_cmd += ' '+dcd

# append this command to the log file
with open (log_file, "a") as f_log:
        f_log.write('\n')
        f_log.write('Remove DCDs command   : '+rm_dcd_cmd+'\n')

# issue the deletion command
subprocess.call(rm_dcd_cmd, shell=True)

# write a normal termination marker to the log file
with open(log_file, "a") as f_log:
        f_log.write('\n'+'NORMAL TERMINATION'+'\n')
