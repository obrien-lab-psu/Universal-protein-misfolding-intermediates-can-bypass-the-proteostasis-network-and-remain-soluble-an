#!/bin/usr/python

import os
import sys
import numpy as np
import subprocess
from random import randint
import time
from datetime import datetime

      #
########
######### Set up timing functionality
########
      #
# Note on timing. start_time is the time since epoch at program initiation in seconds
#                 current_time is time since epoch at that moment in seconds
#                 elapsed_time is the elapsed since the start of the program in seconds
#                 
#                 elapsed_time = current_time - start_time
#
#                 When restarting a run, I want the script to look in the log file for the
#                 elapsed time after the most-recently completed dynamics iteration. This value
#                 is the total time previously run, and I don't want to run it again. Subtracting 
#                 this previous elapsed time from start_time ensures that the true elapsed time 
#                 since the start of the INITIAL run of the program is maintained. 

# record the start time of the script
start_time = time.time()

# 30 days, 24 hours per day, 60 minutes per hour, 60 seconds per minute
seconds_in_30_days = 30. * 24. * 60. * 60.

      #
########
######### Purpose
########
      #

# This python script makes a series of calls to a CHARMM script in order
# run a specified number of steps of dynamics, write a restart file, and then 
# restart using that file and run another specified number of steps of dynamics.

      #
########
######### Check input command for proper number and type of arguments
########
      #

# check number of arguments
if len(sys.argv) != 3:

        print 'USAGE: python python_scripts/post_trans_manager.py [control file] [trajectory index]'
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
######### Default parameter values
########
      #

# pdbid
prot         = '2amj' 

# trajectory index
traj         = sys.argv[2]

# charmm executable for simulations
charmm_exec  = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm'

# integration time step in uits of ps
timestp      = 0.015

# frictional coefficient, units of ps**-1
fbsolu       = 0.050

# temperature for dynamics
batht        = 310

# script to use for post-translational simulations
post_script  = 'charmm_scripts/test_post_trans_dwell.inp'

# topology file
top          = '2amj/inpfiles/2amj_rnc_complex.top'

# parameter file
param        = '2amj/inpfiles/2amj_rnc_complex.prm'

# size of full-length protein
max_res      = 193

# directory for simulation log files
log_file_dir = 'None'

# directory for simulation output
output_dir   = 'None'

# number of integration time steps to simulate in each dynamics run
steps        = 10000000

# the iteration on which to start the script. If we want to restart
# the run then we need to be able to start at a different location
iteration    = 1

# flag for whether we are restarting this run
restart      = '0'

# flag that indicates whether we are going to be using a post-translational
# conformation or the xstal structure conformation for this run. '0' indicates
# we will use a post-translational conformation found in the output_dir
# with the same trajectory index as indicated in the post_trans_manager.py command
# '1' means that the crystal structure coordinates found at xstal_cor and xstal_psf 
# will be used
xstal        = '0'

# xstal coordinates
xstal_cor    = ''

# xstal protein structure file
xstal_psf    = ''

      #
########
######### Read in cntrl file information
########
      #

if os.path.exists(sys.argv[1]) == False:

        print 'Requested control file '+sys.argv[1]+' does not exist.'
        sys.exit()

# load the cntrl file into memory
f_cntrl = open(sys.argv[1])
cntrl = f_cntrl.readlines()
f_cntrl.close()

for line in cntrl:

        # if this line is a comment
        if line[0] == '#':

                # go to the next line
                continue

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

        elif var_name == 'post_script':

                post_script = var_value

        elif var_name == 'top':

                top = var_value

        elif var_name == 'param':

                param = var_value

        elif var_name == 'steps':

                steps = var_value

        elif var_name == 'max_res':

                max_res = var_value

        elif var_name == 'iteration':

                iteration = int(var_value)

        elif var_name == 'restart':

                restart = var_value

                # make sure restart = {'1', '0'}
                if restart != '1' and restart != '0':

                        print 'A restart value of '+restart+' is not allowed; choose from 1 or 0.'
                        sys.exit()

        elif var_name == 'log_file_dir':

                log_file_dir = var_value

                # make sure log_file_dir exists
                if os.path.exists(log_file_dir) == False:

                        print 'Requested log file directory '+log_file_dir+' does not exist.'
                        sys.exit()

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

        elif var_name == 'xstal':

                xstal = var_value

                if xstal != '1' and xstal != '0':

                        print 'xstal is a binary flag and can only have a value of 1 or 0, so'
                        print 'a value of '+xstal+' is not allowed.'
                        sys.exit()

        elif var_name == 'xstal_cor':

                xstal_cor = var_value

        elif var_name == 'xstal_psf':

                xstal_psf = var_value

        else:

                print 'Variable name '+var_name+' is not recognized as a variable by '+sys.argv[0]
                sys.exit()

# generate the log file path as well
# as paths for dynamics and elong/mini information
# use a different labeling scheme depending on whether
# we are using the xstal structure or a post-translational
# conformation
if log_file_dir[-1] == '/':

        if xstal == '1':

                log_file = log_file_dir+traj+'_xstal.log'

        else:

                log_file = log_file_dir+traj+'_post_trans.log'

elif log_file_dir[-1] != '/':

        if xstal == '1':

                log_file = log_file_dir+'/'+traj+'_xstal.log'

        else:

                log_file = log_file_dir+'/'+traj+'_post_trans.log'

else:

        print 'The log file path could not be generated.'
        sys.exit()

post_out = log_file.split('.')[0]+'.out'

# check to see if this log_file already exists
if (os.path.exists(log_file) == True) and restart != '1':

        # and if it does, print an error message and die
        print 'The log file you attempted to use for this run already exists; as a'
        print 'precaution, you must manually delete or rename that log file so you'
        print 'do not accidentally delete something important!'
        sys.exit()

if restart == '1':

        # append a message to the log file saying so
        with open(log_file, "a") as f_log:
                f_log.write('\n'+'############################'+'\n')
                f_log.write('# RESTARTING AT I = '+str(iteration).ljust(7)+'#'+'\n')
                f_log.write('############################'+'\n\n')        

        # now read in the log file line by line and see what the last elapsed_time call
        # was and make that the new 'start_time' of the script so we keep working on the
        # same 30 days and don't start a new 30 day window
        f_log_file = open(log_file)
        temp_log_file = f_log_file.readlines()
        f_log_file.close()

        # read through line by line
        for line in temp_log_file:

                if 'Elapsed time' in line:

                        # split the line on white space
                        line = line.split()

                        # and update the start time in units of seconds
                        # to be the last recorded elapsed time from the run
                        #start_time = start_time - float(line[3])
                        temp_time = float(line[3])

        start_time = start_time - temp_time
        with open(log_file, "a") as f_log:
                f_log.write('\n'+'Setting internal timer to ' + '%.1f' %temp_time + ' seconds.'+'\n\n')


# make sure the cntrl file included an updated log_file_dir
if log_file_dir == 'None':

        print 'No log file directory specified in control file.'
        sys.exit()

# make sure the cntrl file included an updated output_dir
if output_dir == 'None':

        print 'No output directory specified in control file.'
        sys.exit()

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
######### Final quality control to make sure correct combination of arguments has been entered
########
      #

# make sure that if xstal = '1' we have updated xstal_cor and xstal_psf
if xstal == '1':

        if xstal_cor == '':

                with open(log_file, "a") as f_log:
                        f_log.write('When running with xstal = 1 you MUST specify an appropriate starting coordinate file.'+'\n')
                        f_log.write('i.e., you must provide a file path for the variable xstal_cor in your control file.'+'\n')

                sys.exit()

        if xstal_psf == '':

                with open(log_file, "a") as f_log:
                        f_log.write('When running with xstal = 1 you MUST specify an appropriate starting protein structure file.'+'\n')
                        f_log.write('i.e., you must provide a file path for the variable xstal_psf in your control file.'+'\n')

                sys.exit()

# make sure that if xstal = '0' neither xstal_cor nor xstal_psf have been specified
if xstal == '0':

        if xstal_cor != '':

                with open(log_file, "a") as f_log:
                        f_log.write('When running with xstal = 0 you CANNOT specify a value for xstal_cor in your control file.'+'\n')

                sys.exit()

        if xstal_psf != '':

                with open(log_file, "a") as f_log:
                        f_log.write('When running with xstal = 0 you CANNOT specify a value for xstal_psf in your control file.'+'\n')

                sys.exit()

      #
########
######### Setup and begin restart for post-translational dynamics loop 
########
      #

# add a descriptive header to the log file
with open(log_file, "a") as f_log:
        f_log.write('############################'+'\n')
        f_log.write('#  BEGIN POST-TRANS LOOP   #'+'\n')
        f_log.write('############################'+'\n\n')

# determine the psf to use for this simulation
# this file will be the same throughout all iterations. Use the
# psf file output by termination_manager.py if this is not a xstal
# structure run
if xstal == '1':

        psf       = xstal_psf

        # also, use this label for the output
        job_label = 'xstal'

if xstal == '0':

        psf       = output_dir+traj+'_r'+max_res+'_trans_term_'+prot+'.psf'
        
        # use this label for output trajectories
        job_label = 'post_trans'

# if the psf does not exist then die
if os.path.exists(psf) == False:

        with open(log_file, "a") as f_log:
                f_log.write('Requested psf file '+psf+' does not exist.')

        sys.exit()

# this pointer just keeps the while loop running, it will never actually change
# to false, the script will simply die when it reaches 30 days of walltime and I 
# decide to kill it
keep_running = True

# loop over a counter with an arbitrarily large upper bound
# these simulations will be run for a 30-day walltime and different 
# proteins run at different speeds, so it needs to be sufficiently
# large for the fastest-running proteins.
while keep_running == True:

        # generate a name for the output dcd file for this run
        trajname    = output_dir+traj+'_i'+str(iteration)+'_'+job_label+'_'+prot+'.dcd'

        # generate a name for the output restart file for this run
        res_out     = output_dir+traj+'_i'+str(iteration)+'_'+job_label+'_'+prot+'.res'

        # output cor file for this iteration
        cor_out     = output_dir+traj+'_i'+str(iteration)+'_'+job_label+'_'+prot+'.cor'

        # the same psf is used for all simulations but the cor file will change
        # in the first iteration, use the cor output by the termination script
        if iteration == 1:
        
                if xstal == '0':

                        # input cor file will be, like the input psf, taken from the output of translation termination script
                        cor_in         = psf.split('psf')[0]+'cor'

                else:

                        cor_in = xstal_cor

                # generate CHARMM command
                post_trans_cmd = charmm_exec+' < '+post_script+' traj='+traj+' rand='+str(random_with_N_digits(9))+' top='+top+' param='+\
                                 param+' psf='+psf+' cor='+cor_in+' num='+str(iteration)+' trajname='+trajname+\
                                 ' resout='+res_out+' timestp='+timestp+' batht='+batht+' outcor='+cor_out+\
                                 ' steps='+steps+' >> '+post_out

        # in all other iterations, use the cor file output by the previous iteration
        else:

                # grab the .cor and .res files written at the end of the previous iteration
                cor_in = output_dir+traj+'_i'+str(iteration-1)+'_'+job_label+'_'+prot+'.cor'

                res_in = output_dir+traj+'_i'+str(iteration-1)+'_'+job_label+'_'+prot+'.res'

                # generate CHARMM command, which is slightly different
                post_trans_cmd = charmm_exec+' < '+post_script+' traj='+traj+' rand='+str(random_with_N_digits(9))+' top='+top+' param='+\
                                 param+' psf='+psf+' cor='+cor_in+' num='+str(iteration)+' trajname='+trajname+\
                                 ' resout='+res_out+' timestp='+timestp+' batht='+batht+' outcor='+cor_out+\
                                 ' resin='+res_in+' steps='+steps+' >> '+post_out

                # no need to reseed the random number generator for restarts, use the ISEED from the restart file
                #post_trans_cmd = charmm_exec+' < '+post_script+' traj='+traj+' top='+top+' param='+\
                #                 param+' psf='+psf+' cor='+cor_in+' num='+str(iteration)+' trajname='+trajname+\
                #                 ' resout='+res_out+' timestp='+timestp+' batht='+batht+' outcor='+cor_out+\
                #                 ' resin='+res_in+' steps='+steps+' >> '+post_out

                # note in the log file that this trajectory is being restarted
                with open(log_file, "a") as f_log:
                        f_log.write('Restarting run from   : '+res_in+'\n')

        # append the CHARMM simulation command to the log file
        with open(log_file, "a") as f_log:
                f_log.write('Post-trans command    : '+post_trans_cmd+'\n')
        
        # and then run it on the command line
        subprocess.call(post_trans_cmd, shell=True)

        # append the CHARMM simulation command to the log file
        with open(log_file, "a") as f_log:
                f_log.write('Iteration complete at : '+str(datetime.now())+'\n')

        # check to see if the previous dynamics run gave normal termination by checking the post_trans.out file
        temp_text   = os.popen('tail -n 50 '+post_out).read().split('\n')

        normal_stop = False

        # check line by line for the normal termination text
        for line in temp_text:

                if 'NORMAL TERMINATION BY NORMAL STOP' in line:

                        #print 'Normal termination detected.'
                        normal_stop = True

                        if iteration != 1:

                                # remove the input restart file from this most recent iteration - it is no longer required
                                #subprocess.call('rm '+res_in, shell=True)                        

                                pass

                        # exit this for loop, we have done what we came to do
                        break

        # if normal termination is not detected 
        if  normal_stop == False:
                
                # report this in the log file and then die
                with open(log_file, "a") as f_log:
                        f_log.write('\n'+'Abnormal Termination Detected'+'\n')

                sys.exit()

        # increment the iteration counter
        iteration   += 1

        # check to see if we want to terminate the script following the previous iteration
        current_time = time.time()

        # calculate the time since the start of this iteration of the script
        elapsed_time = current_time - start_time

        # append the CHARMM simulation command to the log file
        with open(log_file, "a") as f_log:
                f_log.write('Elapsed time          : '+str(elapsed_time)+' s, '+str(elapsed_time)+'/'+str(seconds_in_30_days)+'='+str((elapsed_time/seconds_in_30_days)*100.0)+'%'+'\n\n')

        print 'Trajectory '+str(traj)+' iteration '+str(iteration)+' done.'

        # if 30 days or more have elapsed
        if elapsed_time >= seconds_in_30_days:

                # print normal termination to the log file
                with open(log_file, "a") as f_log:
                        f_log.write('Normal Termination'+'\n')

                # and then kill the script   
                sys.exit()

        # otherwise proceed to the next iteration of the script
        else:

                pass
