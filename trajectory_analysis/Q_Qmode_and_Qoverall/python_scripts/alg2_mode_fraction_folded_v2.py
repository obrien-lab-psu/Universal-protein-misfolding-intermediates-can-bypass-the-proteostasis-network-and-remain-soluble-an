#!/bin/usr/python

# INTRO

"""
                                   
 Originally written 01/24/2019     
                                   
 Daniel A. Nissley, Ph.D.          
                                   
 O'Brien Lab                       
                                   
 The Pennsylvania State University 

Purpose:  To compute the time-dependent fraction of folded domain/interface based on either Q or chi
          time series produced from simulations of protein synthesis, ejection, and post-translational
          dynamics. 

          This version computes the fraction folded based on comparison to a mode calculated from 
          xstal structure simulations using calc_mode_xstal_values.py

          These values are loaded from xstal_stats.inp

          Updated 02/11/19 to always use bins spread from 0 to 1 with 0.02 widths (will work for all normalized 
                           order parameters, may need to be doctored for other order parameters like RMSD)

          Updated 02/12/19 to bin data based on time edges not number of points. Was a PITA.

          Updated 02/15/19 to use a slightly different algorithm for folded state identification: after identifying
                           an assessment window at which Q_mode(traj) >= <Q_mode(ref)>, it checks the next (tol - 1) 
                           assessment windows. Previously those subsequent windows were all required to be >=
                           <Q_mode(ref)>, but in this version of the algorithm I only require that they be >= the 
                           folding_boundary = <Q_mode(ref)> - gamma*sigma(Q_mode(ref))

"""
#####################################################################################################################

# LOAD MODULES

import os
import sys
import numpy as np
from datetime import datetime
#from decimal import Decimal

# record the start time of the program
start = datetime.now()

#####################################################################################################################

# USAGE STATEMENT

if len(sys.argv) != 5:

        print 'python mode_fraction_folded.py <PDB ID> <structure; d1, dom1_dom5, etc.> <Qi or chi> <trajectory index>'

        sys.exit()

#####################################################################################################################

# function to compute mode based on histogram of data slice
def get_mode(data, bin_width):

        # histogram the input data; always spread bins from 0 to 1 
        bin_counts, bin_edges = np.histogram(data, bins=np.arange(0.0, 1.0 + bin_width, bin_width))

        # if all points in this bin are zero np.histogram returns an empty counts array
        if bin_counts.size == 0:

                # so return a mode of zero
                return 0.0

        # get left edge of most-probable bin
        left_edge             = bin_edges[np.argmax(bin_counts)]

        # get right edge of most-probable bin
        right_edge            = bin_edges[np.argmax(bin_counts)+1]

        # compute the mode as the bin center
        mode                  = (left_edge + right_edge)/2.

        # check to make sure mode is physically reasonable
        if mode > 1. or mode < 0.0:

                print 'A mode of '+'%.4f' %mode+' is not allowed; Q and chi are only'+\
                      'defined on the interval 0 <= (Q or chi) <= 1'

                sys.exit()

        # if so return its value
        else:

                return mode

#####################################################################################################################

# command line arguments

# PDB ID in lower case
pdb         = sys.argv[1].lower()

# struc; d1, dom1_dom3, etc.
struc       = sys.argv[2] 

# Qi or chi
jobtype     = sys.argv[3]

# if fraction native contacts to be used
if jobtype == 'Qi':

        data_column = 3

# if structural overlap function to be used
elif jobtype == 'chi':

        data_column = 4

# if the specified jobtype is not allowed
else:

        print jobtype+' is not an accepted value for jobtype; choose from chi and Qi.'
        sys.exit()

# integer trajectory index
traj        = sys.argv[4]

# other pointers

# number of frames to check following F or U event
tol         = 3

# width of window in ns to consider for each mode computation
window_size = 15. 

# histogram fixed bin width
binwidth    = 0.02

# multiplicative scaling factor; threshold for folding is
# xs_mode - std_scale*xs_std
std_scale   = 3.0

# load xstal reference data from xstal_stats.inp
f_xstal_ref         = open('xstal_stats.inp')
xstal_ref           = f_xstal_ref.readlines()
f_xstal_ref.close()

found_data  = False

for line in xstal_ref:

        line = line.split('\n')[0]

        xs_pdb, xs_struc, xs_jobtype, xs_mode, xs_std = line.split()

        if xs_pdb == pdb and xs_struc == struc and xs_jobtype == jobtype:

                found_data = True

                break

xs_mode, xs_std = float(xs_mode), float(xs_std)

# make sure something was found for this structure in xstal_stats.dat
if found_data == False:

        print 'XSTAL data for '+pdb+' '+struc+' '+jobtype+' not found in xstal_stats.inp'
        sys.exit()

# deallocate memory
xstal_ref    = None

# get trajectory data for this pdb id, structure, and trajectory
temp_data    = np.loadtxt('analysis/'+pdb+'/'+traj+'_'+struc+'_Qi.out')

data         = temp_data[:, data_column]

time         = temp_data[:, 2]

# deallocate
temp_data    = None

# make an array of negative ones sufficiently large to hold all values for any protein
mode_list    = np.ones((2000000, 1))*-1.

bin_centers  = np.ones((2000000, 1))*-1.

left_edges   = np.arange(0.0, 2000000, 0.15, dtype=np.float64)

right_edges  = left_edges + np.float64(15.0)

mode_count   = 0

index1       = 0

index2       = 0

for j in range (0, len(left_edges)):

        left_edge  = left_edges[j]

        right_edge = right_edges[j]

        index1_found = False

        # if this bin would extend past the data
        if right_edge > time[-1]:

                # then we are done
                break

        # loop through the trajectory data
        for i in range (index1, len(data)):

                # if we have found the first point for this bin
                if time[i] >= left_edge and time[i] < right_edge and index1_found == False:

                        # save the index
                        index1 = i

                        index1_found = True

                # or the last point for this bin
                if time[i] > right_edge:

                        # save the index; in this case we have the first index NOT in the current bin
                        index2 = i

                        # ... which is convenient for slicing the array
                        mode_list[mode_count] = get_mode(data[index1:index2], binwidth)

                        # get the center of the bin and save it for later
                        bin_centers[mode_count] = (left_edge + right_edge)/2.

                        #print '%.3f' %((left_edge*2.)+right_edge), '%.3f' %mode_list[mode_count], '%.9f' %time[index1], '%.9f' %time[index2-1], 
                        #print left_edge, right_edge, datetime.now() - start

                        mode_count += 1

                        break

# start trajectory in the unfolded state, as a 1-residue nascent chain is unfolded
# by definition. 2kfw and 3ofo were started at 50 residues, should still hold.
state      = 0

# folding boundary is mean minus one standard deviation
folding_boundary = xs_mode - (std_scale*xs_std)

# generate output file path
outpath = 'analysis/'+pdb+'/'+traj+'_'+struc+'_'+jobtype+'_fF_alg2_'+str(int(std_scale))+'sig_'+str(tol)+'tol.txt'

# and make sure it doesn't already exist
if os.path.exists(outpath) == True:

        print 'Output file already exists.'

        sys.exit()

# loop through mode_list time series and convert to conformational state time series
for s in range (0, len(mode_list)-tol):

        if mode_list[s] < 0.0:

                break

        # if protein unfolded at previous frame
        if state == 0:

                # compare current mode to xs_mode
                if mode_list[s] >= xs_mode:

                        # tentatively change state to 1
                        state = 1

                        # check the next four frames
                        for a in range(s, s+tol):

                                # make sure it remains folded
                                if mode_list[a] < folding_boundary:
        
                                        # or change the state back
                                        # to zero, the fold was not
                                        # stable
                                        state = 0

        # if protein folded at previous frame
        elif state == 1:

                # compare current mode to folding_boundary,
                # which is defined as xs_mode - xs_std

                # if current mode is less than the
                # folding boundary
                if mode_list[s] < folding_boundary:

                        # tentatively change state to 0
                        state = 0

                        # see if the protein goes back above 
                        # folding_threshold in next tol frames
                        for a in range (s, s+tol):

                                # if domain/interface returns to the 
                                # native state within tol frames
                                if mode_list[a] >= folding_boundary:

                                        # then it was never really unfolded
                                        state = 1

        # if state somehow attains a value other than 1 or 0
        else:

                print 'A value of '+str(state)+' is not allowed for the variable \'state\'.'

                sys.exit()

        pline = '%.9f' %bin_centers[s]+'\t'+str(state)+'\t'+'%.5f' %mode_list[s]+'\t'+'%.5f' %xs_mode+\
              '\t'+'%.5f' %(mode_list[s]-xs_mode)+'\t'+'%.5f' %folding_boundary

        with open (outpath, "a") as ofile:

                ofile.write(pline+'\n')

print 'NORMAL TERMINATION'
