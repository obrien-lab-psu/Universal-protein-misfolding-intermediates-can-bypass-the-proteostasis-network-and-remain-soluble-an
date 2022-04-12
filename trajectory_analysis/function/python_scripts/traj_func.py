#!/bin/usr/python

from datetime import datetime
startTime = datetime.now()
import os
import sys
import numpy as np
import subprocess
from numpy import array
from MDAnalysis import *
import math

"""
 PURPOSE

 Load single- and multi-domain protein functional residue files, parse them, and carry out
 functional residue analysis for each protein.

 Order of operations for this script
 (1) Get command line arguments: [1] database [2] cntrl file for analysis [3] trajectory 
 (2) Load cntrl file. Will need the domain definitions, the location of the data files,
     and information about chi analysis (sigma, sigma_adjust, etc.)
 (3) Load the appropriate database file and find the functional residues for each domain 
     in each protein - easy for single-domain proteins, harder for multi-domain proteins

"""

# usage statement
if len(sys.argv) != 4:

        print 'USAGE: python traj_func.py [1] [2]'
        print '[1]  : path to appropriate .cntrl file'
        print '[2]  : trajectory index for this calculation'
        print '[3]  : calculation type; choose from sm, intf, and both'
        sys.exit()

##############################################################################################

### DEFINE FUNCTIONS

def read_file(filename):

        f_temp = open(filename)
        temp   = f_temp.readlines()
        f_temp.close()
        return temp

def struc_overlap_dmap(coor, fres):

        # coor -> the set of coordinates for the resids identified to be tied to function 
        #         for the domain currently undergoing analysis. MDAnanlysis automatically
        #         sorts these with the lowest resid at the top
        # fres -> the set of residue ids corresponding to the coordinates we have loaded for this frame
        #  
        # N.b., atom selections within MD analysis are automatically ordered (verified in VMD)
        fres = fres.split()
        resids = []
        for resid in fres:
                resids.append(int(resid))
        resids.sort()
        dmap = np.zeros((len(coor), len(coor)))
        ndists = 0

        for i in range (0, len(coor)):
                for j in range (i+1, len(coor)):
                        if resids[j] > (resids[i]+1):
                                dist = np.sqrt((coor[i,0]-coor[j,0])**2.+ \
                                               (coor[i,1]-coor[j,1])**2.+ \
                                               (coor[i,2]-coor[j,2])**2.)
                                dmap[i, j] = dist
                                dmap[j, i] = dist
                                ndists    += 1
                                #print 'dmap: ', i, j, resids[i], resids[j], dist, ndists

        # return the distance map for and number of non-trivial contacts 
        return dmap, ndists

def struc_overlap(ref_dmap, frame_dmap, ref_ndists):

        # ref_dmap   -> output of func_struc_overlap for reference coordinates
        # frame_dmap -> output of func_struc_overlap for current ts coordinates
        # ref_ndists -> number of distances computed for the reference state

        diff_matrix  = np.zeros((len(ref_dmap), len(ref_dmap)))
        diff_matrix.fill(sigma*sigma_adjust)
        diff_matrix  = diff_matrix - abs(ref_dmap - frame_dmap)
        dists_formed = 0

        for i in range (0, len(frame_dmap)):
                for j in range (i+1, len(frame_dmap)):
                        if ref_dmap[i,j] > 0.0:
                                if diff_matrix[i,j] >= 0.0:
                                        dists_formed += 1
                        #print 'chi: ', i, j, ref_dmap[i,j], diff_matrix[i,j], sigma*sigma_adjust, dists_formed
        return float(dists_formed)/float(ref_ndists)

##############################################################################################

### LOAD COMMAND LINE ARGUMENTS AND SET UP SOME OTHER IMPORTANT VARIABLES

db          = read_file('database_results/Database_122_16_July.txt')
cntrl       = read_file(sys.argv[1])
traj        = sys.argv[2]
calc        = sys.argv[3]
xstal       = False
if calc not in ['sm', 'intf', 'both']:
        print calc+' is not recognized as a valid calculation method - choose from sm, intf, and both'
        sys.exit()
##############################################################################################

### LOAD CNTRL FILE DATA

for line in cntrl:
        # if this line is a comment
        if line[0] == '#':
                continue
        line = line.split('=')
        if len(line) != 2:
                print 'The following line in '+sys.argv[1]+' does not obey the expected'
                print 'control file format of var_name = var_value'
                print line
                sys.exit()
        # save the left part as the variable/flag name
        var_name  = line[0].strip()
        # save the right part as the variable/flag value
        var_value = line[1].strip()
        # get path to the reference coordinates
        if var_name == 'ref_coor':
                ref_coor = var_value
                if os.path.exists(ref_coor) == True:
                        pass
                else:
                        print 'Requested reference coordinate set does not exist.'
                        sys.exit()
        #get path to the reference protein structure file
        elif var_name == 'ref_psf':
                ref_psf = var_value
                if os.path.exists(ref_psf) == True:
                        pass
                else:
                        print 'Requested reference protein structure file does not exist.'
                        sys.exit()
        # get the integration time step, assumed to be in units of picoseconds
        elif var_name == 'tstep':
                tstep = float(var_value)
        # get the numer of integration time steps between frames
        elif var_name == 'nsavc':
                steps_per_frame = float(var_value)
        #get path to directory containing co-translational trajectory data
        #elif var_name == 'cot_data_dir':
        #        cot_data_dir = var_value
                # make sure cot_data_dir ends in '/'
        #        if cot_data_dir[-1] != '/':
        #                cot_data_dir += '/'
                # make sure this path exists
        #        if os.path.exists(cot_data_dir) != True:
        #                print 'The selected data directory for your co-translational data, '+cot_data_dir+' does not exist.'
        #                sys.exit()
        # get path to directory containing termination trajectory data
        #elif var_name == 'term_data_dir':
        #        term_data_dir = var_value
                # make sure term_data_dir ends in '/'
        #        if term_data_dir[-1] != '/':
        #                term_data_dir += '/'
                # make sure this path exists
        #        if os.path.exists(term_data_dir) != True:
        #                print 'The selected data directory for your termination data, '+term_data_dir+' does not exist.'
        #                sys.exit()
        # get path to directory containing post-translational trajectory data
        elif var_name == 'pt_data_dir':
                pt_data_dir = var_value
                # make sure pt_data_dir ends in '/'
                if pt_data_dir[-1] != '/':
                        pt_data_dir += '/'
                # make sure this path exists
                if os.path.exists(pt_data_dir) != True:
                        print 'The selected data directory for your termination data, '+pt_data_dir+' does not exist.'
                        sys.exit()
        # get path to directory containing trajectory data
        elif var_name == 'analysis_dir':
                analysis_dir = var_value
                # make sure analysis_dir ends in '/'
                if analysis_dir[-1] != '/':
                        analysis_dir += '/'
                # make sure this path exists
                if os.path.exists(analysis_dir) != True:
                        print 'The selected analysis directory, '+analysis_dir+' does not exist.'
                        sys.exit()
        # get the integer value corresponding to the maxmimum length of the protein
        elif var_name == 'max_res':
                max_res = int(var_value)
        # get the integer value correcsponding to the first nascent chain length to analyze
        elif var_name == 'init_res':
                init_res = int(var_value)
        # get path to CHARMM executable to use for binary to ascii conversions
        elif var_name == 'charmm_exec':
                charmm_exec = var_value
        # flag for whether or not psf and cor files are in binary
        elif var_name == 'binary':
                binary = var_value
        # path to CHARMM script for converting binary to ASCII psf
        elif var_name == 'psf_conv_inp':
                psf_conv_inp = var_value
        # path to CHARMM script for converting binary to ASCII cor
        elif var_name == 'cor_conv_inp':
                cor_conv_inp = var_value
        # value to use for sigma in structural overlap function
        # equivalent to the covalent CG bond distance
        elif var_name == 'sigma':
                sigma = float(var_value)
        # adjustment value to apply to sigma to take into account thermal fluctuations
        elif var_name == 'sigma_adjust':
                sigma_adjust = float(var_value)
        # determine whethere or not this run is for analysis of simulations initiated from 
        # the crystal structure
        elif var_name == 'xstal':
                if var_value == '1':
                        xstal = True
        # for all other cases
        else:
                pass

# grab the pdb id from the analysis_dir path; these are all the same for my analysis so it works
pdbid = analysis_dir.split('/')[-2]

##############################################################################################

### PARSE DATABASE AND FIND INFORMATION FOR THIS PROTEIN

found     = False
temp_fres = ''
fres      = ''

# find the appropriate section based on pdb id matching
for line in db:

        line = line.split(':')
        if line == ['\n']:
                continue
        name, entry = line[0].strip(), line[1].strip()

        if name == 'PDB used in model' and entry == pdbid.upper():
                found = True
        if found == True and calc == 'both' and name == 'Combined active residues (mapped)':
                temp_fres = entry
                break
        if found == True and calc == 'sm' and name == 'Small molecule active residues (mapped)':
                temp_fres = entry
                break
        if found == True and calc == 'intf' and name == 'Interfacial active residues (mapped)':
                temp_fres = entry
                break

if temp_fres == '':
        print 'No database entry found for '+pdbid+' when using calculation mode '+calc
        sys.exit()

# grab the resids from the functional residue definition and format them for
# use as an atom selection command
temp_fres = temp_fres.split(',')

for item in temp_fres:
        item = item.split()
        if item != []:
                res = int(item[0])
                fres += str(res) + ' '

if fres == '0 ':
        print 'No functional residues are identified in the database for '+pdbid+' when using calculation mode '+calc
        sys.exit()
if len(temp_fres) == 1:
        print 'One functional residue, '+str(temp_fres[0])+', found for '+pdbid+' when using calculation mode '+calc
        print 'Analysis cannot be performed with a single residue - aborting.'
        sys.exit()

print 'pdb:', pdbid
print 'Calculation mode:', calc
print 'Functional residues from database:', fres

# reference universe
ref_univ        = Universe(ref_psf, ref_coor, format='CRD')

# get master list of all dcd files to be analyzed
if xstal == True:
        posttr_dcds = os.popen('ls -1v '+pt_data_dir+traj+'_*'+'xstal'+'*dcd').read().split()
elif xstal == False:
        posttr_dcds = os.popen('ls -1v '+pt_data_dir+traj+'_*'+'post_trans'+'*dcd').read().split()
else:

        print 'Unable to determine simulation type from xstal = '+str(xstal)
        sys.exit()

### LOOP OVER DCD FILES AND PERFORM ANALYSIS
if xstal == True:
        outpath = analysis_dir+traj+'_func_'+calc+'_chi_xstal.txt'

else:
        outpath = analysis_dir+traj+'_func_'+calc+'_chi.txt'

if os.path.exists(outpath) == True:
        print outpath+' already exists - delete it to continue'
        sys.exit()

print 'Data written out to:', outpath

# all xstal simulations have an ascii psf available already
if xstal == True:
        traj_psf = ref_psf

# most post-trans simulations were run in binary and so we need to convert the trans_term psf to 
# ascii for use with MDAnalysis
else:

        # path to original (binary) psf
        traj_psf = os.popen('ls -1v '+pt_data_dir+traj+'_*trans_term_*psf' ).read().split()[0]

        if binary == '1':

                # carry out the conversion (place the new PSF in analysis_dir)
                conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+traj_psf+\
                           ' outpsf='+analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                # run the command as a child process
                subprocess.call(conv_cmd, shell=True)

                # get the name of the new PSF file
                traj_psf = analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

current_step = 5000

ref_atomgroup = ref_univ.select_atoms('segid A and (resid '+fres+')')

ref_dmap, ref_ndists = struc_overlap_dmap(ref_atomgroup.positions, fres)
# loop over all dcds file for this analysis run 
for dcd in posttr_dcds:

        current_univ = Universe(traj_psf, dcd)
        print 'Starting calculations for dcd file: '+dcd
        current_atomgroup = current_univ.select_atoms('segid A and (resid '+fres+')')

        for ts in current_univ.trajectory:

                # get the dmap for this frame
                curr_dmap = struc_overlap_dmap(current_atomgroup.positions, fres)[0]

                # and compute chi for this frame using the reference dmap and ndists
                chi_frame = struc_overlap(ref_dmap, curr_dmap, ref_ndists)

                # write out time in ns and chi_functional for this frame
                pline = '%.9f' %(float(current_step)*0.001*0.015)+'\t'+'%.5f' %chi_frame

                with open (outpath, "a") as out:
                        out.write(pline+'\n')
                        
                current_step += 5000

# remove ascii protein structure files
#if binary == '1' and xstal == False:
#        subprocess.call('rm '+analysis_dir+traj+'_*ascii.psf', shell=True)

# print normal termination demarkator
print 'NORMAL TERMINATION'
