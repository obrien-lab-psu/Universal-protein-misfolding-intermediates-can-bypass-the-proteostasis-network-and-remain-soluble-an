#!/bin/usr/python

             #
            # #   
      ######
################# ## Purpose
           # 
            # #   
             #  

"""
The purpose of this script is to calculate Q for each domain and interface for a given protein 
during the synthesis, termination, and post-translational dynamics phases of simulation. 

Can be used with either normal simulation data or with runs started from the xstal structure

12/14/2018 - Updated to include structural overlap function from Guo and Thirumalai (1994) Biopolymers Eq. 13
             This version only contains functionality to calculate chi for domains, must be updated to be 
             general for interfaces as well.

12/14/2018 - Updated to include structural overlap function for interfaces. This implementation calculates the
             pairwise distances for only those amino acid pairs that are not in the same domain (so it does not
             consider the pairwise distance for amino acids involved in the interface but in the same domain)

             Check to make sure domain-wise structural overlap function works correctly for non-contiguous domains
             If domain 1 is 1-20;31-40, I think the current calculation will not calculate the 20-31 distance because
             in terms of array indices it sees them as "bonded" to one another as they have sequential array indices
             FIX THIS

             I suspect a related issue may be present in the co-translational calculation. When at a nascent chain length
             shorter than the full-length protein I make a truncated ref_chi_dmap called local_chi_dmap which contains
             pairwise distances only for this residues present in the nascent chain. However, when structural_overlap_function
             calculates the delta function is will simply see the zeroes in ref_chi_dmap and still calculate the pairwise distance;
             if the ref_chi_dmap is also zero at that point it may miscount. I need to print these data out and carefully inspect
             the output to see if this is occurring. Could potentially fix this by making a contact_map that defines all of the 
             residue pairs for which the calculation is carried out in the native state. That way you can make sure you are evaluating
             the correct set of residue distances.

12/17/2018   Further updates to chi functionality; changed calculation of chi(frame) to take only to distance maps as inputs,
             one that represents the native state distances and one that represents the current, frame-dependent distances

             Need to update: domain-wise calculation to work correctly for domains with non-contiguous definitions
                             interface calculation to work for inter- and intra- domain interface-involved contacts

             Updated "get_maps" to "get_maps_updated"; only important/usable for domain-wise Q analysis, this new function
             correctly considers contacts between non-contiguous portions of domains (for those domains composed of non-contiguous
             sections of a given protein's primary structure

01/21/2019   Updated chi section to use a corrected algorithm for the continuous synthesis section of the code. This will not be 
             needed for the translation termination or post-translational dynamics portion of the analysis

             The domain-wise chi calculation appears to be working correctly for both contiguous and non-contiguous domains, but
             will require some further verification. For now, moving on. Script passed basic self-consistency check - analysis of
             a trajectory harmonically restrained to the native state produced a chi of 1 for all frames.

             Need to update interface analysis as well and then run for all proteins. Currently running for all single-domain 
             proteins.

02/25/2019   Updated to correct a mistake in structural_overlap_intf_dmap that caused miscounting when some residues constituting
             the domain listed first in the interface definition input file come after some or all of the residues listed second
             in the interface definition file. For example, this bug would cause a problem for the dom3_dom4 interface of 2hnh

             3 1079 1160
             4 558 931

             But not for dom1_dom2 interface of 1cli

             1 1 170
             2 171 345

             I have corrected structural_overlap_intf_dmap so that it now compares the absolute difference between the residues
             under consideration. I do not believe this bug influences the Q{domain}, Q{interface}, or chi{domain} calculations.

             This bug means that I need to rerun interface analysis for 22 proteins: {1d2f, 1duv, 1ger, 1gqe, 1ng9, 1p7l, 1qf6,
             1svt, 1t4b, 1uuf, 1xru, 1xvi, 2hg2, 2hnh, 2kfw, 2kx9, 2ptq, 2wiu, 3ofo, 3pco, 4hr7, 4kn7}. Single-domain proteins 
             are complately unaffected.

08/07/2020   Updated to comment out the chi analysis and fixed a bug that leads to incorrect calculations of co-translational
             order parameters when calculations are started at a nascent chain length other than 1. Only influences {2kfx, 3gn5, and 4g36} 

"""

             #
            # #   
      ######
################# ## Set up modules and record start time
           # 
            # #   
             #  

from datetime import datetime

# record the start time of the program
startTime = datetime.now()

# load other modules
import os
import sys
import numpy as np
import subprocess
from numpy import array
from MDAnalysis import *
import math

if len(sys.argv) != 3:

        print 'USAGE: python trajectory_analysis.py <.cntrl file> <traj index> '
        print 'ALL ARGUMENTS ARE REQUIRED.'
        sys.exit()

# function to open files and read them into memory
def read_file(f_file):

        f_temp = open(f_file)
        temp   = f_temp.readlines()
        f_temp.close()
        return temp

# function to produce distance and contact maps based on an input array
# that is assumed to contain (x,y,z) positions in Angstroms
def get_maps(ref_coor):

        # make an empty distance map of sufficient size for this domain
        dmap           = np.zeros((len(ref_coor), len(ref_coor)))

        # make an empty contact map of sufficient size for this domain
        cmap           = np.zeros((len(ref_coor), len(ref_coor)))

        # counter for total number of contacts within the protein
        total_contacts = 0

        # loop over the set of all coordinates in this domain
        for m in range (0, len(ref_coor)-4):

                # loop over set of residues that could potentially
                # be in (non-trivial) contact with residue m                        
                for k in range (m+4, len(ref_coor)):

                        # calculate the distance
                        dist = ((ref_coor[m,0] - ref_coor[k,0])**2. + \
                                (ref_coor[m,1] - ref_coor[k,1])**2. + \
                                (ref_coor[m,2] - ref_coor[k,2])**2.)**0.5

                        # compare distance to distance cutoff value
                        if dist <= dist_cut:

                                # populate contact map
                                cmap[m,k] = 1
                                cmap[k,m] = 1

                                # populate distance map
                                dmap[m,k]  = dist
                                dmap[k,m]  = dist

                                total_contacts += 1

                        #print m+1, k+1, dist, cmap[m,k], dmap[m,k]

        return dmap, cmap, total_contacts

# updated get_maps function that correctly handles the boundaries between
# different sections of non-contiguous domains
def get_maps_updated(ref_coor, domain_def, num_res_this_domain):

        # make empty dmap and cmap with one extra row and column
        dmap = np.zeros((len(ref_coor), len(ref_coor)))
        cmap = np.zeros((len(ref_coor), len(ref_coor)))

        # make 1D array to hold resid's of all residues in this domain
        # with their number in the original protein numbering scheme
        dom_array = np.zeros((num_res_this_domain, 1))

        # counter for total contacts
        total_contacts = 0

        # convert domain definition residue indices into a list
        temp = domain_def.replace('resid','').replace(':', ' ').replace(';', ' ').split()

        count = 0

        # add an extra column and an extra row to label residues
        for i in range (0, len(temp), 2):

                # start resid  - > end resid
                start, end = int(temp[i]), int(temp[i+1])

                for j in range (start, end+1):

                        dom_array[count] = j

                        count += 1

        # loop over the set of all coordinates in this domain
        for m in range (0, len(ref_coor)):

                # loop over set of residues that could potentially
                # be in (non-trivial) contact with residue m                        
                for k in range (m, len(ref_coor)):

                        if int(dom_array[k]) >= int(dom_array[m])+4:

                                # calculate the distance
                                dist = ((ref_coor[m,0] - ref_coor[k,0])**2. + \
                                        (ref_coor[m,1] - ref_coor[k,1])**2. + \
                                        (ref_coor[m,2] - ref_coor[k,2])**2.)**0.5

                                # compare distance to distance cutoff value
                                if dist <= dist_cut:

                                        # populate contact map
                                        cmap[m,k] = 1
                                        cmap[k,m] = 1

                                        # populate distance map
                                        dmap[m,k]  = dist
                                        dmap[k,m]  = dist

                                        total_contacts += 1

                                #print m, k, dist, cmap[m,k], dmap[m,k], dom_array[m], dom_array[k]

        return dmap, cmap, total_contacts

# function to calculate the number of type i and ii contacts in reference state
def ref_type_i_and_ii_contacts(cmap, dmap, sec_elems):

        # cmap      = contact map for reference state 
        # dmap      = distance map for reference state
        # sec_elems = array of secondary structure definitions

        type_i_ref = 0

        # get type i contacts in this reference state
        for i in range (0, len(sec_elems)):

                for m in range (int(sec_elems[i,1])-1, int(sec_elems[i,2])):

                        for k in range (m+1, int(sec_elems[i,2])):

                                # if these residues are in contact in the native state
                                if cmap[m,k] == 1:

                                        # increment the counter for Type I contacts
                                        type_i_ref += 1

                                #print 'TYPE I:', i, m, k, cmap[m,k], type_i_ref

        type_ii_ref= 0

        # get type ii contacts in this reference state
        for i in range (0, len(sec_elems)):

                for m in range(int(sec_elems[i,1])-1, int(sec_elems[i,2])):

                        for f in range (i+1, len(sec_elems)):

                                for k in range (int(sec_elems[f,1])-1, int(sec_elems[f,2])):

                                        # if these residues are in contact in the native state
                                        if cmap[m,k] == 1:

                                                # increment the counter for Type II contacts
                                                type_ii_ref += 1

                                        #print 'TYPE II:', i, f, m, k, cmap[m,k], type_ii_ref

        return type_i_ref, type_ii_ref

# function to calculate the Q at a particular frame of a DCD file
def frame_type_i_and_ii_contacts(local_cmap, dmap, frame_coor, sec_elems):

        # local_cmap = contact map containing only those contacts possible at the current nascent chain length
        # dmap       = reference (CRITICAL!) distance map for comparisons
        # frame_coor = current coordinates of all extant CG interaction sites in the model
        # sec_elems  = array of secondary structure definitions

        # calculate the Type I contacts for this frame
        type_i_frame = 0

        # loop over secondary structural elements
        for i in range (0, len(sec_elems)):

                # loop over residues in secondary structural element i
                for m in range (int(sec_elems[i,1])-1, int(sec_elems[i,2])):

                        for k in range (m+1, int(sec_elems[i,2])):

                                # if these residues are in contact in the native state
                                if local_cmap[m,k] == 1:

                                        # then compute the distance between them at this frame
                                        dist = ((frame_coor[m,0] - frame_coor[k,0])**2. + \
                                                (frame_coor[m,1] - frame_coor[k,1])**2. + \
                                                (frame_coor[m,2] - frame_coor[k,2])**2.)**0.5

                                        # and if the distance is less than or equal to 
                                        # dmap(m,k)*dist_adjust then the contact is 
                                        # considered to be formed at this frame
                                        if dist <= dmap[m,k]*dist_adjust:

                                                # increment the counter for frame-wise Type I contacts
                                                type_i_frame += 1

        # calculate the Type II contacts for this frame
        type_ii_frame = 0

        for i in range (0, len(sec_elems)):

                for m in range(int(sec_elems[i,1])-1, int(sec_elems[i,2])):

                        for f in range (i+1, len(sec_elems)):

                                for k in range (int(sec_elems[f,1])-1, int(sec_elems[f,2])):

                                        # if these residues are in contact in the native state
                                        if local_cmap[m,k] == 1:

                                                # then compute the distance between them at this frame
                                                dist = ((frame_coor[m,0] - frame_coor[k,0])**2. + \
                                                        (frame_coor[m,1] - frame_coor[k,1])**2. + \
                                                        (frame_coor[m,2] - frame_coor[k,2])**2.)**0.5

                                                # and if the distance is less than or equal to 
                                                # dmap(m,k)*dist_adjust then the contact is 
                                                # considered to be formed at this frame
                                                if dist <= dmap[m,k]*dist_adjust:

                                                        # increment the counter for frame-wise Type II contacts
                                                        type_ii_frame += 1

        return type_i_frame, type_ii_frame

# define function to calculate interface contacts in reference state
def ref_interface_contacts(cmap, dmap, interface_array):

        interface_cont_ref = 0

        # compute total number of interfacial contacts for this interface in the reference state
        for i in range(0, len(interface_array)):

                for m in range(int(interface_array[i,1])-1, int(interface_array[i,2])):

                        for f in range (i+1, len(interface_array)):

                                for k in range (int(interface_array[f,1])-1, int(interface_array[f,2])):

                                        if cmap[m,k] == 1 and int(interface_array[i,0]) != int(interface_array[f,0]):

                                                interface_cont_ref += 1

                                        else:

                                                pass

                                        #print i, m+1, f, k+1, cmap[m,k], interface_cont_ref

        return interface_cont_ref

# define function to calculate interface contacts in a given frame
def frame_interface_contacts(local_cmap, dmap, frame_coor, interface_array):

        interface_cont_i = 0

        # compute total number of interfacial contacts
        for i in range(0, len(interface_array)):

                for m in range(int(interface_array[i,1])-1, int(interface_array[i,2])):

                        for f in range (i+1, len(interface_array)):

                                for k in range (int(interface_array[f,1])-1, int(interface_array[f,2])):

                                        # only bother computing a distance if the local_cmap indicates a contact is necessary
                                        if local_cmap[m,k] == 1 and int(interface_array[i,0]) != int(interface_array[f,0]):

                                                # then compute the distance between them at this frame
                                                dist = ((frame_coor[m,0] - frame_coor[k,0])**2. + \
                                                        (frame_coor[m,1] - frame_coor[k,1])**2. + \
                                                        (frame_coor[m,2] - frame_coor[k,2])**2.)**0.5

                                                # and if the distance is less than or equal to 
                                                # dmap(m,k)*dist_adjust then the contact is 
                                                # considered to be formed at this frame
                                                if dist <= dmap[m,k]*dist_adjust:

                                                        # increment the counter for frame-wise
                                                        # interface contacts

                                                        interface_cont_i += 1
        
        return interface_cont_i

# define function to calculate the distances required by the structural overlap function from Guo and Thirumalai (1994) Biopolymers Eq. 13
# Note that I am omitted the "1 - " factor, such that this metric approaches unity as a protein folds and 0 as it unfolds (completely)
def structural_overlap_dmap(coor, domain_def, num_res_this_dom):

        # make 1D array to hold resid's of all residues in this domain
        # with their number in the original protein numbering scheme
        dom_array = np.zeros((num_res_this_dom, 1))

        # convert domain definition residue indices into a list
        temp = domain_def.replace('resid','').replace(':', ' ').replace(';', ' ').split()

        count = 0

        # counter for total number of pairwise distances considered for this set of coordinates
        dists_possible = 0

        # populate dom_array with resid's
        for i in range (0, len(temp), 2):

                # start resid  - > end resid
                start, end = int(temp[i]), int(temp[i+1])

                #print i, start, end

                for j in range (start, end+1):

                        dom_array[count] = j

                        count += 1

                        #print dom_array[count-1]

        # compute reference structure distances note that this dmap 
        # is distinct from the dmap utilized for Q calculations
        dmap = np.zeros((len(coor), len(coor)))

        # loop over i
        for i in range (0, len(coor)):

                # loop over j 
                for j in range (0, len(coor)):

                        # only consider a pair of residues if they are
                        # not bonded to one another in the protein's
                        # primary amino acid sequence
                        if int(dom_array[j]) > (int(dom_array[i]) + 1):

                                dists_possible += 1

                                # calculate the reference distance
                                dist      = ((coor[i,0] - coor[j,0])**2. + \
                                             (coor[i,1] - coor[j,1])**2. + \
                                             (coor[i,2] - coor[j,2])**2.)**0.5

                                # save the distance to dmap (locally defined)
                                dmap[i,j], dmap[j,i] = dist, dist

                                #print 'chi dmap', i+1, j+1, int(dom_array[i]), int(dom_array[j]), '%.5f' %dmap[i,j], dists_possible

        # output the distance map
        return dmap, dists_possible

# define function that, given a reference and time-dependent coordinate set, will calculate
# the value of the structural overlap function chi from Guo and Thirumalai (1994) Biopolymers Eq. 13
# Note that ref_map and frame_map must be the same size and shape, so that may necessitate making a local_ref
# map that matches the dimensions of the frame_map at a particular nascent chain length
def structural_overlap_function(ref_dmap, frame_dmap, dists_possible):

        # epsilon is the cutoff in the delta function; it is the 
        # covalent bond distance (sigma) multiplied by a factor of 0.2
        # if sigma_adjust == 0.0 then this function tracks the progression
        # of the trajectory towards a single perfectly folded microstate
        diff_matrix     = np.zeros((len(ref_dmap), len(ref_dmap)))

        # fill the entire matrix will 'sigma*sigma_adjust'
        diff_matrix.fill(sigma*sigma_adjust)

        # compute eps minus abs(r_ij(k) - r_ij(ref))
        diff_matrix = diff_matrix - abs(ref_dmap - frame_dmap)

        dists_formed = 0

        # loop over the ref_dmap 
        for i in range (0, len(frame_dmap)):

                for j in range (i+1, len(frame_dmap)):

                        # non-zero elements in the reference dmap
                        # indicate that a contact is possible;
                        # this enforces the selection rules, which are
                        # correctly written into structural_overlap_dmap
                        if ref_dmap[i,j] > 0.0:

                                # then we have a possible contact

                                # now, check to see if the delta function
                                # is zero or greater than 0
                                if diff_matrix[i,j] >= 0.0:

                                        # if the diff matrix is greater than zero
                                        # then we need to count this contact for 
                                        # this domain at this frame
                                        dists_formed += 1 
                        
                        #print 'chi', i, '\t', j,'\t',  '%.5f' %frame_dmap[i,j],'\t', '%.5f' %ref_dmap[i,j],'\t', '%.5f' %diff_matrix[i,j], '\t', dists_formed

        return float(dists_formed)/float(dists_possible)

def structural_overlap_intf_dmap(coor, interface_array, num_res_curr):

        # compute reference structure distances note that this dmap 
        # is distinct from the dmap utilized for Q calculations
        dmap = np.zeros((len(coor), len(coor)))

        # initialize counter for total number of pairwise distances considered in calculation
        contacts_possible = 0

        # get pairwise distances for interface map
        for i in range(0, len(interface_array)):

                for m in range(int(interface_array[i,1])-1, int(interface_array[i,2])):

                        for f in range (i+1, len(interface_array)):

                                for k in range (int(interface_array[f,1])-1, int(interface_array[f,2])):

                                        # only bother computing a distance if the local_cmap indicates a contact is necessary
                                        #if (int(interface_array[i,0]) != int(interface_array[f,0])) and (k > m+1):

                                        # modify this so we check to make sure the absolute difference in the location of 
                                        # the residues being compared in the primary sequence is > 1, i.e. they are not
                                        # bonded to one another
                                        if (int(interface_array[i,0]) != int(interface_array[f,0])) and (abs(m-k) > 1):

                                                # additional conditional to verify final condition:
                                                # this contact must be possible at the current nascent chain length!
                                                #if (m+1) and (k+1) <= num_res_curr:
                                                if (m+1) <= num_res_curr and (k+1) <= num_res_curr:

                                                        # then compute the distance between them at this frame
                                                        dist = ((coor[m,0] - coor[k,0])**2. + \
                                                                (coor[m,1] - coor[k,1])**2. + \
                                                                (coor[m,2] - coor[k,2])**2.)**0.5

                                                        # record the pairwise distance betwen m and k
                                                        dmap[m,k], dmap[k,m] = dist, dist

                                                        # increment counter for number of pairwise distances considered in calculation
                                                        contacts_possible += 1

                                        #print 'intf', m, k, m+1, k+1, '%.5f' %dmap[m,k], contacts_possible

        return dmap, contacts_possible

             #
            # #   
      ######
################# ## Initialize variables and load cntrl file
           # 
            # #   
             #  

# get path to control file for this analysis run
f_cntrl        = sys.argv[1]

# trajectory index for this job
traj           = sys.argv[2]

# read cntrl file into memory
cntrl          = read_file(f_cntrl)

# dictionary to hold domain defintions
domain_defs    = {}

# initialize interface_inp as an empty string
interface_inp  = ''

# initialize pointer to psf binary -> ASCII conversion CHARMM script
psf_conv_inp   = ''

# initialize pointer to cor binary -> ASCII conversion CHARMM script
cor_conv_inp   = ''

# we do not want to skip domains by default
skip_domain    = '0'

# we do not want to skip interfaces by default
skip_interface = '0'

# flag for whether or not the data to be analyzed were initialized from the crystal structure
xstal = False

# loop over the cntrl file
for line in cntrl:

        # if this line is a comment
        if line[0] == '#':

                # skip it
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

        # get path to the reference protein structure file
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

        # get path(s) to secondary structure input files
        elif var_name == 'sec_struc_def':

                sec_struc_def = var_value

        # get path(s) to interface input files
        elif var_name == 'interface_inp':

                interface_inp = var_value

        # get value of the cutoff used to determine whether or not two interaction sites are
        # in contact with one another
        elif var_name == 'dist_cut':

                dist_cut = np.float64(var_value)

        # get value of adjustment factor used to account for thermal fluctuations
        elif var_name == 'dist_adjust':

                dist_adjust = np.float64(var_value)

        # get resid indices for dom1
        elif var_name == 'dom1':

                domain_defs['dom1'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom2
        elif var_name == 'dom2':

                domain_defs['dom2'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom3
        elif var_name == 'dom3':

                domain_defs['dom3'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom4
        elif var_name == 'dom4':

                domain_defs['dom4'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom5
        elif var_name == 'dom5':

                domain_defs['dom5'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom6
        elif var_name == 'dom6':

                domain_defs['dom6'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom7
        elif var_name == 'dom7':

                domain_defs['dom7'] = 'resid '+var_value.replace(';', ' ')

        # get resid indices for dom8
        elif var_name == 'dom8':

                domain_defs['dom8'] = 'resid '+var_value.replace(';', ' ')

        # get value for the total number of domains in the protein
        elif var_name == 'num_dom':

                num_dom = int(var_value)

        # get path to directory containing co-translational trajectory data
        elif var_name == 'cot_data_dir':

                cot_data_dir = var_value

                # make sure cot_data_dir ends in '/'
                if cot_data_dir[-1] != '/':

                        cot_data_dir += '/'

                # make sure this path exists
                if os.path.exists(cot_data_dir) != True:

                        print 'The selected data directory for your co-translational data, '+cot_data_dir+' does not exist.'

                        sys.exit()

        # get path to directory containing termination trajectory data
        elif var_name == 'term_data_dir':

                term_data_dir = var_value

                # make sure term_data_dir ends in '/'
                if term_data_dir[-1] != '/':

                        term_data_dir += '/'

                # make sure this path exists
                if os.path.exists(term_data_dir) != True:

                        print 'The selected data directory for your termination data, '+term_data_dir+' does not exist.'

                        sys.exit()

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

        # get binary flag for whether or not domain analysis is to be carried out
        elif var_name == 'skip_domain':

                skip_domain = var_value

                if skip_domain != '0' and skip_domain != '1':

                        print 'skip_domain must be either 1 or 0, '+skip_domain+' is not an allowed value.'

                        sys.exit()

        # get binary flag for whether or not interface analysis is to be carried out
        elif var_name == 'skip_interface':

                skip_interface = var_value

                if skip_interface != '0' and skip_interface != '1':

                        print 'skip_interface must be either 1 or 0, '+skip_domain+' is not an allowed value.'

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

        elif var_name == 'xstal':

                if var_value == '1':

                        xstal = True

        # for all other cases
        else:

                print var_name + ' is not currently recognized as an input parameter.'
                sys.exit()

             #
            # #   
      ######
################# ## Run some checks and make sure we do not have conflicting information / flags
           # 
            # #   
             #  

# change skip domain to '1' for all proteins; should only be used for reanalysis of multi-domain proteins
# with updated interface chi function that works correctly
#skip_domain = '1'

# perform some error checks; if we have binary data then we need to supply paths to
# CHARMM scripts to carry out binary -> ASCII conversion
if binary == '1':

        if psf_conv_inp == '':

                print 'You must supply a path to a CHARMM script to carry out binary -> ASCII'
                print 'conversion for psf files if binary = 1'

                sys.exit()

        if cor_conv_inp == '':

                print 'You must supply a path to a CHARMM script to carry out binary -> ASCII'
                print 'conversion for cor files if binary = 1'

                sys.exit()

# initialize dictionary relating resid to domain
map_resid_to_domain = {}

# populate map_resid_to_domain based on domain definitions stored in domain_defs
for dom in range (1, num_dom+1):

        temp = domain_defs['dom'+str(dom)]
        temp = temp.strip('resid')
        temp = temp.split()

        for item in temp:

                # get starting and ending resids for this
                # segment of the current domain
                r1 = int(item.split(':')[0])
                r2 = int(item.split(':')[1])

                # loop over domain segment resids
                for res in range (r1, r2+1):

                        # populate dictionary
                        map_resid_to_domain[res] = dom

# generate a list of the secondary structure input files domain(s)
sec_struc_list = sec_struc_def.split(';')

# make sure each of the secondary structure definition files exists
for sec_def_file in sec_struc_list:

        temp = os.path.exists(sec_def_file.strip())

        # if file does not exist
        if temp == False:

                print 'Secondary structure definition file '+sec_def_file.strip()+' does not exist.'
                sys.exit()

# check to make sure we have one secondary structure input file for each domain
if len(sec_struc_list) != num_dom:

        print 'Number of secondary structure files does not match number of domains;'
        print len(sec_struc_list), '!=', num_dom

        sys.exit()

# if interface analysis is required
if skip_interface == '0':

        # generate list of interface input files for domain interfaces
        # make sure that interface_inp contains some information
        if interface_inp != '':

                interface_inp_list = interface_inp.split(';')

                # make sure each interface input file exists
                for interface_def in interface_inp_list:

                        temp = os.path.exists(interface_def.strip())

                        # if file does not exist
                        if temp == False:

                                print 'Interface definition file '+interface_def.strip()+' does not exist.'
                                sys.exit()

        else:

                if num_dom != 1:

                        print "Interface analysis requested, but no interface definition files were loaded."
                        sys.exit()

             #
            # #   
      ######
################# ## Calculate fraction of native contacts for each domain
           # 
            # #   
             #  

# define the reference coordinate universe for this protein
ref_univ       = Universe(ref_psf, ref_coor, format='CRD')

# switch to tell us whether or not we have already completed a loop through the set of 
# dcd files for this protein; used to determine whether or not we need to make 
# binary -> ASCII conversions
already_looped = 0

# loop over each of the domains and perform fraction of native contact analysis
for dom in range (1, num_dom + 1):

        # skip domain calculation if instructed to do so in cntrl file
        if skip_domain == '1':

                print 'Skipping domain analysis as requested in .cntrl file.'

                break

        print '\n'+'Beginning co-translational Q calculations for domain '+str(dom)+'\n'

        # make an AtomGroup for this domain based on the domain definitions
        domain_ref_atomgroup = ref_univ.select_atoms(domain_defs['dom'+str(dom)])

        # calculate the distance and contact maps for the native state coordinates
        # N.B., we need to use the updated "get_maps_updated" function to generate
        # contact and distance maps for domains; get_maps has a bug that causes it 
        # to generate incomplete contact maps for domains with non-contiguous defs
        dmap, cmap, total_contacts = get_maps_updated(domain_ref_atomgroup.positions,\
                                     domain_defs['dom'+str(dom)], len(domain_ref_atomgroup.positions))

        #load the secondary structure input file for this domain
        temp_sec_elems       = read_file(sec_struc_list[dom-1])

        # make an empty numpy array
        sec_elems            = np.zeros((len(temp_sec_elems), 3))

        # counter
        count                = 0

        # and parse it into an array
        for elem in temp_sec_elems:

                # split the line on white space
                elem                = elem.split()

                # populate secondary structure element matrix
                sec_elems[count, 0] = elem[0]
                sec_elems[count, 1] = elem[1]
                sec_elems[count, 2] = elem[2]

                # increment the counter for number of lines
                count              += 1

        # determine the number of type I contacts in the reference state
        type_i_ref, type_ii_ref = ref_type_i_and_ii_contacts(cmap, dmap, sec_elems)

        # make the domain reference dmap for structural_overlap_function
        # requires the reference coordinates, the domain definition for the current domain in terms of resids, and the total
        # number of residues in this domain (total, including all non-contiguous segments). Also, record the total number of pairwise
        # distances possible for the reference state
        #ref_chi_dmap, ref_chi_dists = structural_overlap_dmap(domain_ref_atomgroup.positions, domain_defs['dom'+str(dom)],len(domain_ref_atomgroup.positions))

        # if there are no contacts in the reference state
        if type_i_ref == 0 and type_ii_ref == 0:

                # then move onto the next domain
                # this occurs for PDB ID: 1FTS which has 
                # an intrinisically disordered N-terminal domain
                continue

                     #
                    # #   
              ######
        ################# ## Loop over DCD files and over frames in each DCD file
                   # 
                    # #   
                     #  

        # need to keep track of the overall time for this trajectory
        # initialize current_time at the value of steps_per_frame, as
        # data are printed only after nsavc steps have been simulated, 
        # not at step 0
        current_step                     = steps_per_frame

        total_steps_after_last_nc_len    = 0

        # generate a path for output for this protein, trajectory, and domain
        props_out_path                   = analysis_dir+traj+'_d'+str(dom)+'_Qi.out'

        if xstal == True:

                props_out_path = analysis_dir+traj+'_d'+str(dom)+'_Qi_xstal.out'

        # check to see if this path exists already and die if it does
        if os.path.exists(props_out_path) == True:

                print 'The procedurally generated output file for this protein and trajectory of '+props_out_path+' already exists.'
                print 'Manually delete it to continue.'

                sys.exit()

        # phase simulation is in; 1 -> synthesis
        #                         2 -> termination
        #                         3 -> post-translational
        simulation_phase = 1

	
        #if xstal == False:
        #        # print a time point for time zero; needed for binning
        #        outline = '%.0f' %simulation_phase+'\t'+'%.0f' %float(1.0)+'\t'+'%.9f' %float(0.0)+'\t'+'%.4f' %float(0.0)+'\t'+'%.5f' %float(0.0)+'\n'
        #        # write time point zero to outfile
        #        with open (props_out_path, "a") as props:
        #                props.write(outline)

        # counter for the number of residues present in the nascent chain for this domain
        num_res_this_dom                 = 0

        # added 08/07/2020 by Dan Nissley
        if init_res != 1:
                for r in range (1, init_res): # only want to check up to residue BEFORE init_res, 
                                              # because a check is performed for init_res in the loop over 
                                              # nascent chain length
                        if map_resid_to_domain[r] == dom:
                                num_res_this_dom += 1

        # loop over the set of dcd files for this trajectory and job type        
        for nc_length in range (init_res, max_res+1):

                # there is no continupus synthesis phase for the xstal simulations
                # so skip this section as necessary
                if xstal == True:

                        break

                print 'Starting calculation for domain '+str(dom)+' nascent chain length '+str(nc_length)+' ... '

                # determine the domain in which the new C-terminal bead falls
                if map_resid_to_domain[nc_length] == dom:

                        # if the newest bead falls in this domain, increment counter
                        num_res_this_dom += 1

                # get the cmap and dmap for only the contacts possible at this nascent chain length
                local_cmap            = np.zeros((len(cmap), len(cmap)))

                # only populate local_cmap to the index corresponding to num_res_this_dom
                for a in range (0, num_res_this_dom):

                        for b in range (0, num_res_this_dom):

                                # local contact map for fraction native contacts
                                local_cmap[a, b]     = cmap[a, b]
                                local_cmap[b, a]     = cmap[b, a]

                # generate list of the the dcd file for this nascent chain length and trajectory
                dcd_file_list         = os.popen('ls -1v '+cot_data_dir+str(traj)+'_r'+str(nc_length)+'_steps*dcd').read().split()

                # make sure there is exactly 1 DCD file for this trajectory and nc length
                if len(dcd_file_list) != 1:

                        print 'ERROR IN DCD FILE SELECTION.'
                        print 'CHECK TRAJ '+str(traj)+' NC LENGTH '+str(nc_length)

                        sys.exit()

                # get dcd_file as the only element of dcd_file_list
                dcd_file              = dcd_file_list[0]

                # determine the number of steps we expect to find in this dcd file
                temp                  = dcd_file.split('/')[-1]

                steps_this_nc_len     = int(temp.split('_')[2].split('steps')[1])

                # check to see if this DCD file is empty
                file_size_zero        = os.stat(dcd_file).st_size == 0

                # counter for number of frames in this dcd file
                frames_this_dcd       = 0

                # if the DCD file contains data
                if file_size_zero != True:

                        # generate a path to a PSF file
                        psf_file = dcd_file.split('dcd')[0]+'psf'

                        # check to see if a binary to ascii conversion is necessary
                        if binary == '1' and already_looped == 0:

                                # if so, carry out the conversion (place the new PSF in analysis_dir)
                                conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+psf_file+\
                                           ' outpsf='+analysis_dir+psf_file.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                                # run the command as a child process
                                subprocess.call(conv_cmd, shell=True)

                                # get the name of the new PSF file
                                psf_file = analysis_dir+psf_file.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                print 'ASCII psf '+psf_file+' generated.'

                        # if we do not need to make the binary to ascii conversion
                        elif binary == '1' and already_looped == 1:

                                # get the name of the new PSF file
                                psf_file = analysis_dir+psf_file.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                        # otherwise
                        else:

                                # do nothing
                                pass

                        # load DCD file into PSF to create Universe
                        current_univ = Universe(psf_file, dcd_file)

                        # select SEGID A atoms in the current domain
                        domain_atomgroup = current_univ.select_atoms('segid A and ('+domain_defs['dom'+str(dom)]+')')

                        # loop over frames in the trajectory file
                        for ts in current_univ.trajectory:

                                # calculate the Type I and II contacts for this frame
                                type_i_frame, type_ii_frame = frame_type_i_and_ii_contacts(local_cmap, dmap, domain_atomgroup.positions, sec_elems)

                                # calculate pairwise distances for this frame
                                #temp_frame_chi_dmap, frame_chi_dists = structural_overlap_dmap(domain_atomgroup.positions, domain_defs['dom'+str(dom)],\
                                #                                       len(domain_ref_atomgroup.positions))

                                # reshape temp_frame_chi_dmap so that it is the same shape as the ref_chi_dmap
                                #frame_chi_dmap = np.zeros((len(domain_ref_atomgroup.positions), len(domain_ref_atomgroup.positions)))

                                # populate frame_chi_dmap with data from temp_frame_chi_dmap, but let all elements not populated at this
                                # nascent chain length be zero 
                                #for i in range (0, len(temp_frame_chi_dmap)):

                                #        for j in range (0, len(temp_frame_chi_dmap)):

                                #                frame_chi_dmap[i,j] = temp_frame_chi_dmap[i,j]

                                #                frame_chi_dmap[j,i] = temp_frame_chi_dmap[j,i]

                                # compute the structural overlap function for this domain
                                #chi_frame = structural_overlap_function(ref_chi_dmap, frame_chi_dmap, ref_chi_dists)

                                # manage time incrementation and print data for this frame
                                if frames_this_dcd == 0:

                                        current_step += (total_steps_after_last_nc_len - current_step + steps_per_frame)

                                else:

                                        current_step += steps_per_frame

                                frames_this_dcd += 1

                                # generate print line and append to the end of the data file
                                outline = '%.0f' %simulation_phase+'\t'+'%.0f' %nc_length +'\t'+ '%.9f' %(float(current_step)*tstep*0.001)+'\t'\
                                         +'%.4f' %(float(type_i_frame+type_ii_frame)/float(type_i_ref+type_ii_ref))+'\n'#+'\t'+'%.5f' %chi_frame+'\n'

                                with open (props_out_path, "a") as props:
                                        props.write(outline)

                else:

                        current_step += (total_steps_after_last_nc_len - current_step + steps_per_frame)

                total_steps_after_last_nc_len += steps_this_nc_len

                # append final time point from synthesis phase
                if nc_length == int(max_res):

                        print '\n'+'Appending final time point... '

                        # grab the final coordinate file for this protein and trajectory
                        final_coor      = os.popen('ls -1v '+cot_data_dir+str(traj)+'_r'+str(nc_length)+'_steps*cor').read().split()[0]

                        # as well as the final psf file
                        final_psf       = final_coor.split('cor')[0]+'psf'

                        # check to see if a binary to ascii conversion is necessary
                        if binary == '1' and already_looped == 0:

                                # if so, carry out the conversion (place the new PSF in analysis_dir)
                                conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+final_psf+\
                                           ' outpsf='+analysis_dir+final_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                                # run the command as a child process
                                subprocess.call(conv_cmd, shell=True)

                                # get the name of the new PSF file
                                final_psf = analysis_dir+final_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                # carry out the cor file conversion as well
                                conv_cmd = charmm_exec+' < '+cor_conv_inp+' psf='+final_psf+' cor='+final_coor+\
                                           ' outcor='+analysis_dir+final_coor.split('/')[-1].split('.cor')[0]+'_ascii.cor >> junk.out'

                                # run the conversion command
                                subprocess.call(conv_cmd, shell=True)

                                # repath final_coor to direct to the ASCII version
                                final_coor = analysis_dir+final_coor.split('/')[-1].split('.cor')[0]+'_ascii.cor'

                                print 'ASCII psf '+final_psf+' generated.'
                                print 'ASCII cor '+final_coor+' generated.'

                        # if we do not need to make the binary to ascii conversion
                        elif binary == '1' and already_looped == 1:

                                # get the name of the new PSF file
                                final_psf = analysis_dir+final_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                # as well as the name of the new coor file
                                final_coor = analysis_dir+final_coor.split('/')[-1].split('.cor')[0]+'_ascii.cor'

                        # otherwise
                        else:

                                # do nothing
                                pass

                        # load the universe for this final time point
                        final_univ      = Universe(final_psf, final_coor, format='CRD')

                        # make an AtomGroup for this domain based on the domain definitions
                        final_atomgroup = final_univ.select_atoms('segid A and ('+domain_defs['dom'+str(dom)]+')')

                        # dcalculate type i and ii contacts for this final structure
                        type_i_frame, type_ii_frame = frame_type_i_and_ii_contacts(local_cmap, dmap, final_atomgroup.positions, sec_elems)

                        # calculate pairwise distances for this frame
                        # note that temp_frame_chi_dmap is not necessary because the ref_dmap is the same size as the frame_dmap
                        # at the final nascent chain length
                        #frame_chi_dmap, frame_chi_dists = structural_overlap_dmap(final_atomgroup.positions, domain_defs['dom'+str(dom)],\
                        #                                  len(domain_ref_atomgroup.positions))

                        # compute the structural overlap function for this domain
                        #chi_frame = structural_overlap_function(ref_chi_dmap, frame_chi_dmap, ref_chi_dists)

                        # the time at the end of the simulation has already been calculated as total_steps_after_last_nc_len, so we can just
                        # go ahead and calculate Q and print the final time point
                        outline = '%.0f' %simulation_phase+'\t'+'%.0f' %nc_length +'\t'+ '%.9f' %(float(total_steps_after_last_nc_len)*tstep*0.001)+'\t'+\
                                  '%.4f' %(float(type_i_frame+type_ii_frame)/float(type_i_ref+type_ii_ref))+'\n'#'\t'+'%.5f' %chi_frame+'\n'

                        with open (props_out_path, "a") as props:
                                props.write(outline)

                print 'Calculation for domain '+str(dom)+' nascent chain length '+str(nc_length)+' done; run time at completion:',
                print  datetime.now() - startTime,'\n'

        # initialize counter for the total time
        current_step = total_steps_after_last_nc_len

################################################################

        simulation_phase = 2

        # only perform termination Q calculation if NOT xstal data
        if xstal == False:

                print 'Beginning translation termination calculation for domain '+str(dom)+'\n'

                # perform translation termination Q analysis
                # need to increment the time by an extra frame worth of steps due to binary conversion issue
                # can use the same psf for termination Q analysis as was used for the final nascent chain length
                # calculation        

                # get the dcd file path
                dcd_file = os.popen('ls -1v '+term_data_dir+traj+'_termi_merged_*.dcd').read().split()[0]

                print 'Domain '+str(dom)+' contacts being calculated for: ', dcd_file, '\n'

                # account for off-by-one error due to binary psf's that caused the 
                # first DCD frame not to be written to file
                if binary == '1':

                        frames_this_dcd = 2
                
                        current_step += 10000

                else:

                        frames_this_dcd  = 1

                        current_step += 5000

                # load the current dcd file into current_univ
                current_univ     = Universe(final_psf, dcd_file)

                # make an AtomGroup for this domain based on the domain definitions
                # in the .cntrl file; remember that the ribosome is present here, so we need to specifically select Segid A
                domain_atomgroup = current_univ.select_atoms('segid A and ('+domain_defs['dom'+str(dom)]+')')

                # loop over frames
                for ts in current_univ.trajectory:

                        # calculate the Type I and II contacts for this frame
                        type_i_frame, type_ii_frame = frame_type_i_and_ii_contacts(cmap, dmap, domain_atomgroup.positions, sec_elems)

                        # calculate the chi dmap for this frame
                        #frame_chi_dmap, frame_chi_dists = structural_overlap_dmap(domain_atomgroup.positions,domain_defs['dom'+str(dom)],\
                        #                                  len(domain_ref_atomgroup.positions))

                        # compute the structural overlap function for this domain
                        #chi_frame = structural_overlap_function(ref_chi_dmap, frame_chi_dmap, ref_chi_dists)

                        # the time at the end of the simulation has already been calculated as total_steps_after_last_nc_len, so we can just
                        # go ahead and calculate Q and print the final time point
                        outline = '%.0f' %simulation_phase+'\t'+'%.0f' %max_res +'\t'+ '%.9f' %(float(current_step)*tstep*0.001)+'\t'+\
                                  '%.4f' %(float(type_i_frame+type_ii_frame)/float(type_i_ref+type_ii_ref))+'\n'#'\t'+'%.5f' %chi_frame+'\n'

                        with open (props_out_path, "a") as props:
                                props.write(outline)

                        current_step += steps_per_frame

############################################################################

        print 'Beginning post-translational calculation for domain '+str(dom)+'\n'

        simulation_phase = 3

        if xstal == True:

                traj_psf = ref_psf

        else:

                traj_psf = os.popen('ls -1v '+pt_data_dir+traj+'_*trans_term_*psf' ).read().split()[0]

        if binary == '1' and already_looped == 0 and xstal == False:

                # carry out the conversion (place the new PSF in analysis_dir)
                conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+traj_psf+\
                           ' outpsf='+analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                # run the command as a child process
                subprocess.call(conv_cmd, shell=True)

                # get the name of the new PSF file
                traj_psf = analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

        elif binary == '1' and already_looped == 1 and xstal == False:

                traj_psf = analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

        else:

                pass

        current_step += 5000.

        # make a list of the DCD files to analyze for this trajectory
        if xstal == True:
        
                dcd_file_list        = os.popen('ls -1v '+pt_data_dir+traj+'_*'+'xstal'+'*dcd').read().split()

        else:

                dcd_file_list        = os.popen('ls -1v '+pt_data_dir+traj+'_*'+'post_trans'+'*dcd').read().split()

        # loop over the set of dcd files for this trajectory and job type        
        for dcd_file in dcd_file_list:

                print 'Domain '+str(dom)+' contacts being calculated for: ', dcd_file

                frames_this_dcd  = 1

                # load the current dcd file into current_univ
                current_univ     = Universe(traj_psf, dcd_file)

                # make an AtomGroup for this domain based on the domain definitions
                # in the .cntrl file
                domain_atomgroup = current_univ.select_atoms(domain_defs['dom'+str(dom)])

                # loop over frames
                for ts in current_univ.trajectory:

                        # calculate the Type I and II contacts for this frame
                        type_i_frame, type_ii_frame = frame_type_i_and_ii_contacts(cmap, dmap, domain_atomgroup.positions, sec_elems)

                        # calculate the chi dmap for this frame
                        #frame_chi_dmap, frame_chi_dists = structural_overlap_dmap(domain_atomgroup.positions,domain_defs['dom'+str(dom)],\
                        #                                  len(domain_ref_atomgroup.positions))

                        # compute the structural overlap function for this domain
                        #chi_frame = structural_overlap_function(ref_chi_dmap, frame_chi_dmap, ref_chi_dists)

                        outline = '%.0f' %simulation_phase+'\t'+'%.0f' %max_res +'\t'+ '%.9f' %(float(current_step)*tstep*0.001)+'\t'+\
                                  '%.4f' %(float(type_i_frame+type_ii_frame)/float(type_i_ref+type_ii_ref))+'\n'#'\t'+'%.5f' %chi_frame+'\n'

                        with open (props_out_path, "a") as props:
                                props.write(outline)

                        current_step += steps_per_frame

                print 'Domain '+str(dom)+' calculation done for '+dcd_file+'; runtime at completion:',datetime.now() - startTime 

        # we have already made all of the binary to ascii conversions necessary for this protein, so
        # indicate that by changing the value of already_looped from 0 to 1
        already_looped = 1

#############################################################

# check to see if there are any interfaces we need to analyze
if skip_interface == '1':

        print 'No interface analysis requested, skipping to end of the script.'+'\n'

# or, if there was only one domain, then there is no interface analysis possible
elif num_dom == 1:

        print 'No interface analysis requested, skipping to end of the script.'+'\n'

# if there are then start doing it
else:

        print '\n'+'Beginning interface analysis.'

        # reset the simulation phase to 1, as this is for co-translational simulations again
        simulation_phase = 1

                     #
                    # #   
              ######
        ################# ## Calculate co-translational fraction of native contacts for each interface
                   # 
                    # #   
                     #  

        # make an AtomGroup for the entire protein by selecting all of segid A
        protein_ref_atomgroup = ref_univ.select_atoms('segid A')

        # get distance map, contact map, and total contacts for reference state
        # note that we can use the old function for this calculation - interface analysis
        # was not influence by the bug in the algorithm
        dmap, cmap, total_contacts = get_maps(protein_ref_atomgroup.positions)

        # loop over all interfaces for this protein
        for interface in interface_inp_list:

                # load the interface structure file
                temp_interface_array = read_file(interface)

                # counter for number of lines
                count = 0

                # determine how many interface definition segments there are
                for line in temp_interface_array:

                        count += 1

                # make an array of the correct size
                interface_array = np.zeros((count, 3))

                # reset line counter
                count = 0

                # populate interface_array
                for line in temp_interface_array:

                        line = line.split()

                        # populate interface_array
                        interface_array[count, 0] = line[0]
                        interface_array[count, 1] = line[1]
                        interface_array[count, 2] = line[2]

                        count += 1

                # compute reference state number of contacts
                interface_cont_ref = ref_interface_contacts(cmap, dmap, interface_array)

                # compute reference state chi dmap and count number of pairwise distances formed in native state
                #int_ref_chi_dmap, int_ref_contacts_formed = structural_overlap_intf_dmap(protein_ref_atomgroup.positions, interface_array, \
                #                                            len(protein_ref_atomgroup.positions))

                # if there are no contacts at this interface 
                if interface_cont_ref == 0:

                        # print a message and proceed to the next interface
                        print 'No contacts found for '+interface

                        continue

                current_step                     = steps_per_frame

                total_steps_after_last_nc_len    = 0

                # generate a path for output for this protein, trajectory, and domain
                props_out_path                   = analysis_dir+traj+'_'+interface.split('/')[-1].split('.')[0]+'_Qi.out'

                if xstal == True:

                        props_out_path = analysis_dir+traj+'_'+interface.split('/')[-1].split('.')[0]+'_xstal_Qi.out'

                # check to see if this path exists already and die if it does
                if os.path.exists(props_out_path) == True:

                        print 'The procedurally generated output file for this protein and trajectory of '+props_out_path+' already exists.'
                        print 'Manually delete it to continue.'

                        sys.exit()

                # print a time point for time zero; needed for binning and thus not needed for XSTAL calculations
                #if xstal == False:

                #        outline = '%.0f' %simulation_phase+'\t'+'%.0f' %float(1.0)+'\t'+'%.9f' %float(0.0)+'\t'+'%.4f' %float(0.0)+\
                #                  '\t' +'%.5f' %float(0.0)+'\n'

                #        with open (props_out_path, "a") as props:
                #                props.write(outline)

                # loop over the set of dcd files for this trajectory and job type        
                for nc_length in range (init_res, max_res+1):

                        # do not perform this analysis if xstal == True, no synthesis DCDs to analyze
                        if xstal == True:

                                break

                        print 'Starting calculation for interface '+interface.split('/')[-1].split('.')[0]+' nascent chain length '+str(nc_length)+' ... '

                        # get the cmap and dmap for only the contacts possible at this nascent chain length
                        local_cmap            = np.zeros((len(protein_ref_atomgroup.positions), len(protein_ref_atomgroup.positions)))

                        # only populate local_cmap to the index corresponding to nc_length
                        for a in range (0, nc_length):

                                for b in range (0, nc_length):

                                        local_cmap[a, b] = cmap[a, b]
                                        local_cmap[b, a] = cmap[b, a]

                        # generate list of the dcd files for this nascent chain length and trajectory
                        dcd_file_list         = os.popen('ls -1v '+cot_data_dir+str(traj)+'_r'+str(nc_length)+'_steps*dcd').read().split()

                        # make sure there is exactly 1 DCD file for this trajectory and nc length
                        if len(dcd_file_list) != 1:

                                print 'ERROR IN DCD FILE SELECTION.'
                                print 'CHECK TRAJ '+str(traj)+' NC LENGTH '+str(nc_length)

                                sys.exit()

                        # grab single element from dcd_file_list
                        dcd_file              = dcd_file_list[0]

                        # determine the number of frames we expect to find in this dcd file
                        temp                  = dcd_file.split('/')[-1]

                        steps_this_nc_len     = int(temp.split('_')[2].split('steps')[1])

                        # check to see if this DCD file is empty
                        file_size_zero        = os.stat(dcd_file).st_size == 0

                        # counter for number of frames in this dcd file
                        frames_this_dcd       = 0

                        # if the DCD file contains data
                        if file_size_zero != True:

                                # generate a path to a PSF file
                                psf_file = dcd_file.split('dcd')[0]+'psf'

                                # if we need ascii psfs and we skipped the domain section such that
                                # the required psfs have not been created yet
                                if binary == '1' and skip_domain == '1':

                                        # if so, carry out the conversion (place the new PSF in analysis_dir)
                                        conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+psf_file+\
                                                   ' outpsf='+analysis_dir+psf_file.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                                        # run the command as a child process
                                        subprocess.call(conv_cmd, shell=True)

                                        # get the name of the new PSF file
                                        psf_file = analysis_dir+psf_file.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                        print 'ASCII psf '+psf_file+' generated.'

                                # if we need binary psfs but they have already been generated
                                elif binary == '1' and skip_domain == '0':

                                        # just update the psf path to direct to the ascii version
                                        psf_file = analysis_dir+psf_file.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                # in all other cases
                                else:

                                        # do nothing
                                        pass

                                # load DCD file into PSF to create Universe
                                current_univ = Universe(psf_file, dcd_file)

                                # select all SEGID A atoms
                                protein_atomgroup = current_univ.select_atoms('segid A')

                                # loop over frames in the trajectory file
                                for ts in current_univ.trajectory:

                                        # get total interface contacts in this frame
                                        interface_cont_i = frame_interface_contacts(local_cmap, dmap, protein_atomgroup.positions, interface_array)

                                        # facilitate interface chi analysis by reshaping the current coordinates so they are in the top left
                                        # corner of a matrix with (max_res x max_res) dimensions otherwise populated with zeros
                                        temp_coor = np.zeros((len(protein_ref_atomgroup.positions), 3))

                                        for i in range (0, len(protein_atomgroup.positions)):

                                                # note well, these are Cartesian coordinates and NOT a square matrix
                                                # needs to have 3 columns and N rows
                                                temp_coor[i,0] = protein_atomgroup.positions[i,0]
                                                temp_coor[i,1] = protein_atomgroup.positions[i,1]
                                                temp_coor[i,2] = protein_atomgroup.positions[i,2]

                                        # calculate the chi dmap for this frame
                                        #int_frame_chi_dmap, int_frame_contacts_formed = structural_overlap_intf_dmap(temp_coor, interface_array, nc_length)

                                        # compute the value of chi by comparing the current and reference distance maps
                                        #intf_chi_frame = structural_overlap_function(int_ref_chi_dmap, int_frame_chi_dmap, int_ref_contacts_formed)

                                        # manage time incrementation and print data for this frame
                                        if frames_this_dcd == 0:

                                                current_step += (total_steps_after_last_nc_len - current_step + steps_per_frame)

                                        else:

                                                current_step += steps_per_frame

                                        frames_this_dcd += 1

                                        outline = '%.0f' %simulation_phase+'\t'+'%.0f' %nc_length +'\t'+ '%.9f' %(float(current_step)*tstep*0.001)+'\t'+'%.4f'\
                                                  %(float(interface_cont_i)/float(interface_cont_ref))+'\n'#'\t'+'%.5f' %intf_chi_frame+'\n'

                                        with open (props_out_path, "a") as props:
                                                props.write(outline)

                        else:

                                current_step += (total_steps_after_last_nc_len - current_step + steps_per_frame)

                        total_steps_after_last_nc_len += steps_this_nc_len

                        # append final synthesis time point
                        if nc_length == int(max_res):

                                print 'Appending final time point... '

                                # grab the final coordinate file for this protein and trajectory
                                final_coor      = os.popen('ls -1v '+cot_data_dir+str(traj)+'_r'+str(nc_length)+'_steps*cor').read().split()[0]

                                # as well as the final psf file
                                final_psf       = final_coor.split('cor')[0]+'psf'

                                # if we need ascii psfs and we skipped the domain section such that
                                # the required psfs have not been created yet
                                if binary == '1' and skip_domain == '1':

                                        # if so, carry out the conversion (place the new PSF in analysis_dir)
                                        conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+final_psf+\
                                                   ' outpsf='+analysis_dir+final_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                                        # run the command as a child process
                                        subprocess.call(conv_cmd, shell=True)

                                        # get the name of the new PSF file
                                        final_psf = analysis_dir+final_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                        # carry out the cor file conversion as well
                                        conv_cmd = charmm_exec+' < '+cor_conv_inp+' psf='+final_psf+' cor='+final_coor+\
                                                   ' outcor='+analysis_dir+final_coor.split('/')[-1].split('.cor')[0]+'_ascii.cor >> junk.out'

                                        # run the conversion command
                                        subprocess.call(conv_cmd, shell=True)

                                        # repath final_coor to direct to the ASCII version
                                        final_coor = analysis_dir+final_coor.split('/')[-1].split('.cor')[0]+'_ascii.cor'

                                        print 'ASCII psf '+final_psf+' generated.'
                                        print 'ASCII cor '+final_coor+' generated.'

                                # if the conversion has already been made
                                elif binary == '1' and skip_domain == '0':

                                        # just update the psf path to direct to the ascii version
                                        final_psf = analysis_dir+final_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                                        # as well as the cor path 
                                        final_coor = analysis_dir+final_coor.split('/')[-1].split('.cor')[0]+'_ascii.cor'

                                # in all other cases
                                else:

                                        # do nothing
                                        pass

                                # load the universe for this final time point
                                final_univ      = Universe(final_psf, final_coor, format='CRD')

                                # make an AtomGroup
                                final_atomgroup = final_univ.select_atoms('segid A')

                                # get contacts at this final frame
                                interface_cont_i = frame_interface_contacts(local_cmap, dmap, final_atomgroup.positions, interface_array)

                                # calculate the chi dmap for this frame
                                #int_frame_chi_dmap, int_frame_contacts_formed = structural_overlap_intf_dmap(final_atomgroup.positions, interface_array, max_res)

                                # compute the value of chi by comparing the current and reference distance maps
                                #intf_chi_frame = structural_overlap_function(int_ref_chi_dmap, int_frame_chi_dmap, int_ref_contacts_formed)

                                # the time at the end of the simulation has already been calculated as total_steps_after_last_nc_len, so we can just
                                # go ahead and calculate Q and print the final time point
                                outline = '%.0f' %simulation_phase+'\t'+'%.0f' %nc_length +'\t'+ '%.9f' %(float(total_steps_after_last_nc_len)*tstep*0.001)+'\t'+\
                                          '%.4f' %(float(interface_cont_i)/float(interface_cont_ref))+'\n' #'\t'+'%.5f' %intf_chi_frame+'\n'

                                with open (props_out_path, "a") as props:
                                        props.write(outline)

                        print 'Calculation for '+interface.split('/')[-1].split('.')[0]+' interface nascent chain length '+str(nc_length)+' done; run '+\
                                               'time at completion:', datetime.now() - startTime, '\n'

##################################################################

                simulation_phase = 2

                # initialize counter for the total time
                current_step = total_steps_after_last_nc_len

                if xstal == False:

                        print 'Beginning translation termination calculation for interface '+interface.split('/')[-1].split('.')[0]

                        # perform translation termination Q analysis
                        # need to increment the time by an extra frame worth of steps due to binary conversion issue
                        # can use the same psf for termination Q analysis as was used for the final nascent chain length
                        # calculation        

                        # get the dcd file path
                        dcd_file = os.popen('ls -1v '+term_data_dir+traj+'_termi_merged_*.dcd').read().split()[0]

                        print 'Interface '+interface.split('/')[-1].split('.')[0]+' contacts being calculated for: ', dcd_file, '\n'

                        # account for off-by-one error due to binary psf's that caused the 
                        # first DCD frame not to be written to file
                        if binary == '1':

                                frames_this_dcd = 2

                                current_step += 10000

                        else:

                                frames_this_dcd  = 1

                                current_step += 5000

                        # load the current dcd file into current_univ
                        current_univ     = Universe(final_psf, dcd_file)

                        # make an AtomGroup for the entire protein
                        protein_atomgroup = current_univ.select_atoms('segid A')

                        # loop over frames
                        for ts in current_univ.trajectory:

                                # calculate interface contacts for this frame
                                # get total interface contacts in this frame
                                interface_cont_i = frame_interface_contacts(cmap, dmap, protein_atomgroup.positions, interface_array)

                                # calculate the chi dmap for this frame
                                #int_frame_chi_dmap, int_frame_contacts_formed = structural_overlap_intf_dmap(protein_atomgroup.positions, interface_array, max_res)

                                # compute the value of chi by comparing the current and reference distance maps
                                #intf_chi_frame = structural_overlap_function(int_ref_chi_dmap, int_frame_chi_dmap, int_ref_contacts_formed)

                                outline = '%.0f' %simulation_phase+'\t'+'%.0f' %nc_length +'\t'+ '%.9f' %(float(current_step)*tstep*0.001)+'\t'+'%.4f'\
                                        %(float(interface_cont_i)/float(interface_cont_ref))+'\n'#'\t'+'%.5f' %intf_chi_frame+'\n'

                                with open (props_out_path, "a") as props:
                                        props.write(outline)

                                current_step += steps_per_frame

                # perform post-translational Q analysis
                # need to use a different psf for this stage of the analysis

#############################################################

                simulation_phase = 3

                if xstal == True:

                        traj_psf = ref_psf

                else:

                        traj_psf = os.popen('ls -1v '+pt_data_dir+traj+'_*trans_term_*psf' ).read().split()[0]

                if binary == '1' and already_looped == 0 and xstal == False:

                        # carry out the conversion (place the new PSF in analysis_dir)
                        conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+traj_psf+\
                                   ' outpsf='+analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf >> junk.out'

                        # run the command as a child process
                        subprocess.call(conv_cmd, shell=True)

                        # get the name of the new PSF file
                        traj_psf = analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                elif binary == '1' and already_looped == 1 and xstal == False:

                        traj_psf = analysis_dir+traj_psf.split('/')[-1].split('.psf')[0]+'_ascii.psf'

                else:

                        pass

                current_step += 5000.

                # make a list of the DCD files to analyze for this trajectory
                if xstal == True:

                        dcd_file_list        = os.popen('ls -1v '+pt_data_dir+traj+'_*'+'xstal'+'*dcd').read().split()

                else:

                        dcd_file_list        = os.popen('ls -1v '+pt_data_dir+traj+'_*'+'post_trans'+'*dcd').read().split()

                # loop over the set of dcd files for this trajectory and job type        
                for dcd_file in dcd_file_list:

                        frames_this_dcd  = 1

                        # load the current dcd file into current_univ
                        current_univ     = Universe(traj_psf, dcd_file)

                        # make an AtomGroup for this domain based on the domain definitions
                        # in the .cntrl file
                        protein_atomgroup = current_univ.select_atoms('segid A')

                        # loop over frames
                        for ts in current_univ.trajectory:

                                # calculate interface contacts for this frame
                                # get total interface contacts in this frame
                                interface_cont_i = frame_interface_contacts(cmap, dmap, protein_atomgroup.positions, interface_array)

                                # calculate the chi dmap for this frame
                                #int_frame_chi_dmap, int_frame_contacts_formed = structural_overlap_intf_dmap(protein_atomgroup.positions, interface_array, max_res)

                                # compute the value of chi by comparing the current and reference distance maps
                                #intf_chi_frame = structural_overlap_function(int_ref_chi_dmap, int_frame_chi_dmap, int_ref_contacts_formed)

                                outline = '%.0f' %simulation_phase+'\t'+'%.0f' %nc_length +'\t'+ '%.9f' %(float(current_step)*tstep*0.001)+'\t'+'%.4f'\
                                        %(float(interface_cont_i)/float(interface_cont_ref))+'\n'#'\t'+'%.5f' %intf_chi_frame+'\n'

                                with open (props_out_path, "a") as props:
                                        props.write(outline)

                                current_step += steps_per_frame

                        print 'Interface '+interface.split('/')[-1].split('.')[0]+' calculation done for '+dcd_file+'; runtime at completion:',datetime.now() - startTime

# do some cleanup - if temporary ascii psfs were created
if binary == '1' and xstal == False:

        # delete them now for this trajectory
        subprocess.call('rm '+analysis_dir+traj+'_*ascii.*', shell=True)

# finish the calculations and be done with it
print 'All requested calculations done; required time:', datetime.now() - startTime

print '\n'+'NORMAL TERMINATION'

# prepare output line
if xstal != False:

        print_line = 'XS TRAJ '+traj+' DONE.'+'\n'

else:

        print_line = 'PT TRAJ '+traj+' DONE.'+'\n'

# print to file a statement that this trajectory completed normally
with open (analysis_dir+'done.log', "a") as tempf:

        tempf.write(print_line)
