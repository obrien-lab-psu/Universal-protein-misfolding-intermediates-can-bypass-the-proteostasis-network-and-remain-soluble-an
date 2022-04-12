# load modules
import MDAnalysis as MDA
import sys, os
import numpy as np
import subprocess

# Original version: 11-August-2020
# Current version{ 8-September-2020

# GOAL: write a function that, given a set of secondary structures, a DCD file,
#       and a nascent chain length between 1 and the maximum length of the protein,
#       will print the fraction of native contacts

# getQ(ref_cor ref_psf, max_length, current_psf, current_length, ssinp, cor_file=current_file)
usage = 'python python_scripts/overall_Q.py [ref .cor] [ref .psf] [max length] [current .psf] [current length] [sec. struc. inp] [trajectory file] [.cor? True or False] [outfile path] [binary? True or False] [multidom? True or False]'

#print sys.argv
#print len(sys.argv)
if len(sys.argv) != 12:
        print usage
        sys.exit()

ref_cor = sys.argv[1]
ref_psf = sys.argv[2]
max_length = int(sys.argv[3])
current_psf = sys.argv[4]
current_length = int(sys.argv[5])
ssinp = sys.argv[6]
current_file = sys.argv[7]
temp = sys.argv[8]
if temp == 'True':
        coor = True
else:
        coor = False
if current_length > max_length:
        print ('current_length of', current_length, 'is impossible with a total length of', max_length)
        sys.exit()
outfile = sys.argv[9]
if os.path.exists(outfile):
        print (outfile, 'already exists - it will be overwritten!')
        #sys.exit()
temp = sys.argv[10]
if temp == 'True':
        binary = True
else:
        binary = False
temp = sys.argv[11]
if temp == 'True':
        multidom = True
else:
        multidom = False

def getQ(ref_cor_file, ref_psf_file, max_length, psf_file, length, ssinp, outfile, dcd_file=False, cor_file=False):
    
    # ref_cor_file    : str, path to reference state coordinate file in CHARMM format
    # ref_psf_file    : str, path to reference state protein structure file in CHARMM format
    # max_length      : int, the number of amino acid residues in the full-length protein
    # psf_file        : str, path to protein structure file into which dcd_file can be loaded
    # length          : int, the current nascent chain length
    # ssinp           : str, path to the secondary structure input file for the native state of this protein
    # outfile         : str, path to which output time series will be written
    # dcd_file        : (optional, default=False), path to DCD file to be analyzed or False if cor file is to be analyzed instead
    # cor_file        : (optional, default=False), path to COR file to be analyzed or False if dcd file is to be analyzed instead
    
    threshold = 8.0  # threshold distance in Angstroms for determining whether a contact is formed
    adjust = 1.2     # multiplicative adjustment factor allowing for slightly longer distances to still be
                     # counted as contacts due to thermal fluctuations at finite temperature
    ref_contacts = 0 # counter for the number of contacts formed in the reference state
    
    # make sure all of the input files exist and run a few quick checks for input parameter correctness
    if dcd_file != False and cor_file != False:
        print ('You can enter either a dcd file or a cor file for analysis, not both!')
        print ('DCD:', dcd_file, '\t', 'COR:', cor_file)
    if dcd_file != False:
        if os.path.exists(dcd_file) != True:
            print (dcd_file, 'does not exist.')
            sys.exit()
    if cor_file != False:
        if os.path.exists(cor_file) != True:
            print (cor_file, 'does not exist.')
            sys.exit()
    if os.path.exists(ref_cor_file) != True:
        print (ref_cor_file, 'does not exist.')
        sys.exit()
    if os.path.exists(ref_psf_file) != True:
        print (ref_psf_file, 'does not exit.')
        sys.exit()
    if os.path.exists(ssinp) != True:
        print (ssinp, 'does not exist.')
        sys.exit()
    if length > max_length:
        print ('The entered current length of', length, 'is greater than the entered maximum length of', max_length)
        sys.exit()
    
    if '1sg5' in ref_cor_file:
        shift = -2
    elif '1rqj' in ref_cor_file:
        shift = -21
    else:
        shift = 0 
    ss_residues = getSSresidues(ssinp, shift)
    #print ss_residues
    #sys.exit()
    # Note: I do not think that I need to bother with identifying the contiguous
    #       sections of secondary structure (i.e., individual helices and strands)
    #       because I just want to compute the Q over all of them. 
    # 8-Sept-2020: this is incorrect because it results in counting residues that were
    # not counted in the original Q analysis. See, for example, 1a6j stride output. There is a 310 Helix from
    # 16-18 which is itself too small to count (3 residues) but directly abuts a strand from 19-23. The code as
    # previously written just saw a structure from 16 to 23, which it counted, when my previous analysis would
    # only use 19-23 and not 16-18. It is not clear to me which is more correct, but I have adjusted this code
    # so that it precisely matches the way in which this issue is handled in the original program
    
    # set up Universe object for reference system
    ref_universe = MDA.Universe(ref_psf_file, ref_cor_file, format='crd')
    ref_coor     = ref_universe.select_atoms('segid A').positions
    
    # generate reference contact and distance maps
    # Note: we only need to make these as large as the current system being analyzed,
    #       and in fact we need to be careful to only count contacts that are physically
    #       possible at the current nascent chain length
    
    # these maps are references starting from zero, but the secondary structure definitions start from 1;
    # be careful to apply a minus 1 adjustment to the residue indices in ss_residues before entering data
    # into these arrays to preserve correct indexing
    cmap = np.zeros((max_length, max_length))
    dmap = np.zeros((max_length, max_length))
    
    # we only need to check contacts between residues that are both in ss_residues
    # and also have at least i -> i+3 difference in the primary sequence
    #temp = open('test.dat', 'w')
    for i in range (0, len(ss_residues)):
        for j in range (i+1, len(ss_residues)):
            #temp.write(str(ss_residues[i]) + '\t' + str(ss_residues[j]) + '\n')
            if ss_residues[j] - ss_residues[i] > 3: # enforce primary sequence difference criterion
                if ss_residues[j] <= max_length and ss_residues[i] <= max_length:
                    dist = d(ref_coor[ss_residues[j]-1, :], ref_coor[ss_residues[i]-1, :])
                    if dist <= threshold:
                        cmap[ss_residues[i]-1, ss_residues[j]-1] = 1
                        cmap[ss_residues[j]-1, ss_residues[i]-1] = 1
                        dmap[ss_residues[i]-1, ss_residues[j]-1] = dist
                        dmap[ss_residues[j]-1, ss_residues[i]-1] = dist
                        #print 'CONTACT:', ss_residues[i], ss_residues[j]
                        ref_contacts += 1
                        #temp.write('Contact: ' + '%.0f' %ss_residues[i] + '\t' + '%.0f' %ss_residues[j] + '\t' +/
                        #            '%.0f' %cmap[ss_residues[i]-1, ss_residues[j]-1] + '\t' + '%.5f' %dmap[ss_residues[i]-1, ss_residues[j]-1] + '\n')
    #print ref_contacts
    #print 'Number of reference contacts:', ref_contacts
    #sys.exit()
    # Now, determine which contacts are currently formed in each frame of the current DCD file
    ofile = open(outfile, 'w')
    if dcd_file != False and cor_file == False:
    
        # set up Universe for current system to be analyzed
        universe     = MDA.Universe(psf_file, dcd_file)
        selection    = universe.select_atoms('segid A')
        f = 1
        for frame in universe.trajectory:
            coor = selection.positions
            contacts = 0
            for i in range (0, len(ss_residues)):
                for j in range (i+1, len(ss_residues)):
                    if ss_residues[j] - ss_residues[i] > 3:
                        if ss_residues[j] <= length and ss_residues[i] <= length:
                            dist = d(coor[ss_residues[j]-1, :], coor[ss_residues[i]-1, :])
                            if dist <= dmap[ss_residues[i]-1, ss_residues[j]-1]*adjust and cmap[ss_residues[i]-1, ss_residues[j]-1] == 1:
                                contacts += 1
                                #temp.write(str(frame)+'\t'+'%.5f' %(float(contacts)/float(ref_contacts))+'\n')
            #print contacts, ref_contacts
            #print (str(f)+'\t'+'%.5f' %(float(contacts)/float(ref_contacts)))
            ofile.write(str(f)+'\t'+'%.5f' %(float(contacts)/float(ref_contacts))+'\n')
            f += 1
    
    elif dcd_file == False and cor_file != False:
        # set up Universe for current system to be analyzed
        universe     = MDA.Universe(psf_file, cor_file, format='CRD')
        selection    = universe.select_atoms('segid A')
        f = 1
        coor = selection.positions
        contacts = 0
        for i in range (0, len(ss_residues)):
            for j in range (i+1, len(ss_residues)):
                if ss_residues[j] - ss_residues[i] > 3:
                    if ss_residues[j] <= length and ss_residues[i] <= length:
                        dist = d(coor[ss_residues[j]-1, :], coor[ss_residues[i]-1, :])
                        if dist <= dmap[ss_residues[i]-1, ss_residues[j]-1]*adjust and cmap[ss_residues[i]-1, ss_residues[j]-1] == 1:
                            contacts += 1
                            #temp.write(str(frame)+'\t'+'%.5f' %(float(contacts)/float(ref_contacts))+'\n')
        #print (str(f)+'\t'+'%.5f' %(float(contacts)/float(ref_contacts)))
        ofile.write(str(f)+'\t'+'%.5f' %(float(contacts)/float(ref_contacts))+'\n')
        
    else:
        print ('How have you achieved this?')
        sys.exit()
        
    return
    
def d(a, b):
    
    # a: numpy array corresponding to x, y, z coordinates for one atom/interaction site
    # b: numpy array corresponding to x, y, z coordinates for a second atom/interaction site
    
    return np.sqrt(((a[0]-b[0])**2.)+((a[1]-b[1])**2.)+((a[2]-b[2])**2.))

def getSSresidues(ssinp, shift):
        # ssinp: path to stride output file for this protein
        # shift: integer constant to adjust residue indices from LOC entries
        #        in the stride output file so that they index from 1
        # returns: list of residues in secondary structures

        temp = open(ssinp)
        stride = temp.readlines()
        temp.close()

        ss_residues = []

        for line in stride:
                line = line.split()
                if line[0] == 'LOC':
                        if 'Helix' in line[1]:
                                r1, r2 = int(line[3]), int(line[6])
                                if r2 - r1 + 1 >= 4:
                                        for i in range (r1, r2+1):
                                                ss_residues.append(i+shift)
                        elif 'Strand' in line[1]:
                                r1, r2 = int(line[3]), int(line[6])
                                if r2 - r1 + 1 >= 4:
                                        for i in range (r1, r2+1):
                                                ss_residues.append(i+shift)
                        else:
                                pass
        ss_residues.sort()
        return ss_residues

if binary:
        # need to make an ASCII copy of the psf in ascii_psfs
        psf_conv_inp = 'cont_syn/charmm_scripts/simple_binary_to_ascii_psf.inp'
        if multidom:
                charmm_exec = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm_multidom'
        else:
                charmm_exec = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm'
        outpsf = 'ascii_psfs/'+current_psf.split('/')[-1]
        conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+current_psf+\
                   ' outpsf='+outpsf+' >> junk2.out' 
        subprocess.call(conv_cmd, shell=True)
        current_psf = outpsf
        #print conv_cmd
        # if we are analyzing a coordinate file then we also need to make an ascii copy of the coordinate file
        if coor == True:
                cor_conv_inp = 'cont_syn/charmm_scripts/simple_binary_to_ascii_cor.inp'
                outcor = 'ascii_psfs/'+current_file.split('/')[-1]
                conv_cmd = charmm_exec+' < '+cor_conv_inp+' psf='+current_psf+' cor='+current_file+\
                           ' outcor='+outcor+' >> junk2.out'
                current_file = outcor
                subprocess.call(conv_cmd, shell=True)
                #print conv_cmd

# function to find contiguous sets of residues in a list of residue indices
# from https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
#def residueRanges(nums):
#        nums = sorted(set(nums))
#        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
#        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
#        return list(zip(edges, edges))

if coor == True:
        getQ(ref_cor, ref_psf, max_length, current_psf, current_length, ssinp, outfile, cor_file=current_file)
else:
        getQ(ref_cor, ref_psf, max_length, current_psf, current_length, ssinp, outfile, dcd_file=current_file)

#getQ('1cli_chain_a_rebuilt_mini_ca.cor', '1cli_chain_a_rebuilt_mini_ca.psf', 345,
#     'dcd_files/1_r200_steps1263034_1cli.psf', 200, '1cli_sec_struc.txt', dcd_file='dcd_files/1_r200_steps1263034_1cli.dcd')

#getQ('1cli_chain_a_rebuilt_mini_ca.cor', '1cli_chain_a_rebuilt_mini_ca.psf', 345,
#     'dcd_files/1_r200_steps1263034_1cli.psf', 200, '1cli_sec_struc.txt', cor_file='dcd_files/1_r200_steps1263034_1cli.cor')

#getQ('1cli_chain_a_rebuilt_mini_ca.cor', '1cli_chain_a_rebuilt_mini_ca.psf', 345,
#     'dcd_files/1_r201_steps99717_1cli.psf', 201, '1cli_sec_struc.txt', dcd_file='dcd_files/1_r201_steps99717_1cli.dcd')

#getQ('1cli_chain_a_rebuilt_mini_ca.cor', '1cli_chain_a_rebuilt_mini_ca.psf', 345,
#     'dcd_files/1_r201_steps99717_1cli.psf', 201, '1cli_sec_struc.txt', cor_file='dcd_files/1_r201_steps99717_1cli.cor')

"""
Next steps in implementation

(1) Need to write a program to generate commands for this function and then run them all in parallel
(2) The output file names will need to record the number of steps
(3) Another program can then be used to merge the individual files into one time series taking into account
    all of the issues with concatenating the time series
"""
