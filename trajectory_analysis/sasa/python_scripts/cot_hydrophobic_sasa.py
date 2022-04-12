import os, sys
import numpy as np
import MDAnalysis as MDA
import numpy as np
import subprocess
# PURPOSE: to compute the co-translational SASA
# for xstal trajectories this is not a true co-translational SASA, it is the SASAS of hydrophobic residues
# in the chain at different lengths

# Updated 21-September-2020 to handle issues arising from xstal_comp values equal to zero, which occurs when
# analyzing trajectories initiated from nascent chain length 1, as the reference folded trajs will have 0.0 SASA for
# all nascent chain lengths until they emerge from the tunnel. If folded_val == 0.0 and current val also ~0.0, pass to zero

# function to open files and read them into memory
def read_file(f_file):
        f_temp = open(f_file)
        temp   = f_temp.readlines()
        f_temp.close()
        return temp

# this version of id_hydrophob_res uses MDAnalysis rather than custom code to determine amino acid identities - should be more robust
# added 2-September-2020
def id_hydrophob_res(f_pdb):

        # f_pdb -> file path to protein databank file for the protein
        # currently being studied

        if os.path.exists(f_pdb) != True:
                print (f_pdb, 'does not exit')
                sys.exit()

        universe = MDA.Universe(f_pdb)

        # list of all 20 standard amino acids by 3-letter code
        # note HIS -> HSE in all multi-domain PDBs, but some single-
        # domain PDBs that did not require rebuilding in CHARMM have HIS
        amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                       'GLU', 'GLN', 'GLY', 'HSE', 'HIS',
                       'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                       'PRO', 'SER', 'THR', 'TRP', 'TYR',
                       'VAL']

        # initialize list of hydrophobic residue types

        # this list is not based on any hydrophobicity scale, it is based on
        # Ed's general feelings on the issue 
        hydro_list   = ['ILE', 'VAL', 'LEU', 'PHE',
                        'CYS', 'MET', 'ALA', 'GLY', 'TRP']

        # list of hydrophobic residues for this pdb file
        hydro_resids = []

        residues = universe.select_atoms('all').residues.resids
        offset = residues[0] # needed for 1sg5 and 1rqj which are not indexed from 1
        for res in residues:
                aa = universe.select_atoms('resid '+str(res)).resnames[0]
                if aa not in amino_acids:
                        print (aa, 'is not recognized as an amino acid type')
                        sys.exit()
                if aa in hydro_list:
                        hydro_resids.append(res-offset+1)
                        print (res-offset+1, aa, 'hydrophobic')
                else:
                        print (res-offset+1, aa, 'NOT hydrophobic')
                        pass

        # return the list of resids corresponding to hydrophobic residues
        return hydro_resids

def str2bool(s):
        if s == 'True':
                return True
        else:
                return False

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

# usage statement
usage = 'python python_scripts/cot_hydrophobic_sasa.py [cntrl file] [comp; Boolean] [trajectory index]'

if len(sys.argv) != 4:
        print usage
        sys.exit()

cntrl = read_file(sys.argv[1])

pdb = sys.argv[1].split('/')[-1].split('_')[0]

# read the cntrl file and get information about binary to ascii conversion
cntrl = read_file(sys.argv[1])
for line in cntrl:
        line = line.split('=')
        if line[0].strip() == 'binary':
                temp = line[1].strip()
                if temp == '0':
                        binary = False
                elif temp == '1':
                        binary = True
                else:
                        print 'Could not determine whether or not data are in binary'
                        sys.exit()
        elif line[0] == 'psf_conv_inp':
                psf_conv_inp = line[1].strip()
        else:
                pass

# if comp == True, assume we have already computed the <SASA> as a function of nascent chain length for the folded-state simulations
# this means that we can compute the zeta metric for the synthesis trajectories
comp = str2bool(sys.argv[2])

traj = sys.argv[3]

md_proteins = ['1cli','1d2f','1duv','1ef9','1fts',
               '1fui','1ger','1glf','1gqe','1gyt',
               '1gz0','1ksf','1ng9','1p7l','1qf6',
               '1svt','1t4b','1u0b','1uuf','1w78',
               '1xru','1xvi','2fym','2h1f','2hg2',
               '2hnh','2id0','2kfw','2kx9','2ptq',
               '2qcu','2qvr','2r5n','2wiu','2ww4',
               '3brq','3gn5','3m7m','3nxc','3ofo',
               '3pco','3qou','4dcm','4dzd','4e8b',
               '4fzw','4hr7','4im7','4iwx','4kn7', '4g36']

# get path to pdb format file for this protein from atomistic_pdbs/
pdb_file               = os.popen('ls -1v atomistic_pdbs/'+pdb+'_*.pdb').read().split()[0]

# pass PDB through function to ID hydrophobic residues by their resid (note this is for the full chain!)
analysis_resids       = id_hydrophob_res(pdb_file)

# determine directory for dcd files in cont_syn/
dcd_dir = 'singledom'
if pdb in md_proteins:
        dcd_dir = 'multidom'

# get first and last residue for this analysis
init = int(os.popen('grep init_res cntrl_files/'+pdb+'_analysis.cntrl').read().split()[2])
nres = int(os.popen('grep max_res cntrl_files/'+pdb+'_analysis.cntrl').read().split()[2])

# Bukau has found proteins do not interact strongly with TF until ~100 amino acids in length
init = 100

# residues are considered to be sterically accessible to TF if their x-coordinate is >= xref
xref = 100.0
frame_count = 0
steps_per_frame = 5000.
total_steps_after_last_nc_len = 0
current_step = steps_per_frame

# determien output file name
if comp:
        # set up output file and make sure it doesn't already exist
        outpath = 'analysis/'+pdb+'/'+traj+'_cot_hydrophobic_sasa_comp.txt'
        if os.path.exists(outpath):
                print outpath, 'already exists. It will be overwritten!'
                #sys.exit()
else:
        outpath = 'analysis/'+pdb+'/'+traj+'_cot_hydrophobic_sasa_folded.txt'
        if os.path.exists(outpath):
                print outpath, 'already exists. It will be overwritten!'
                #sys.exit()
ofile = open(outpath, 'w')

"""
if comp:
        # add in a zero time point to facilitate binning
        # the initial value should almost certainly be 0%, so I am just assuming that here (will not matter for binning)
        ofile.write('0.000000000'+'\t'+str(init)+'\t'+'0       0.000000000     0.000000000\n')
else:
        ofile.write('0.000000000'+'\t'+str(init)+'\t'+'0       0.000000000\n')
"""

#if pdb == '1glf':
#        map_1glf = {}
#        temp = open('1glf_trajs.txt').readlines()
#        for traj in range (1, 50+1):
#                map_1glf[str(traj)] = temp[traj-1].strip()

# loop over nascent chain length
for r in range (init, nres+1):

        # load coordinate data for this nascent chain length
        if pdb == '4g36':
                dcd = os.popen('ls -1v cont_syn/continuous_synthesis/proteins/4g36/output/'+traj+'_r'+str(r)+'_steps*dcd').read().split()[0]
        else:
                dcd = os.popen('ls -1v cont_syn/'+dcd_dir+'/'+pdb+'/output/'+traj+'_r'+str(r)+'_steps*dcd').read().split()[0]

        # determine whether or not this file is empty
        file_size_zero = os.stat(dcd).st_size == 0

        # determine number of steps run for this nascent chain length
        steps_this_nc_len = int(dcd.split('/')[-1].split('_')[2].strip('steps'))

        if file_size_zero != True:

                psf = dcd.split('.')[0]+'.psf'

                if binary:
                        # need to make an ASCII copy of the psf in ascii_psfs
                        psf_conv_inp = 'cont_syn/charmm_scripts/simple_binary_to_ascii_psf.inp'
                        if pdb in ['1ng9', '1ksf', '1qf6', '2hnh', '2id0', '2r5n', '3pco', '4kn7', '1h16', '4g36']:
                                charmm_exec = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm_multidom'
                        else:
                                charmm_exec = '/storage/home/dan182/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm'
                        outpsf = 'ascii_psfs/'+psf.split('/')[-1]
                        conv_cmd = charmm_exec+' < '+psf_conv_inp+' psf='+psf+\
                                   ' outpsf='+outpsf+' >> junk2.out'
                        print conv_cmd
                        print outpsf
                        print psf
                        subprocess.call(conv_cmd, shell=True)
                        psf = outpsf

                # get Universe object and selection
                u = MDA.Universe(psf, dcd)
                protein = u.select_atoms('segid A')

                # load SASA data for this nascent chain length
                f_asa = os.popen('ls -1v analysis/'+pdb+'/asa/'+traj+'_r'+str(r)+'_steps*').read().split()[0]
                #if pdb == '1glf' and comp == True:
                #        f_asa = os.popen('ls -1v analysis/'+pdb+'/asa/'+map_1glf[traj]+'_r'+str(r)+'_steps*').read().split()[0]
                asa = np.loadtxt(f_asa)

                # get the xstal structure reference value of <A(subset)> for this nascent chain length
                if comp:
                        folded_val = '0'
                        #raw = read_file('1cli_xstal_cot_hydrophobic_sasa.dat')
                        raw= read_file('folded_cot_hydrophobic_sasa.dat')
                        for line in raw:
                                line = line.split()
                                if line[0] == pdb and int(line[1]) == r:
                                        folded_val = np.float64(line[2])
                        if folded_val == '0':
                                print 'Failed to find correct folded_val for', pdb, 'r='+str(r)
                                sys.exit()

                # get the list of hydrophobic residues added to the chain by this nascent chain length
                resids = []
                for resid in analysis_resids:
                        if resid <= r:
                                resids.append(resid)
                print 'hydrophobic r='+str(r), resids
                # make sure the lengths of the trajectory and asa data match
                if len(np.shape(asa)) == 1:
                        if len(u.trajectory) != 1:
                                print 'Length Mismatch: dcd file', dcd, 'has', len(u.trajectory), 'coordinate sets, but the asa file', f_asa, 'has 1 line'
                                sys.exit()
                else:
                        if len(u.trajectory) != len(asa):
                                print 'Length Mismatch: dcd file', dcd, 'has', len(u.trajectory), 'coordinate sets, but the asa file', f_asa, 'has', len(asa),'lines.'
                                sys.exit()
                # loop over the frames in the DCD file
                j = 0
                for frame in u.trajectory:
                        # keep track of time
                        if j == 0:
                                current_step += (total_steps_after_last_nc_len - current_step + steps_per_frame)
                        else:
                                current_step += steps_per_frame
                        frame_count += 1
                        resids_to_use = []
                        coor = protein.positions
                        # for each hydrophobic residue at this length
                        for resid in resids:
                                if coor[resid-1,0] >= xref:
                                        resids_to_use.append(resid)
                        print 'hydrophobic outside tunnel r='+str(r), resids_to_use
                        # now we need to deal with the different possible shapes of the asa array and compute the sum
                        # asa files with a single line
                        if len(np.shape(asa)) == 1:
                                print 'SINGLE FRAME DCD FILE! r =',r, dcd, f_asa
                                subset_sasa = np.sum(asa[resids_to_use], dtype=np.float64)
                                j += 1
                        # all other asa files
                        elif len(np.shape(asa)) == 2:
                                subset_sasa = np.sum(asa[j, resids_to_use], dtype=np.float64)
                                j += 1
                        else:
                                print 'Shape of input data array is not recognized.'
                                sys.exit()
                        pline = '%.9f' %(float(current_step)*0.015*0.001)+'\t'+str(r)+'\t'+str(frame_count)+\
                                '\t'+'%.9f' %subset_sasa #+'\t'+'%.9f' %(((subset_sasa/xstal_val)-1.0)*100.)+'\n'
                        if comp:
                                if folded_val == np.float64('0.000000000'):
                                        if np.abs(subset_sasa - folded_val) < 1E-9:
                                                pline += '\t'+'%.9f' %(0.0)
                                        else:
                                                # if the function cannot be computed or passed to zero for a reasonable, physically
                                                # justifiable reason, print a warning and then go to the next frame (N.B., no data
                                                # point is printed in this situation!
                                                print 'FAILED TO GET REAL NUMBER VALUE!'
                                                print pline, folded_val
                                                continue
                                                #sys.exit()
                                else:
                                        pline += '\t'+'%.9f' %(((subset_sasa/folded_val)-1.0)*100.)
                        ofile.write(pline+'\n')

        # if this dcd file (and, therefore, the asa file as well) is empty
        if file_size_zero:
                current_step += (total_steps_after_last_nc_len - current_step + steps_per_frame)
        total_steps_after_last_nc_len += steps_this_nc_len
