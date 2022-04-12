import os, sys
import numpy as np
import MDAnalysis as MDA
import numpy as np

# PURPOSE: to unify the various SASA calculations used in this project previously into a single program
#          to streamline analysis moving forward. Note that Ian's composition of the zeta metric will be used

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

# function to get residues implicated in DnaK interactions
def id_DnaK_res(pdb):

        # pdb: str, PDB ID of the current protein being analyzed

        dnak_resids = []

        resid_info = read_file('chaperones/DnaK/Limbo_output_parsed.dat')

        # get list of resids to be included in the analysis
        for line in resid_info:
                line = line.split()
                if line[0] == pdb:
                        if 'nothing' in line:
                                print 'No residues identified for this protein and analysis type.'
                                print 'NORMAL TERMINATION'
                                sys.exit()
                        for j in range (1, len(line)):
                                #print line[j].replace(',', '').split('-')
                                temp = line[j].replace(',', '').split('-')
                                start, end = int(temp[0]), int(temp[1])
                                for resid in range (start, end+1):
                                        # AVOID OVERCOUNTING
                                        if resid not in dnak_resids:
                                                dnak_resids.append(resid)
                        break

        return dnak_resids

# function to get residues implicated in aggregation
def id_Agg_res(pdb):

        # pdb: str, PDB ID of the current protein being analyzed

        agg_resids = []

        resid_info = read_file('aggregation/AMYLPRED2.dat')

        # get list of resids to be included in the analysis
        for line in resid_info:
                line = line.split()
                if line[0] == pdb:
                        if 'nothing' in line:
                                print 'No residues identified for this protein and analysis type.'
                                print 'NORMAL TERMINATION'
                                sys.exit()
                        for j in range (1, len(line)):
                                #print line[j].replace(',', '').split('-')
                                temp = line[j].replace(',', '').split('-')
                                start, end = int(temp[0]), int(temp[1])
                                for resid in range (start, end+1):
                                        # AVOID OVERCOUNTING
                                        if resid not in agg_resids:
                                                agg_resids.append(resid)
                        break

        return agg_resids

def str2bool(s):
        if s == 'True':
                return True
        else:
                return False

# usage statement
usage = 'python python_scripts/traj_sasa_master.py [cntrl file] [trajectory index] [mode; choose from hydrophobic, DnaK, Agg] [co-translational?]'

if len(sys.argv) != 5:
        print usage
        sys.exit()

cntrl = read_file(sys.argv[1])

pdb = sys.argv[1].split('/')[-1].split('_')[0]

xstal = '0' # whether or not this is a xstal simulation

traj = sys.argv[2]

mode = sys.argv[3]

cot = str2bool(sys.argv[4])

if mode not in ['hydrophobic', 'DnaK', 'Agg']:
        print mode, 'is not recognized as an acceptable mode of analysis'
        sys.exit()

# loop over the cntrl file and get parameters as needed
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

        # save the left part as the variable/flag name and the right part
        # as the variable/flag value
        var_name, var_value  = line[0].strip(), line[1].strip()

        # determine whether or not this run is for a xstal structure simulation
        if var_name == 'xstal':
                if var_value == '1':
                        xstal = '1'
        # get the length of the entire protein as an integer
        elif var_name == 'max_res':
                max_res = int(var_value)
        else:
                pass

# if xstal == '1' cot MUST be False
if xstal == '1' and cot == True:
        print 'If xstal == 1 cot MUST be False'
        sys.exit()

# if xstal == '1' and xstal_comp == False then we will compute A(subset) only
# if xstal == '1' and xstal_comp == True then we will compute A(subset) as well as A(subset)/<A(subset)>
# if xstal == '0' and xstal_comp == True then we will compute A(subset) as well as A(subset)/<A(subset)>
# if xstal == '0' and xstal_comp == False then we will die because that is a useless calculation
xstal_comp = True

# we have not parsed the cntrl file and loaded the command line arguments, so now we can determine
# the set of residues that we are going to analyze
if mode == 'hydrophobic':
        # get path to pdb format file for this protein from atomistic_pdbs/
        pdb_file               = os.popen('ls -1v atomistic_pdbs/'+pdb+'_*.pdb').read().split()[0]
        # pass PDB through function to ID hydrophobic residues by their resid
        analysis_resids       = id_hydrophob_res(pdb_file)

elif mode == 'DnaK':
        analysis_resids = id_DnaK_res(pdb)

elif mode == 'Agg':
        analysis_resids = id_Agg_res(pdb)

else:
        print ('Failed to determine the correct mode for this analysis')
        sys.exit()

# get list of data files to be analyzed for this protein and trajectory; 
# naming depends on whether or not xstal == '1' or '0'
if xstal == '1' and cot == False:

        # make a list of data files
        data_files       = os.popen('ls -1v analysis/'+pdb+'/asa/'+traj+\
                           '_i*_xstal_'+pdb+'_*sasa.txt').read().split()

# if this is not a xstal run we need to do some other junk, too
elif xstal == '0' and cot == False:

        # make a list of data files
        data_files       = os.popen('ls -1v analysis/'+pdb+'/asa/'+traj+\
                           '_i*_post_trans_'+pdb+'_*sasa.txt').read().split()

elif xstal == '0' and cot == True:

        print 'cot = True functionality is not yet implemented in this program'
        sys.exit()

else:
        print 'Failed to generate a list of files for analysis'
        sys.exit()

print 'Analysis mode: ', mode
print 'Residues considered for analysis: ', analysis_resids

# set up output file path and get path tot he file containing time-averaged values for the xstal simulations
if xstal == '1' and mode == 'DnaK':
        if xstal_comp:
                outpath = 'analysis/'+pdb+'/'+traj+'_DnaK_sasa_xstal_comp.txt'
                xstal_comp_inp = 'xstal_DnaK_sasa.dat'
        else:
                outpath = 'analysis/'+pdb+'/'+traj+'_DnaK_sasa_xstal.txt'
elif xstal == '0' and mode == 'DnaK':
        if xstal_comp:
                outpath = 'analysis/'+pdb+'/'+traj+'_DnaK_sasa.txt'
                xstal_comp_inp = 'xstal_DnaK_sasa.dat'
        else:
                print 'If xstal == 0 xstal_comp must == True'
                sys.exit()
elif xstal == '1' and mode == 'Agg':
        if xstal_comp:
                outpath = 'analysis/'+pdb+'/'+traj+'_agg_sasa_xstal_comp.txt'
                xstal_comp_inp = 'xstal_agg_sasa.dat'
        else:
                outpath = 'analysis/'+pdb+'/'+traj+'_agg_sasa_xstal.txt'
elif xstal == '0' and mode == 'Agg':
        if xstal_comp:
                outpath = 'analysis/'+pdb+'/'+traj+'_agg_sasa.txt'
                xstal_comp_inp = 'xstal_agg_sasa.dat'
        else:
                print 'If xstal == 0 xstal_comp must == True'
                sys.exit()
elif xstal == '1' and mode == 'hydrophobic':
        if xstal_comp:
                outpath = 'analysis/'+pdb+'/'+traj+'_hydrophobic_sasa_xstal_comp.txt'
                xstal_comp_inp = 'xstal_hydrophobic_sasa.dat'
        else:
                outpath = 'analysis/'+pdb+'/'+traj+'_hydrophobic_sasa_xstal.txt'
elif xstal == '0' and mode == 'hydrophobic':
        if xstal_comp:
                outpath = 'analysis/'+pdb+'/'+traj+'_hydrophobic_sasa.txt'
                xstal_comp_inp = 'xstal_hydrophobic_sasa.dat'
        else:
                print 'If xstal == 0 xstal_comp must == True'
                sys.exit()
else:
        print 'Output file path could not be generated from given inputs.'
        sys.exit()

if os.path.exists(outpath) == True:
        print 'Output file '+outpath+' already exists - delete it manually to continue.'
        sys.exit()

print 'Data will be written to', outpath

# if xstal_comp == True, then we need to get the value of <A(subset)> to use for this protein and mode
if xstal_comp:
        xstal_comp_val = 0
        raw = read_file(xstal_comp_inp)
        for line in raw:
                line = line.split()
                if line[0] == pdb:
                        xstal_comp_val = np.float64(line[1])
                        break
        if xstal_comp_val == '0':
                print 'Failed to find correct xstal_comp_val for', pdb, 'in', xstal_comp_inp
                sys.exit()
        print 'The value <A(subset)> will be set to', xstal_comp_val

# loop over the SASA data files for this trajectory and analysis method
file_cnt = 1
current_time = 5000.
ofile = open(outpath, 'w')
for data_file in data_files:

        data = np.loadtxt(data_file.strip())
        # some files contain only a single line, which is handled by this conditional
        if len(np.shape(data)) == 1:
                subset_total = np.sum(data[analysis_resids], dtype=np.float64)
                pline = str(file_cnt)+'\t'+'%.9f' %(current_time*0.015*0.001)+'\t'+'%.9f' %subset_total
                if xstal_comp:
                        pline += '\t'+'%.9f' %(((subset_total/xstal_comp_val)-1.)*100.)
                current_time += 5000.
                ofile.write(pline+'\n')
        # other files are a set of 2000 or more lines to be analyzed
        elif len(np.shape(data)) == 2:
                for j in range (0, len(data)):
                        subset_total = np.sum(data[j, analysis_resids], dtype=np.float64)
                        pline = str(file_cnt)+'\t'+'%.9f' %(current_time*0.015*0.001)+'\t'+'%.9f' %subset_total
                        if xstal_comp:
                                pline += '\t'+'%.9f' %(((subset_total/xstal_comp_val)-1.)*100.)
                        current_time += 5000.
                        ofile.write(pline+'\n')
        else:
                print 'Shape of input data array is not recognized.'
                sys.exit()
        file_cnt += 1

print 'NORMAL TERMINATION'
