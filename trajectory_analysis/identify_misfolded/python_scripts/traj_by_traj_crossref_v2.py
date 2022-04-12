import os
import sys
import numpy as np
from datetime import datetime
start = datetime.now()

### PURPOSE - this program parses various input files in order to determine which trajectories fall into various categories
###           this new version is meant to be loaded as a module to facilitate various analyses

### CLASSES

# this class organizes all decisions made for a particular trajectory about whether or not 
# it is misfolded, likely to interact with chaperones, aggregate, be degraded, and be non-functional
class trajectory:

        # Getting a protein trajectory that is true for all values in {misfolded, TF, DnaK, GroEL, agg, deg. fun}
        # are classified as being misfolded, likely to avoid proteostasis machinery, and to have perturbed function

        def __init__(self, pdb, traj, q_mode, G_mode, misfolded):
                self.pdb = pdb # string; the pdb id for this trajectory
                self.traj = traj # string; the unique trajectory index for this simulation. Run {1, ..., 50} for each pdbid
                self.q_mode = q_mode # Boolean; iff True, Q_mode analysis says trajectory is misfolded
                self.G_mode = G_mode # Boolean; iff True, G_mode analysis says trajectory is misfolded
                self.misfolded = misfolded # Boolean; iff True, trajectory is concluded to be misfolded

                self.TF = None # Boolean; iff True, trajectory predicted to bypass TF
                self.DnaK = None # Boolean; iff True, trajectory predicted to bypass DnaK
                self.GroEL = None # Boolean; iff True, trajectory predicted to bypass GroEL
                self.agg = None # Boolean; iff True, trajectory predicted not to aggregate
                self.deg = None # Boolean; iff True, trajectory predicted not to be targeted for rapod degradation
                self.func = None # Boolean; iff True, trajectory is predicted to have perturbed function
                self.mean_overall_Q = None # float; the mean overall Q in the final 10 ns of trajectory
                self.percent_change_overall_Q = None # float; the percent change in overall Q between mean in final 10 ns and native state mean

        def addGroEL(self, GroEL):
                self.GroEL = GroEL

        def addTF(self, TF):
                self.TF = TF

        def addDnaK(self, DnaK):
                self.DnaK = DnaK

        def addAgg(self, agg):
                self.agg = agg

        def addDeg(self, deg):
                self.deg = deg

        def addFunc(self, func):
                self.func = func
        # commented out 20-Dec-2020 by Dan Nissley
        #def addOverallQ(self, mean_overall_Q, percent_change):
        #        self.mean_overall_Q = mean_overall_Q
        #        self.percent_change_overall_Q = percent_change

### FUNCTIONS

# read files into memory as a list of lines
def readFile(f_path):
        if os.path.exists(f_path):
                f_temp = open(f_path)
                temp = f_temp.readlines()
                f_temp.close()
                return temp
        else:
                print f_path, 'does not exist.'
                sys.exit()

# to map the shifted trajectory indices for 1d2f and 1glf to match the 1-50 inclusive
# scheme used for all other simulation labels
def map_traj_IDs_1d2f_or_1glf(pdb, traj):

        # pdb: string; pdbid, either '1glf' or '1d2f'
        # traj: string; trajectory index in ORIGINAL non-sequential order

        # Note well, this function is only necessary because Ian used the original labels while I
        # remapped them to be {1, ... , 50} sequential

        if pdb == '1glf':
                new_labels = readFile('1glf_trajs.txt')
        elif pdb == '1d2f':
                new_labels = readFile('1d2f_trajs.txt')
        else:
                print 'Not needed for pdbid', pdb
                return
        label_dict = {}
        for i in range (0, len(new_labels)):
                label_dict[new_labels[i].strip()] = str(i+1)
        return label_dict[traj]

# run this ONCE to generate overall Q information for misfolded trajectories
def get_overall_Q_stats(trajectories, outfile):

        # trajectories: a list of trajectory class objects, one for each trajectory to be analyzed
        # outfile     : string, path to which output data will be written

        # read in the reference values, the mean overall Q of the native-state simulations
        f_xstal_ref = readFile('xstal_overall_q_means.dat')
        xstal_ref = {}
        for line in f_xstal_ref:
                line = line.split()
                xstal_ref[line[0]] = float(line[1])

        map_traj_to_overall_Q = {}
        ofile = open(outfile, 'w')
        for traj in trajectories:
                if traj.pdb == '4g36':
                        continue
                if traj.misfolded:
                        count = 0
                        infiles = os.popen('ls -1v analysis/'+traj.pdb+'/'+traj.traj+'_*post_trans*pt.txt').read().split()
                        length = int(os.popen('wc -l analysis/'+traj.pdb+'/'+traj.traj+'_*post_trans*pt.txt').read().split()[-2])
                        all_data = np.zeros((length))
                        for infile in infiles:
                                data = np.loadtxt(infile)
                                if len(np.shape(data)) == 1:
                                        all_data[count] = data[1]
                                        count += 1
                                else:
                                        for i in range (0, len(data)):
                                                all_data[count] = data[i,1]
                                                count += 1
                        mean = np.mean(all_data[-134:])
                        percent_change = percChange(xstal_ref[traj.pdb], mean, absol=True)
                        ofile.write(traj.pdb+'\t'+traj.traj+'\t'+'%.9f' %mean+'\t'+'%.9f' %percent_change+'\n')
                        #print traj.pdb, traj.traj, mean, percent_change
                        #traj.addOverallQ(mean, percent_change)
        return

# function that, given a list of pdb ids, formats them nicely for copy and paste into a Microsoft Word table
def format_for_tables(pdb_list):
        pline = ''
        for s in range (0, len(pdb_list)-1):
                pline += pdb_list[s].upper()+', '
        pline += pdb_list[-1].upper()
        return pline

# given a string, return Boolean equivalent
def str2bool(string):
        if string == 'True':
                return True
        else:
                return False

# get percent change between a reference value and val
def percChange(ref, val, absol=False):

        # if absol == True
        if absol:
                # use the absolute value in the numerator
                return np.abs(val-ref)/ref*100.
        else:
                # otherwise, do not use the absolute value
                return (val-ref)/ref*100.

# given a pdbid return the list of domains and interfaces for the protein
def get_strucs(pdbid):

        # pdbid: string, the name of the pdb to be assessed

        # returns: a list containing strings that correspond to the different structures for this protein
        #          these are of the form ['d1', 'd2', 'dom1_dom2]

        # ls -1v analysis/1cli/1_*Qi*fF*xstal.txt

        if pdbid in ['2kfw', '3gn5']:
                return ['d1', 'd2', 'dom1_dom2']

        if pdbid in ['4g36']:
                return ['d1', 'd2', 'd3', 'd4', 'dom1_dom2', 'dom1_dom3', 'dom1_dom4', 'dom2_dom3', 'dom2_dom4']
                #return ['d1', 'd2', 'd3', 'd4', 'dom1_dom2', 'dom1_dom3', 'dom2_dom3', 'dom2_dom4']

        files = os.popen('ls -1v analysis/'+pdbid+'/1_*Qi*fF*xstal.txt').read().split()
        strucs = []
        for f in files:
                f = f.split('/')[-1].split('_')
                if len(f) == 8:
                        strucs.append(f[1])
                elif len(f) == 9:
                        strucs.append(f[1]+'_'+f[2])
                else:
                        print 'Failed to identify whether this is a domain or interface file'
                        sys.exit()
        return strucs

# function to parse q_mode input data - only needs to be run once, then you can read in the data from file to save time
def get_q_mode_parsed(outname, time_cutoff_type, multiple=3.0):

        # outname: string; used as file name for the output
        # time_cutoff_type; string, choose from keys in column_dict - determines which data column in
        #                   mean_props_by_traj/+pdb+'_Qmode.txt' will be used for analysis
        # multiple: float, folding assumed to occur once this threshold is passed

        # output: a file with two columns. First has the form 'pdbid:trajectory index' and the second is True or False
        #         iff True, that trajectory was found to be misfolded by Qmode analysis

        # keys are of the form 1cli:5, pdb:trajectory, all as string
        q_mode_parsed = {}

        # key is pdbid
        q_mode_ref = {}

        # maps the chosen 
        column_dict = {'first100':3, 'first200':4, 'first300':5, 'first400':6, 'first500':7,
                       'last100':8,  'last200':9,  'last300':10,  'last400':11, 'last500':12}

        # get list of proteins to analyze from file
        pdbs = readFile('all_protein_list.txt')

        # load the reference mean Qmode values computed over the crystal structure reference simulations
        xstal_ref = readFile('xstal_stats.inp')
        for line in xstal_ref:
                line = line.split()
                if line[2] == 'Qi':
                        # key is pdbid+':'+'dom1_dom2', values are list containing [<xs>, xs_std]
                        q_mode_ref[line[0]+':'+line[1]] = line[3:5]

        for pdb in pdbs:
                pdb = pdb.strip()
                inpfile = readFile('mean_props_by_traj/'+pdb+'_Qmode.txt')
                strucs = get_strucs(pdb)
                for traj in range (1, 50+1):
                        misfolded = False
                        for struc in strucs:
                                ref_vals = q_mode_ref[pdb+':'+struc]
                                ref_mean = float(ref_vals[0])
                                ref_std = float(ref_vals[1])
                                threshold = ref_mean - (ref_std*multiple)
                                for line in inpfile:
                                        temp = line.split()
                                        if temp[1] == str(traj) and temp[2] == struc:
                                                value = float(temp[column_dict[time_cutoff_type]])
                                                if value < threshold:
                                                        misfolded = True
                        q_mode_parsed[pdb+':'+str(traj)] = misfolded
                print 'Done parsing Q mode information for', pdb, datetime.now() - start
        ofile = open(outname, 'w')
        for key in q_mode_parsed:
                ofile.write(key+' '+str(q_mode_parsed[key])+'\n')
                #print key, q_mode_parsed[key]
        return q_mode_parsed

# for ease of use and saving time - just load the data written to file by get_q_mode_parsed
def get_q_mode_parsed_from_file(fpath):

        # fpath: file to be read in and parsed

        # assumes data is in the exact format written by get_q_mode_parsed

        q_mode_parsed = {}

        temp = readFile(fpath)
        for line in temp:
                line = line.split()
                q_mode_parsed[line[0]] = str2bool(line[1])

        return q_mode_parsed

# function to parse the G metric files from Ian
def get_G_mode_parsed(fpath, G_mode_cutoff=0.1):

        # fpath: the path to the input G data from Ian
        # G_mode_cutoff: threshold for determining misfolded; iff abs(G) > G_mode_cutoff, trajectory is misfolded
        
        # keys are of the form 1cli:5, pdb:trajectory. Values are True if misfolded and False if folded
        G_mode_parsed = {}

        G_mode = readFile(fpath)

        for line in G_mode:
                misfolded = False
                if line[0] == '#':
                        continue
                temp = line.split(',')
                pdb = temp[0].strip()
                traj = temp[1].strip()
                temp2 = []
                for i in temp:
                        temp2.append(i.strip())
                g_vector = temp2[2:7]
                for g in g_vector:
                        if abs(float(g)) > G_mode_cutoff:
                                misfolded = True
                if pdb in ['1glf', '1d2f']:
                        # need to map these trajectory indices using the function for this purpose
                        G_mode_parsed[pdb+':'+map_traj_IDs_1d2f_or_1glf(pdb, traj)] = misfolded
                else:
                        G_mode_parsed[pdb+':'+traj] = misfolded
        return G_mode_parsed

# consider both q_mode_parsed and G_mode_parsed and make a determination for every trajectory for every protein
def find_misfolded(q_mode_parsed, G_mode_parsed):

        # q_mode_parsed: dictionary output by get_q_mode_parsed or get_q_mode_parsed_from_file
        # G_mode_parsed: dictionary output by get_G_mode_parsed

        # returns: list of trajectory class objects, number of misfolded trajectories, num Q misfolded,
        #          num G misfolded, num both misfolded

        q_result = False
        G_result = False
        misfolded = False
        q_misfolded = 0
        G_misfolded = 0
        both_misfolded = 0
        total_misfolded = 0
        trajectories = []

        for i in q_mode_parsed:
                misfolded = False
                if i not in G_mode_parsed:
                        print i, 'not found in G_mode_parsed - keys should be identical to q_mode_parsed'
                        sys.exit()
                pdb, traj_index = i.split(':')[0:2]
                q_result = q_mode_parsed[i]
                G_result = G_mode_parsed[i]
                if q_result:
                        misfolded = True
                        q_misfolded += 1
                if G_result:
                        misfolded = True
                        G_misfolded += 1
                if q_result and G_result:
                        both_misfolded += 1

                if misfolded:
                        total_misfolded += 1
                trajectories.append(trajectory(pdb, traj_index, q_result, G_result, misfolded))

        return trajectories, total_misfolded, q_misfolded, G_misfolded, both_misfolded

# determine the percent of trajectories misfolded for each protein
def get_percent_misfolded_by_protein(trajectories, percentage_list):

        # trajectories: list of trajectory class objects output by find_misfolded
        # percentage_list: list containing percentages, e.g. [90., 80, 50., 20.] - function
        #                  will count how many proteins are misfolded at least that percentage of the time

        # returns: dictionary mapping pdbid to percent misfolded, second dictionary mapping percent to list of pdbs
        #          with at least that percentage of misfolded across 50 trajectories

        # Note: this function ignores 4g36, returning results for the 122 E. coli proteome proteins only

        pdbs = []

        map_pdb_to_percent_misfolded = {}
        map_percent_to_proteins = {}

        for traj in trajectories:
                if traj.pdb == '4g36':
                        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        for pdb in pdbs:
                count = 0
                for traj in trajectories:
                        if traj.pdb == pdb:
                                if traj.misfolded:
                                        count += 1
                # compute the percent misfolded for this protein
                map_pdb_to_percent_misfolded[pdb] = float(count)/50.*100.

        for percentage in percentage_list:
                pdb_list = []
                for pdb in map_pdb_to_percent_misfolded:
                        if map_pdb_to_percent_misfolded[pdb] >= percentage:
                                pdb_list.append(pdb)
                map_percent_to_proteins['%.0f' %percentage] = pdb_list

        return map_pdb_to_percent_misfolded, map_percent_to_proteins

# function to parse DnaK interaction metric(s) and add information to trajectory class objects
def add_DnaK_info(trajectories, time_cutoff_type, dnak_threshold=10.0):
        
        # trajectories: list of trajectory class objects
        # time_cutoff_type: first100, last300, etc. - determines which column to use in input files
        # dnak_threshold: percentage cutoff - if DnaK metric is >dnak_threshold, considered to interact
        #                                     if DnaK metric is <=dnak_threshold, considered not to interact

        # returns: nothing, but modifies trajectory class object list

        # load list of proteins in 122 dataset that are thought to interact with DnaK
        dnak_clients = readFile('chaperones/DnaK/myproteins_that_are_DnaK_clients.txt')

        known_dnak_substrates = []
        for line in dnak_clients:
                known_dnak_substrates.append(line.split(';')[0])

        # note the -1 shift from the q_mode data (bc q_mode data include a column with domain/interface name
        column_dict = {'first100':2, 'first200':3, 'first300':4, 'first400':5, 'first500':6,
                       'last100':7,  'last200':8,  'last300':9,  'last400':10, 'last500':11}

        # keys are pdbid:traj_index, values are Boolean True if not likely to interact with DnaK and Boolean False
        # if interaction is expected
        dnak = {}

        pdbs = []

        for traj in trajectories:
                #if traj.pdb == '4g36':
                #        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        for pdb in pdbs:
                # if pdb not in known substrate list, assume no interactions
                if pdb not in known_dnak_substrates:
                        for traj in range(1, 50+1):
                                dnak[pdb+':'+str(traj)] = True
                # if pdb in the list of known substrates, check exposure of DnaK binding sites
                else:
                        inpdata = readFile('mean_props_by_traj/'+pdb+'_DnaK_sasa.txt')
                        for line in inpdata:
                                if line[0] == '#':
                                        continue
                                temp = line.split()
                                traj = temp[1]
                                value =  float(temp[column_dict[time_cutoff_type]])
                                if value > dnak_threshold:
                                        dnak[pdb+':'+traj] = False # large deviation, likely bound by DnaK
                                else:
                                        dnak[pdb+':'+traj] = True # small deviation, likely avoids DnaK binding

        # make corrections for the three proteins {2jrx, 2v81, 2kfw} that have no predicted DnaK binding sites per Limbo
        # 2jrx and 2kfw are on the Bukau client list while 2v81 is not.
        # thus, assume 2v81 does NOT interact with DnaK but 2jrx and 2kfw do interact
        for traj in range (1, 50+1):
                dnak['2jrx:'+str(traj)] = False
                dnak['2kfw:'+str(traj)] = False
                dnak['2v81:'+str(traj)] = True

        # update trajectory class object list
        for traj in trajectories:
                traj.addDnaK(dnak[traj.pdb+':'+traj.traj])

        # return - the trajectory list has been updated in place, no need to return it
        return

# function to parse GroEL interaction metric(s) and add information to trajectory class objects
def add_GroEL_info(trajectories, time_cutoff_type, groel_threshold=10.0):
        
        # trajectories: list of trajectory class objects
        # time_cutoff_type: first100, last300, etc.  - determines which column to use in input file
        # greol_threshold: percentage cutoff - if GroEL metric is >greol_threshold, considered to interact
        #                                      if GroEL metric is <= groel_threshold, considered not to interact

        # returns: nothing, but modifies list of trajectory class objects

        # load the list of 276 known GroEL substrates from the literature
        groel_clients = readFile('chaperones/GroEL_GroES/myproteins_that_are_substrates.dat')

        known_groel_substrates = []
        for line in groel_clients:
                known_groel_substrates.append(line.split()[2])

        # note the -1 shift from the q_mode data (bc q_mode data include a column with domain/interface name
        column_dict = {'first100':2, 'first200':3, 'first300':4, 'first400':5, 'first500':6,
                       'last100':7,  'last200':8,  'last300':9,  'last400':10, 'last500':11}

        # keys are pdbid:traj_index, values are Boolean True if not likely to interact with GroEL and Boolean False
        # if interaction is expected
        groel = {}

        # make list of pdbs
        pdbs = []
        for traj in trajectories:
                #if traj.pdb == '4g36':
                #        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        for pdb in pdbs:
                # if pdb not a known GroEL substrate assume it does not interact
                if pdb not in known_groel_substrates:
                        for traj in range (1, 50+1):
                                groel[pdb+':'+str(traj)] = True
                # otherwise, check exposed hydrophobic surface area
                else:
                        inpdata = readFile('mean_props_by_traj/'+pdb+'_hydrophobic_sasa.txt')
                        for line in inpdata:
                                if line[0] == '#':
                                        continue
                                temp = line.split()
                                value = float(temp[column_dict[time_cutoff_type]])
                                traj = temp[1]
                                if value > groel_threshold:
                                        groel[pdb+':'+traj] = False # large deviation, likely binds GroEL
                                else:
                                        groel[pdb+':'+traj] = True # small deviation, likely does NOT bind GroEL
        for traj in trajectories:
                traj.addGroEL(groel[traj.pdb+':'+traj.traj])

        # return - the trajectory list has been updated in place, no need to return it
        return

# function to parse TF interaction metric(s) and add information to trajectory class objects
def add_TF_info(trajectories, trigger_factor_threshold=10.0):

        # trajectories: list of trajectory class objects
        # trigger_factor_threshold: percentage cutoff - if trigger factor metric is >trigger_factor_threshold, considered to interact
        #                                             - if <= trigger_factor_threshold, considered not to interact        

        # Note well, trigger factor information can only be computed for proteins that have a folded population, as those
        #            folded trajectories need to be used as the reference state. Also, a few proteins have one or more nascent
        #            chain lengths with no data (e.g., single trajectory in {F} and the random number of integration time steps
        #            selected is <5000, leading to no data for reference state at that nascent chain length).

        # TF prefers to bind proteins at nascent chain lengths ~100 or greater;
        # these proteins are too short, so ignore them
        too_short = ['1ah9', '1dcj', '1dfu', '1fm0', '1h75', '1jns',
                     '1jw2', '1sg5', '1t8k', '1xn7', '2a6q', '2axd',
                     '2hd3', '2jee', '2jrx', '3iv5']

        trigger_factor = {}

        # make list of pdbs
        pdbs = []
        for traj in trajectories:
                #if traj.pdb == '4g36':
                #        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        # initially, set all protein and trajectories to False (i.e., likely bind TF)
        for pdb in pdbs:
                for traj in range (1, 50+1):
                        trigger_factor[pdb+':'+str(traj)] = False

        # account for the short proteins that likely cannot bind trigger factor
        for pdb in too_short:
                for traj in range (1, 50+1):
                        trigger_factor[pdb+':'+str(traj)] = True

        # now, flip trajectories from False to True based on information in input files
        for pdb in pdbs:
                tempfile = 'mean_props_by_traj/'+pdb+'_cot_hydrophobic_sasa.txt'
                if os.path.exists(tempfile):
                        inpdata = readFile(tempfile)
                        for line in inpdata:
                                if line[0] == '#':
                                        continue
                                temp = line.split()
                                value = float(temp[1])
                                traj = temp[0]
                                if value > trigger_factor_threshold:
                                        trigger_factor[pdb+':'+traj] = False
                                else:
                                        trigger_factor[pdb+':'+traj] = True

        for traj in trajectories:
                traj.addTF(trigger_factor[traj.pdb+':'+traj.traj])

        # return - the trajectory list has been updated in place, no need to return it 
        return

# function to parse aggregation metric and add information to trajectory class objects
def add_aggregation_info(trajectories, time_cutoff_type, agg_threshold=10.0):

        # trajectories: list of trajectory class objects
        # time_cutoff_type: first100, last300, etc.  - determines which column to use in input file
        # agg_threshold: percentage cutoff - if agg metric is >agg_threshold, considered to aggregate
        #                                    if agg metric is <= agg_threshold, considered not to aggregate

        # returns: nothing, but modifies list of trajectory class objects

        # note the -1 shift from the q_mode data (bc q_mode data include a column with domain/interface name
        column_dict = {'first100':2, 'first200':3, 'first300':4, 'first400':5, 'first500':6,
                       'last100':7,  'last200':8,  'last300':9,  'last400':10, 'last500':11}

        # keys are pdbid:traj_index, values are Boolean True if not likely to aggregate and Boolean False
        # if interaction is expected
        aggregation = {}

        # make a list of pdbs
        pdbs = []
        for traj in trajectories:
                #if traj.pdb == '4g36':
                #        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        # get information for each protein and trajectory
        for pdb in pdbs:
                inpdata = readFile('mean_props_by_traj/'+pdb+'_agg_sasa.txt')
                for line in inpdata:
                        if line[0] == '#':
                                continue
                        temp = line.split()
                        value = float(temp[column_dict[time_cutoff_type]])
                        traj = temp[1]
                        if value > agg_threshold:
                                aggregation[pdb+':'+traj] = False # large deviation, likely to aggregate
                        else:
                                aggregation[pdb+':'+traj] = True # small deviation, likely not to aggregate

        for traj in trajectories:
                traj.addAgg(aggregation[traj.pdb+':'+traj.traj])

        # return - the trajectory list has been updated in place, do not need to return the list
        return

# function to parse degradation metric and add informaton to trajectory class objects
def add_degradation_info(trajectories, time_cutoff_type, deg_threshold=10.0):

        # trajectories: list of trajectory class objects
        # time_cutoff_type: first100, last300, etc.  - determines which column to use in input file
        # deg_threshold: percentage cutoff - if agg metric is >deg_threshold, considered to degrade
        #                                    if agg metric is <= deg_threshold, considered not to be degraded

        # returns: nothing, but modifies list of trajectory class objects

        # note the -1 shift from the q_mode data (bc q_mode data include a column with domain/interface name
        column_dict = {'first100':2, 'first200':3, 'first300':4, 'first400':5, 'first500':6,
                       'last100':7,  'last200':8,  'last300':9,  'last400':10, 'last500':11}

        # keys are pdbid:traj_index, values are Boolean True if not likely to be degraded and Boolean False
        # if interaction is expected
        degradation = {}

        # make a list of pdbs
        pdbs = []
        for traj in trajectories:
                #if traj.pdb == '4g36':
                #        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        # get information for each protein and trajectory
        for pdb in pdbs:
                inpdata = readFile('mean_props_by_traj/'+pdb+'_hydrophobic_sasa.txt')
                for line in inpdata:
                        if line[0] == '#':
                                continue
                        temp = line.split()
                        value = float(temp[column_dict[time_cutoff_type]])
                        traj = temp[1]
                        if value > deg_threshold:
                                degradation[pdb+':'+traj] = False # large deviation, likely to be degraded
                        else:
                                degradation[pdb+':'+traj] = True # small deviation, likely not to be degraded

        for traj in trajectories:
                traj.addDeg(degradation[traj.pdb+':'+traj.traj])

        # return - the trajectory list has been updated in place, do not need to return the list
        return


# function to parse functional information and add it to trajectory class objects
def add_function_info(trajectories, time_cutoff_type, func_threshold=-10.0):

        # trajectories: list of trajectory class objects
        # time_cutoff_type: last100, last200, etc.
        # func_threshold: if chi function metric is <= func_threshold, function is taken to be perturbed
        #                                           >  func_threshold, function is taken not to be perturbed

        if func_threshold > 0.0:
                print 'This analysis defines functional loci to be perturbed if they have a negative value for the function metric'
                print 'Your selected threshold of', func_threshold, 'is greater than zero - this may lead to erroneous results'
                sys.exit()

        functional = {}

        column_dict = {'last100':2,  'last200':3,  'last300':4,  'last400':5, 'last500':6}

        pdbs = []
        for traj in trajectories:
                #if traj.pdb == '4g36':
                #        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)

        # not all proteins have functional loci - to avoid key errors, initially populate all as False
        for pdb in pdbs:
                for traj in range (1, 50+1):
                        functional[pdb+':'+str(traj)] = False
        for pdb in pdbs:
                inpdata = readFile('mean_props_by_traj/'+pdb+'_chi_func.txt')
                for line in inpdata:
                        if line[0] == '#':
                                continue
                        temp = line.split()
                        value = float(temp[column_dict[time_cutoff_type]])
                        traj = temp[1]
                        if value <= func_threshold:
                                functional[pdb+':'+traj] = True # large magnitude and negative, likely has perturbed function
                        else:
                                functional[pdb+':'+traj] = False # small negative value or positive, likely no perturbed function

        for traj in trajectories:
                traj.addFunc(functional[traj.pdb+':'+traj.traj])

        # no need to return anything, the trajectory class objects have been updated
        return

# function to generate input files with the TF, DnaK, GroEL/GroES, agg, deg, and func metrics for trajectories in the
# misfolded and folded populations
def get_M_and_F_values_to_histogram(trajectories):

        # trajectories: a list of trajectory class objects

        # NOTE WELL: assumes you want to use the last 100 ns data

        # make a list of pdbs and also convert trajectories list into a dictionary
        # keys   -> pdbid+':'+traj_index
        # values -> corresponding trajectory class object

        pdbs = []
        traj_dict = {}
        for traj in trajectories:
                if traj.pdb == '4g36':
                        continue
                if traj.pdb not in pdbs:
                        pdbs.append(traj.pdb)
                traj_dict[traj.pdb+':'+traj.traj] = traj

        # read in all of the analysis files needed
        tf_dict = {} # needs to be populated ahead of time; not every protein has this file
        dnak_dict = {}
        #groel_dict = {}
        agg_dict = {}
        #deg_dict = {}
        hydro_dict = {}
        func_dict = {}

        # some of these dictionaries need to be populated ahead of time
        for pdb in pdbs:
                for traj in range (1, 50+1):
                        func_dict[pdb+':'+str(traj)] = 'None'
                        tf_dict[pdb+':'+str(traj)] = 'None'
                        dnak_dict[pdb+':'+str(traj)] = 'None'

        # gather information for all proteins and trajectories
        for pdb in pdbs:

                # load TF data for all proteins that have any
                temp = 'mean_props_by_traj/'+pdb+'_cot_hydrophobic_sasa.txt'
                if os.path.exists(temp):
                        temp2 = open(temp)
                        tf = temp2.readlines()
                        temp2.close()
                        for line in tf:
                                if '#' in line:
                                        continue
                                line = line.split()
                                traj = line[0]
                                val = line[1]
                                tf_dict[pdb+':'+traj] = val

                # load DnaK data for all trajectories
                temp = 'mean_props_by_traj/'+pdb+'_DnaK_sasa.txt'
                if os.path.exists(temp):
                        temp2 = open(temp)
                        dnak = temp2.readlines()
                        temp2.close()
                        for line in dnak:
                                if '#' in line:
                                        continue
                                line = line.split()
                                traj = line[1]
                                #print pdb, traj
                                val = line[7]
                                dnak_dict[pdb+':'+traj] = val

                # load the aggregation data for all trajectories
                temp = open('mean_props_by_traj/'+pdb+'_agg_sasa.txt')
                agg = temp.readlines()
                temp.close()
                for line in agg:
                        if '#' in line:
                                continue
                        line = line.split()
                        traj = line[1]
                        val = line[7]
                        agg_dict[pdb+':'+traj] = val

                # load post-translational hydrophobic surface area data for all proteins
                temp = open('mean_props_by_traj/'+pdb+'_hydrophobic_sasa.txt')
                hydro = temp.readlines()
                temp.close()
                for line in hydro:
                        if '#' in line:
                                continue
                        line = line.split()
                        traj = line[1]
                        val = line[7]
                        hydro_dict[pdb+':'+traj] = val

                # load functional values for all trajectories
                temp = 'mean_props_by_traj/'+pdb+'_chi_func.txt'
                if os.path.exists(temp):
                        temp2 = open(temp)
                        func = temp2.readlines()
                        temp2.close()
                        for line in func:
                                if '#' in line:
                                        continue
                                line = line.split()
                                traj = line[1]
                                val = line[2]
                                func_dict[pdb+':'+traj] = val

        # open all of the output files we need
        tf_M = open('traj_by_traj_TF_values_misfolded.txt', 'w')
        tf_F = open('traj_by_traj_TF_values_folded.txt', 'w')

        dnak_M = open('traj_by_traj_DnaK_values_misfolded.txt', 'w')
        dnak_F = open('traj_by_traj_DnaK_values_folded.txt', 'w')

        #groel_M = open('traj_by_traj_GroEL_values_misfolded.txt', 'w')
        #groel_F = open('traj_by_traj_GroEL_values_folded.txt', 'w')

        hydro_M = open('traj_by_traj_hydrophobic_values_misfolded.txt', 'w')
        hydro_F = open('traj_by_traj_hydrophobic_values_folded.txt', 'w')

        agg_M = open('traj_by_traj_agg_values_misfolded.txt', 'w')
        agg_F = open('traj_by_traj_agg_values_folded.txt', 'w')

        #deg_M = open('traj_by_traj_deg_values_misfolded.txt', 'w')
        #deg_F = open('traj_by_traj_deg_values_folded.txt', 'w')

        func_M = open('traj_by_traj_func_values_misfolded.txt', 'w')
        func_F = open('traj_by_traj_func_values_folded.txt', 'w')

        for pdb in pdbs:
                for traj in range (1, 50+1):
                        name = pdb+':'+str(traj)
                        if traj_dict[name].misfolded:
                                # get TF value
                                tf_M.write(pdb+' '+str(traj)+' '+tf_dict[name]+'\n')

                                # get DnaK value
                                dnak_M.write(pdb+' '+str(traj)+' '+dnak_dict[name]+'\n')

                                # get post-translational hydrophobic value
                                hydro_M.write(pdb+' '+str(traj)+' '+hydro_dict[name]+'\n')

                                # get agg value
                                agg_M.write(pdb+' '+str(traj)+' '+agg_dict[name]+'\n')

                                # get func value
                                func_M.write(pdb+' '+str(traj)+' '+func_dict[name]+'\n')

                        else:
                                # get TF value
                                tf_F.write(pdb+' '+str(traj)+' '+tf_dict[name]+'\n')

                                # get DnaK value
                                dnak_F.write(pdb+' '+str(traj)+' '+dnak_dict[name]+'\n')

                                # get post-translational hydrophobic value
                                hydro_F.write(pdb+' '+str(traj)+' '+hydro_dict[name]+'\n')

                                # get agg value
                                agg_F.write(pdb+' '+str(traj)+' '+agg_dict[name]+'\n')

                                # get func value
                                func_F.write(pdb+' '+str(traj)+' '+func_dict[name]+'\n')
        return

# generate lists of trajectories in various categories for the construction of Venn diagrams
def get_venn_data(trajectories):

        # trajectories: list of trajectory class objects

        # Used to generate file containing list of trajectories, referred to as
        # "pdbid:trajectory_index", that are misfolded and expected not to interact with:
        # (1) trigger factor
        # (2) DnaK
        # (3) GroEL/GroES
        # (4) aggregation
        # (5) degradation
        # (6) function 

        ofile = open('venn_diagram_input.txt', 'w')
        ofile.write('pdb, traj, TF, DnaK, GroEL, aggregation, degradation, function\n')
        for traj in trajectories:
                if traj.pdb == '4g36':
                        continue
                if traj.misfolded:
                        pline  = traj.pdb+','+traj.traj+','+str(traj.TF)+','+str(traj.DnaK)+','
                        pline += str(traj.GroEL)+','+str(traj.agg)+','+str(traj.deg)+','+str(traj.func)+'\n'
                        ofile.write(pline)
