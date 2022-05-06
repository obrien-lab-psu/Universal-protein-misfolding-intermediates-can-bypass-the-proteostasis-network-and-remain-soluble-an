import os
import sys
from datetime import datetime
sys.path.append('python_scripts')
import traj_by_traj_crossref_v2 as crossref
start = datetime.now()

### GET Q_MODE AND G_MODE DATA FOR EACH PROTEIN AND TRAJECTORY

# generate q_mode_parsed file for last100ns and get the same info as a dictionary
q_mode = crossref.get_q_mode_parsed('traj_by_traj_Qmode_misfolded_trajs.txt', 'last100', multiple=3.0)

# generate G_mode dictionary 
G_mode = crossref.get_G_mode_parsed('entanglement_summary_files/mature_100ns_xstal_corrected_G0-5_only_xstal_corrected.txt', G_mode_cutoff=0.1)

trajectories, N_misfolded, N_Q_misfolded, N_G_misfolded, N_both_misfolded = crossref.find_misfolded(q_mode, G_mode)

# note that these numbers include 4g36's 23 misfolded trajectories
print
print 'Number of misfolded trajectories:', N_misfolded
print 'Number of Q mode misfolded trajectories:', N_Q_misfolded
print 'Number of G misfolded trajectories:', N_G_misfolded
print 'Number of Q and G misfolded trajectories:', N_both_misfolded
print 'A total of', len(trajectories), 'have been considered in this analysis.'
print

### UPDATE TRAJECTORIES WITH TF, DNAK, GROEL, AGGREGATION, DEGRADATION, AND FUNCTION INFORMATION

# add trigger factor information
crossref.add_TF_info(trajectories, trigger_factor_threshold=10.0)
print 'Done adding trigger factor information'

# add DnaK information
crossref.add_DnaK_info(trajectories, 'last100', dnak_threshold=10.0)
print 'Done adding DnaK information'

# add GroEL information
crossref.add_GroEL_info(trajectories, 'last100', groel_threshold=10.0)
print 'Done adding GroEL information'

# add aggregation information
crossref.add_aggregation_info(trajectories, 'last100', agg_threshold=10.0)
print 'Done adding aggregation information'

# add degradation information
crossref.add_degradation_info(trajectories, 'last100', deg_threshold=10.0)
print 'Done adding degradation information'

# add function information
crossref.add_function_info(trajectories, 'last100', func_threshold=-10.0)
print 'Done adding function information'

print

### GENERATE INPUT FILES FOR HISTOGRAMMING METRICS FOR {F} AND {M} POPULATIONS
crossref.get_M_and_F_values_to_histogram(trajectories)
print 'Done generating comparison metric histogram files'
print 

### GENERATE VENN DIAGRAM INPUT FILE
crossref.get_venn_data(trajectories)
print 'Done generating Venn diagram input file'
print

### PERCENT MISFOLDING BY PROTEIN 

# generating percent misfolded for each protein
percentage_list = [100.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 2.0]
map_pdb_to_percent_misfolded, map_percent_to_proteins = crossref.get_percent_misfolded_by_protein(trajectories, percentage_list)

for percentage in percentage_list:
        percentage = '%.0f' %percentage
        print 'A total of', len(map_percent_to_proteins[percentage]), 'proteins are', percentage, 'percent misfolded or more'
        print crossref.format_for_tables(map_percent_to_proteins[percentage])
        print

ofile = open('traj_by_traj_percent_misfolded_by_pdb.txt', 'w')
for pdb in map_pdb_to_percent_misfolded:
        ofile.write(pdb+'\t'+'%.1f' %map_pdb_to_percent_misfolded[pdb]+'\n')

print

### GET OVERALL Q AND PERCENT DIFFERENCE FROM NATIVE STATE VALUES

print 'Generating overall Q and percent change in overall Q data file'
# we cannot supply all of the necessary data files for this function, but we do supply its
# output at file name traj_by_traj_verall_Q_and_perc_change.txt
#crossref.get_overall_Q_stats(trajectories, 'traj_by_traj_overall_Q_and_perc_change_v2.txt')

### GET FOLDING TIME INFORMATION

# See python_scripts/tf_from_double-exponential_fit_v2.py and plots/manuscript/histogram_folding_times.py

### GET NUMBER OF MISFOLDED TRAJECTORIES NOT EXPECTED TO INTERACT WITH TF
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.TF:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and not expected to be bound by TF co-translationally'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET NUMBER OF MISFOLDED TRAJECTORIES NOT EXPECTED TO INTERACT WITH DNAK
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.DnaK:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and not expected to be bound by DnaK'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET THE NUMBER OF MISFOLDED TRAJECTORIES NOT EXPECTED TO INTERACT WITH GROEL
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.GroEL:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and not expected to be bound by GroEL'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET NUMBER OF MISFOLDED TRAJECTORIES NOT EXPECTED TO BIND ANY CHAPERONES
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.TF and traj.DnaK and traj.GroEL:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and not expected to be bound by TF, DnaK, or GroEL'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET NUMBER OF MISFOLDED TRAJECTORIES NOT EXPECTED TO AGGREGATE
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.agg:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and not expected to aggregate'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET NUMBER OF MISFOLDED TRAJECTORIES NOT EXPECTED TO BE DEGRADED
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.deg:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and not expected to be degraded'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET NUMBER OF MISFOLDED TRAJECTORIES WITH PERTURBED FUNCTION
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.func:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and expected to exhibit decreased function'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)
print

### GET NUMBER OF MISFOLDED TRAJECTORIES THAT BYPASS ALL PROTEOSTASIS MACHINERY
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.TF and traj.DnaK and traj.GroEL and traj.agg and traj.deg:
                count += 1
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print count, 'trajectories are misfolded and expected to avoid all proteostasis machinery'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)

### GET NUMBER OF MISFOLDED TRAJECTORIES THAT DO NOT BIND CHAPERONES, AGGREGATE, BE DEGRADED, AND HAVE PERTURBED FUNCTION
count = 0
unique_pdbs = []
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.TF and traj.DnaK and traj.GroEL and traj.agg and traj.deg and traj.func:
                count += 1
                #print traj.traj, traj.pdb
                if traj.pdb not in unique_pdbs:
                        unique_pdbs.append(traj.pdb)
print
print count, 'trajectories are misfolded and expected to avoid all proteostasis machinery and exhibit reduced function'
print 'These trajectories represent', len(unique_pdbs), 'unique proteins:'
print crossref.format_for_tables(unique_pdbs)

### write data file used to generate Figure 4
index = 1
ofile = open('escaping_trajectories.txt', 'w')
for traj in trajectories:
        if traj.pdb == '4g36':
                continue
        if traj.misfolded and traj.TF and traj.DnaK and traj.GroEL and traj.agg and traj.deg and traj.func:
                ofile.write(str(index)+' '+str(traj.traj)+' '+ traj.pdb +' '+ str(traj.G_mode)+'\n')
                index += 1
ofile.close()

### DONE
print
print 'EXECUTION TIME:', datetime.now()-start
