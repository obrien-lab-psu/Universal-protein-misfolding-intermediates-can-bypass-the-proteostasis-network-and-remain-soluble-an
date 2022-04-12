import os, sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
start = datetime.now()
# PURPOSE: extrapolate the folding time of a given protein/domain/interface from the
#          double-exponential fit to the survival probability of the unfolded state as a
#          function of time

# Notes: the folding status of a trajectory is determined from Q mode analysis
#         thus, need to rerun calculations for 3gn5 and 2kfw which were accidentally removed
# this version incorporates G analysis as well

# fit function
def func(t, k1, k2, f1):
        return (f1*np.exp(-k1*t))+((1.0-f1)*np.exp(-k2*t))

# function to parse and load entanglement data
def loadGdata(f_path):
	if os.path.exists(f_path) != True:
		print f_path, 'does not exist!'
		sys.exit()
	fdata = open(f_path)
	data = fdata.readlines()
	fdata.close()
	temp = []
	for line in data:
		if line[0] == '#':
			continue
		if 'frame' in line:
                        continue
                line = line.split()
                temp.append(float(line[2])) # this gets the f2Gtot column
        return np.asarray(temp)

def combine_data(data1, data2):
        # data1 -> qmode data (from 7.5 ns to end by 0.15 ns)
        # data2 -> G data (from 7.5 ns to end by 7.5 ns)

        output = np.zeros((len(data2), 3))

        time = 7.5
        inc = 7.5
        count = 0
        for i in range (0, len(data1), 50):
                output[count, 0] = time
                output[count, 1] = data1[i]
                output[count, 2] = data2[count]
                #print (time, count, output[count, 0], output[count, 1], output[count, 2])
                count += 1
                time += inc

        return output

def make_ref_dict(f_path):

        # f_path: path to the file (sent by Ian) containing the mean values from the xstal structure simulations

        f_raw = open(f_path)
        raw = f_raw.readlines()
        f_raw.close()

        ref_dict = {}

        for line in raw:
                line = line.split()
                pdb = line[2].split('/')[0]
                val = float(line[-1])
                ref_dict[pdb] = val

        return ref_dict

pdbid = sys.argv[1]
struc = sys.argv[2]
n_traj = 50
max_length = 0
ref_dict = make_ref_dict('xstal_f2Gtot_means.out')

for traj in range (1, n_traj+1):

	# load both data series in initial format
        data1 = np.loadtxt('folding_times/'+pdbid+'/'+str(traj)+'_'+struc+'_Qi_fF_alg2_3sig_3tol_pt.txt')[:,1]
        data2 = loadGdata('folding_times/'+pdbid+'/'+str(traj)+'_post_trans_frames-all_conc_G.dat')

	# parse both data series to get values every 7.5 ns for each in a two-column array
	combined = combine_data(data1, data2)

        if len(combined) > max_length:
                max_length = len(combined)

        print 'Done getting length of trajectory', traj

num_folded = np.zeros((max_length))

for traj in range (1, n_traj+1):
        data1 = np.loadtxt('folding_times/'+pdbid+'/'+str(traj)+'_'+struc+'_Qi_fF_alg2_3sig_3tol_pt.txt')[:,1]
        data2 = loadGdata('folding_times/'+pdbid+'/'+str(traj)+'_post_trans_frames-all_conc_G.dat')
        combined = combine_data(data1, data2)
        for i in range (0, len(combined)):
                q_folded = False
                G_folded = False
                # does q_mode information say protein is folded here?
                if combined[i, 1] > 0.0:
                        q_folded = True
                # does G information say protein is folded here?
                if combined[i, 2] <= ref_dict[pdbid]:
                        G_folded = True
                # if both Q and G say protein is folded
                if q_folded and G_folded:
                        for j in range (i, len(num_folded)):
                                num_folded[j] += 1.0        
                        break
        print 'Done with trajectory', traj

surv_prob_U = np.zeros((max_length))
times = np.zeros((max_length))
current_time = 7.5
ofile = open('folding_times/survival_prob_ts/'+pdbid+'_'+struc+'_surv_prob_U_v2.txt', 'w')
for i in range (0, max_length):
        times[i] = current_time
        surv_prob_U[i] = (50. - num_folded[i])/50.
        ofile.write('%.0f' %i +'\t'+ '%.5f' %times[i] +'\t'+ '%.5f' %num_folded[i] +'\t'+ '%.5f' %surv_prob_U[i] + '\n')
        current_time += 7.5

# do the fitting
sigma = np.ones(len(times))
sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
popt, pcov = curve_fit(func, times, surv_prob_U, p0=[1.0, 1.0, 0.2], bounds=(0., [np.inf, np.inf, 1.]), sigma = sigma)
k1 = popt[0]
k2 = popt[1]
f1 = popt[2]
f2 = 1.0 - popt[2]
fit = (f1*np.exp(-k1*times))+((f2)*np.exp(-k2*times))
r = pearsonr(fit, surv_prob_U)

log = open('folding_times/fit_results_v2.log', 'a')
log.write(pdbid+'\t'+struc+'\t'+'%.9e' %k1+'\t'+'%.9e' %k2+'\t'+'%.9e' %f1 +'\t'+'%.9e' %f2+'\t'+'%.9f' %(r[0]**2.)+'\n')

plt.clf()
plt.plot(times, surv_prob_U, label='tf1='+'%.1e' %(1./k1)+', tf2='+'%.1e' %(1./k2))
plt.plot(times, fit)
plt.xlabel('Time since release from ribosome, ns')
plt.ylabel('Survival probability of U')
plt.ylim(-0.05, 1.05)
plt.legend(loc='upper right')
plt.savefig('folding_times/plots/'+pdbid+'_'+struc+'_surv_prob_U_v2.png')

print 'Run time:', datetime.now()-start
