import os
import sys
import numpy as np

pdb  = sys.argv[1] # from all_protein_list.txt
mode = sys.argv[2] # ['both', 'intf', 'sm']
traj = sys.argv[3] # between 1 and 50 except for 1glf and 1d2f

# load trajectory to be processed; 1_func_sm_chi_xstal.txt
ogdata = np.loadtxt('analysis/'+pdb+'/'+traj+'_func_'+mode+'_chi.txt')

# get mean value from file
temp = open('xstal_func_'+mode+'.dat')
xsmeans = temp.readlines()
temp.close()
xsmean = None
for line in xsmeans:
        line = line.split()
        if line[0] == pdb:
                xsmean = np.float64(line[1])   
if xsmean == None:
        print 'This protein does not have data for this analysis mode - DONE.'
        sys.exit()
ofile = 'analysis/'+pdb+'/'+traj+'_func_'+mode+'_chi_processed.txt'
if os.path.exists(ofile) == True:
        print ofile, 'already exists - delete it manually to continue.'
        sys.exit()
# compute new quantity as chi(t)/meanXSchi - 1 * 100%
for i in range (0, len(ogdata)):
        oline = '%.9f' %ogdata[i,0]+'\t'+'%.5f' %ogdata[i,1]+'\t'+'%.5f' %(ogdata[i,1]/xsmean)+'\t'+'%.5f' %((1.0-(ogdata[i,1]/xsmean))*100.)+'\n'
        with open (ofile, "a") as out:
                out.write(oline)
