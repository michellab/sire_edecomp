import numpy as np
import sys

# bash script will input 1 - 299 as sys.argv[1]
residue = int(sys.argv[1])
print "Summing distribution for residue %s" % residue

first_residue = int(sys.argv[2])
last_residue = int(sys.argv[3])

bin_list = np.empty((0,1000))

for j in range(first_residue , (last_residue+1)):
        new = np.loadtxt("all_data_residue_%s/1_combined_energies/residue%s" % (j , residue))
        bin_list = np.vstack((bin_list , new[:,1]))

bin_list = bin_list.T
bin_sum = bin_list.sum(axis=1)
np.savetxt('OUTPUT/1_all_data/residue_%s.dat' % residue , bin_list)
np.savetxt('OUTPUT/2_sum_data/residue_%s.dat' % residue, bin_sum)
