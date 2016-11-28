# -*- coding: utf-8 -*-

import numpy as np

computed_residue = input("Which residue was energy computed for?: ")

KL_list = []
for i in xrange(1, 299):
    if i == computed_residue:
        continue
    else:
        KL = np.loadtxt("residue_%s.dat" % i)
        KL_list.append(KL[1])
        
for i in xrange(1, 298):
    if i != computed_residue:
        continue
    else: 
        KL_list.insert((int(computed_residue))-1 , 0.0)
 
KL_list_5dp = [ '%.5f' % elem for elem in KL_list ]
print KL_list_5dp

np.savetxt('KL_all_%s.dat' % computed_residue , KL_list_5dp, fmt="%s")
