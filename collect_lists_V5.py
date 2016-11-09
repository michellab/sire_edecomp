#!/usr/bin/python
#
# TODO 
# - compute standard deviation of averages
#

import os,sys, pickle
from math import sqrt
import re
import numpy as np

if __name__ == "__main__":
    
    try:
        infolder = sys.argv[1]
    except IndexError:
        print ("Usage is %s energies_folder" % ( sys.argv[0] ))
        sys.exit(-1)
    
    files = os.listdir(infolder)
    files.sort()

    chunk0 = os.path.join(infolder,files[0])
    stream = open(chunk0,'rb')
    energiesDict = pickle.load(stream)
    stream.close()

    number_of_chunks = len(files)
#    print (number_of_chunks)

    # Make sure files are ordered by chunks
    for f in files[1:]:
        chunk = os.path.join(infolder,f)
        stream = open(chunk,'rb')
        data = pickle.load(stream)
        stream.close()

        for key, value in data.items():
            if key == "intraresidues" or key == "interresidues":
                continue
            for v in value:
                energiesDict[key].append(v)
    
    residue_number = infolder.split("-")[1].strip("/")    
 
   # make folder for residue data to be stored
    distribution_outdir = "all_data_residue_%s" % residue_number
    if not os.path.exists(distribution_outdir):
        cmd = "mkdir -p %s" % distribution_outdir
        os.system(cmd)
    

    count = 0
    m = 2
    some_list = []
    while (count < number_of_chunks):
        some_list.append(m)
        print ( m )
        m = m + 3
        count = count + 1

    def hasNumbers(st_in):
        return any(char.isdigit() for char in st_in)

    os.chdir(distribution_outdir)
    for key, value in energiesDict.items():
        if key == "interresidues" or key == "intraresidues":
            continue
        elif hasNumbers(key) == False: 
            continue
        else:
            new_joined_list = []
            n = 0
            while (n < len (some_list)):
                new_joined_list.append(value[some_list[n]])
                n = n + 1
            new_joined_list = [j for i in new_joined_list for j in i]
            x = [float(i) for i in new_joined_list]
            y = np.array(x)
            (n, bins) = np.histogram(y, bins = 100, range = (-2 , 2), normed=True)
            n = n / (sum(n))
            bincentre = 0.5*(bins[1:]+bins[:-1])
            index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
# Add to empty bins only - amount added is equal to 1/(number empty bins) - to allow calculation of KL
            total_bin_addition = 0.000001
            all_bins = len(bincentre)
# To count the number of populated and non populated bins, to allow dividion of the total bin addition
            non_zero = np.count_nonzero(n)
#            import pdb ; pdb.set_trace()
            zero_bins = all_bins - non_zero
           # import pdb ; pdb.set_trace()
            bin_addition = total_bin_addition/float(zero_bins)
            for i in range(len(n)):
                if n[i]==0.0:
                    n[i] = bin_addition
            data = np.vstack((index, n)).T
            filename = re.sub('[^0-9]', '', key)
            filename = ("residue_%s" % (filename))
            np.savetxt(filename, data, fmt=['%d', '%.20f'])


    os.chdir('..')


# original collect-energies.py does not seem to work with a list in the dictionary, so need to add this here and edit to remove the list from the dictionary and complete the analysis. 
