#!/usr/bin/python
#
# TODO 
# - compute standard deviation of averages
#

import os,sys, pickle
from math import sqrt
import re

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

#    import pdb ; pdb.set_trace()


    # Make sure files are ordered by chunks
    for f in files[1:]:
        chunk = os.path.join(infolder,f)
        stream = open(chunk,'rb')
        data = pickle.load(stream)
        stream.close()

       # import pdb ; pdb.set_trace()
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

    # Test using 2 and then 10 workers (10 workers remove comment) - need to find easier way to make "joined_list" 
    os.chdir(distribution_outdir)
    for key, value in energiesDict.items():
        if key == "interresidues" or key == "intraresidues":
            continue
        else:
            joined_lists = (key, value[2]+value[5])  # +value[8]+value[11]+value[14]+value[17]+value[20]+value[23]+value[26]+value[29])
            filename = re.sub('[^0-9]', '', key)
            dist_out = open("residue_%s" % (filename), 'w')
            for item in joined_lists:
                dist_out.write (str(item) + '\n')
            dist_out.close()
    os.chdir('..')
