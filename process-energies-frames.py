#!/usr/bin/python

import os,sys, math, time

try:
    import mdtraj
except ImportError:
    print ("This script depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
    sys.exit(-1)


top = "SYSTEM.top"
crd = "SYSTEM.crd"
dcd = "TRAJ.dcd"

def usage():
    print ("USage is script RESIDUE nworkers local/queue [maxframes] [stepframes]")
    sys.exit(-1)

if __name__ == "__main__":

    try:
        residue = sys.argv[1]
        nworkers = sys.argv[2]
        flag = sys.argv[3]
    except IndexError:
        usage()

    try:
        residue = int(residue)
        nworkers = int(nworkers)
    except:
        usage()

    if flag != "local" and flag != "queue":
        usage()

    maxframes = -1
    try:
        maxframes = sys.argv[4]
    except:
        pass

    step_frame = 1
    try:
        step_frame = sys.argv[5]
    except:
        pass
    step_frame = int(step_frame)

    maxframes = int(maxframes)

    dcd_traj = mdtraj.open(dcd,'r')
    dcdframes = len(dcd_traj)

    nframes = dcdframes

    if maxframes > 0:
        if maxframes < nframes:
            nframes = maxframes
     
    print ("Will process first %s frames out of %s, stepping every %s frame at a time. We have %s workers" % ( nframes, dcdframes, step_frame, nworkers))

    chunks = []

    ratio = nframes/float(nworkers)
    
    ratio_integer = int( math.floor(ratio) )
    ratio_remainder = ratio - ratio_integer
    
    # Add one more frame every extra steps
    if ratio_remainder > 0:
        extra = int( round( 1/ratio_remainder) )
    else:
        extra = 0

    start = 0
    end = -1
    for x in range(0, nworkers):
        end += ( ratio_integer  )
        if ( extra > 0 ):
            if ( (x % extra) == 0 ):
                end += 1
        chunks.append( ( start, end ) )
        start = end + 1

    #print chunks


    count = 1
    for (start_frame, end_frame) in chunks:
        sub_script = """#!/bin/bash
#SBATCH -n 1
#SBATCH -o log-%%J.out
#SBATCH -e log-%%J.err
#SBATCH -p serial
#SBATCH --time 48:00:00 
# Switch to current working directory

source /etc/profile.d/module.sh
module load openmm/6.3
module load cuda/7.5
module load sire/16.1.0_no_avx

export OMP_NUM_THREADS=1

# Sleep a random value between 1s to 10 minutes
RAND=$((1 + RANDOM %% 600))
echo 'sleep for ' $RAND
sleep $RAND
date
python energy-decomposition.py %s %s %s %s %s %s %s
date

""" % (top, crd, dcd, residue, start_frame, end_frame, step_frame)
        stream = open("nrg_chunk_%05d.sh" % count,"w")
        stream.writelines(sub_script)
        stream.close()
        count += 1

    # Time for the OS to write all files...
    time.sleep(2)
    
    count = 1
    for (start_frame, end_frame) in chunks:    
        if flag == "queue":
            cmd = "sbatch nrg_chunk_%05d.sh" % count
        elif flag == "local":
            cmd = "bash nrg_chunk_%05d.sh 1> /dev/null 2> /dev/null  &" % count
        print (cmd)
        os.system(cmd)
        count += 1
        time.sleep(1)
        #sys.exit(-1)
