#!/usr/bin/python
#
# TODO 
# - compute standard deviation of averages
#

import os,sys, pickle
from math import sqrt
import re
import numpy as np
import sys

def properties( list ):

    # Average
    avg = 0.0
    n = len(list)
    totalw = 0

    A = 0.0
    W = 0.0
    Q = 0.0

    for x in range(0 , n , 3):
        value = list[x]
        weight = list[x + 1]
        W += weight

        newA =  A + ( weight / W ) * ( value - A)
        newQ = Q + weight *( value - A ) * ( value - newA )

        A = newA
        Q = newQ

    avg = A
    dev = sqrt( Q / W )

    #print avg, dev

    return avg, dev



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

    max_min_list = []

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
            max_min_list.append(max(x))
            max_min_list.append(min(x))

    
    print ("Max value is: ", max(max_min_list))
    print ("Min value is: ", min(max_min_list))
    
    max_bin = 500.0
    min_bin = -500.0

    print ("Using bin range", str(min_bin), "to",str(max_bin), "Type NO if this is not suitable based on max/min values above. Any other key to continue.")
    if input() =="NO":
        print ("Edit range of histogram.")
        max_bin = int(input("Max bin: "))
        min_bin = int(input("Min bin: "))
    
    print ("Now using bin range " + str(min_bin) + " to " + str(max_bin))

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
            (n, bins) = np.histogram(y, bins = 100, range = (min_bin , max_bin), normed=True)
            n = n / (sum(n))
            bincentre = 0.5*(bins[1:]+bins[:-1])
            index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
# Add to empty bins only - amount added is equal to 1/(number empty bins) - to allow calculation of KL
            total_bin_addition = 0.000001
            all_bins = len(bincentre)
# To count the number of populated and non populated bins, to allow dividion of the total bin addition
            non_zero = np.count_nonzero(n)
            zero_bins = all_bins - non_zero
            bin_addition = total_bin_addition/float(zero_bins)
            for i in range(len(n)):
                if n[i]==0.0:
                    n[i] = bin_addition
            data = np.vstack((index, n)).T
            filename = re.sub('[^0-9]', '', key)
            filename = ("residue_%s_%s" % (filename , key))
        #    np.savetxt(filename, data, fmt=['%d', '%.20f'])

    os.chdir('..')

    residue_number = infolder.split("-")[1].strip("/")

    # Internal energy
    Eires_internal, Eires_internal_dev = properties( energiesDict["E_{ires_internalff}^{internal}"] )
    # Intraresidue non-bonded
    Eires_intraclj,  Eires_intraclj_dev = properties( energiesDict["E_{ires_intraff}^{CLJ}"] )

    Eires_intra = Eires_internal + Eires_intraclj
    Eires_intra_dev = sqrt (  Eires_internal_dev**2 +  Eires_intraclj_dev**2 )

    # ires solvent 
    Ecoulires_solvent,  Ecoulires_solvent_dev = properties( energiesDict["E_{ires:solvent}^{coulomb}"] )
    Eljires_solvent,  Eljires_solvent_dev = properties( energiesDict["E_{ires:solvent}^{LJ}"] )
    Ecljires_solvent,  Ecljires_solvent_dev = properties( energiesDict["E_{ires:solvent}^{CLJ}"] )

    # Intra-residues
    intraresidues = energiesDict["intraresidues"]
    intraresidues.sort()
    Eintrares = {}
    for residue in intraresidues:
        symbol = "E_{intra_ires:res%s}^{coulomb}" % residue
        coul_nrg, coul_nrg_dev = properties( energiesDict[symbol] )
        symbol = "E_{intra_ires:res%s}^{LJ}" % residue
        lj_nrg, lj_nrg_dev = properties( energiesDict[symbol] )
        symbol = "E_{intra_ires:res%s}^{CLJ}" % residue
        clj_nrg, clj_nrg_dev = properties( energiesDict[symbol] )
        #clj_nrg = coul_nrg + lj_nrg
        #clj_nrg_dev = sqrt ( coul_nrg_dev**2 + lj_nrg_dev**2 )
        #print residue, clj_nrg, coul_nrg, lj_nrg
        Eintrares[residue] = ( [ clj_nrg, clj_nrg_dev ], [ coul_nrg, coul_nrg_dev], [ lj_nrg, lj_nrg_dev ] )

    # Inter-residues
    interresidues = energiesDict["interresidues"]
    interresidues.sort()
    Einterres = {}
    for residue in interresidues:
        symbol = "E_{inter_ires:res%s}^{coulomb}" % residue
        coul_nrg, coul_nrg_dev = properties( energiesDict[symbol] )
        symbol = "E_{inter_ires:res%s}^{LJ}" % residue
        lj_nrg, lj_nrg_dev = properties( energiesDict[symbol] )
        symbol = "E_{inter_ires:res%s}^{CLJ}" % residue
        clj_nrg, clj_nrg_dev = properties( energiesDict[symbol] )
        #clj_nrg = coul_nrg + lj_nrg
        #clj_nrg_dev = sqrt ( coul_nrg_dev**2 + lj_nrg_dev**2 )
        #print residue, clj_nrg, coul_nrg, lj_nrg
        Einterres[residue] = ( [ clj_nrg, clj_nrg_dev ], [ coul_nrg, coul_nrg_dev], [ lj_nrg, lj_nrg_dev ] )

    orig_stdout = sys.stdout
    f = open("output.txt", 'w')
    sys.stdout = f   
    print ("#EnergyDecompositionAnalysisforresidue %s (in_kcal/mol)" % residue_number)
    print ("#")
    print ("#InternalEnergy")
    print ("#*************")
    print ("%s\tTernal\t%8.2f\t+-\t%.2f" % (residue_number, Eires_intra,  Eires_intra_dev))
    print ("#*************")
    print ("#Solvent Energy\t    CLJ\t\t\t\t    Coulomb\t\t\t    LJ")
    print ("#*************")
    print ("%s\tSolvent\t%8.2f\t+-\t%.2f\t%8.2f\t+-\t%.2f\t%8.2f\t+-\t%.2f" % (residue_number, Ecljires_solvent, Ecljires_solvent_dev,
                                                                              Ecoulires_solvent, Ecoulires_solvent_dev,  Eljires_solvent, Eljires_solvent_dev))
    print ("#*************")
    print ("#Intramolecularenergies\t    CLJ\t\t\t\t    Coulomb\t\t\t    LJ")
    totalintra = 0.0
    totalintradev = 0.0
    for residue in intraresidues:
        print ("%s\t%s\t\t%8.2f\t+-\t%.2f\t%8.2f\t+-\t%.2f\t%8.2f\t+-\t%.2f" % (residue_number, residue,
                                                                               Eintrares[residue][0][0], Eintrares[residue][0][1],
                                                                               Eintrares[residue][1][0], Eintrares[residue][1][1],
                                                                               Eintrares[residue][2][0], Eintrares[residue][2][1] ))
        totalintra += Eintrares[residue][0][0]
        totalintradev = sqrt ( totalintradev**2 +  Eintrares[residue][0][1]**2 )
    print ("#*************")
    print ("%s\tIntra\t%8.2f\t+-\t%.2f" % (residue_number, totalintra, totalintradev))
    print ("#*************")
    print ("#Intermolecularenergies\t    CLJ\t\t\t\t    Coulomb\t\t\t    LJ")
    totalinter = 0.0
    totalinterdev = 0.0
    for residue in interresidues:
        print ("%s\t%s\t\t%8.2f\t+-\t%.2f\t%8.2f\t+-\t%.2f\t%8.2f\t+-\t%.2f" % (residue_number, residue, #
                                                                               Einterres[residue][0][0], Einterres[residue][0][1],
                                                                               Einterres[residue][1][0], Einterres[residue][1][1],
                                                                               Einterres[residue][2][0], Einterres[residue][2][1]))
        totalinter +=  Einterres[residue][0][0]
        totalinterdev =  sqrt ( totalinterdev**2 +  Einterres[residue][0][1]**2 )
    print ("#*************")
    print ("%s\tInter\t%8.2f\t+-\t%.2f" % (residue_number, totalinter, totalinterdev))
    print ("#*************")
    
    sys.stdout = orig_stdout
    f.close() 
