#!/usr/bin/python
#
# 
# Julien Michel.
# 
# Sep 2016 - update to mdtraj
#
# Jan 2013. 
#
# Compute average interaction energy components of a residue with the rest of the simulated system
#
# TODO
# 
# write average, stddev, number of measurements to pickle to reduce output file size
# 
# Usage 
# 
# To compute the interaction energies of residue number 143
#
# script SYSTEM.top SYSTEM.crd traj.dcd 143
#
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *

import os, sys, pickle
import re

#from MDAnalysis import Universe

try:
    import mdtraj
except ImportError:
    print ("This script depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
    sys.exit(-1)

######  CALCULATION PARAMETERS ##
combining_rules = "arithmetic"
temperature = 298 * kelvin
cutoff = 10.0 * angstrom # Using a reaction field
## the dielectric for the reaction field
rfdielectric=78.3
#####
SOLVENT_RESNAMES = ["WAT","T3P","HOH","T4P"]
IONS_RESNAMES = ["Na+","Cl-"]


def createSystem(molecules, iresnum):

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    all = MoleculeGroup("all")

    imolnum = 0

    for molecule in moleculeList[0:]:
        residues = molecule.residues()
        for x in range(0,residues.count()):
            residue = residues.at(x)
            if residue.number().value() == iresnum:
                imolnum = molecule.number().value()
        
        all.add(molecule)

    # Add these groups to the System
    system = System()

    system.add(all)

    if imolnum == 0:
        print ("FATAL. Could not find a molecule containing residue number %s in the loaded system ! " % (ires))
        sys.exit(-1)

    return system, imolnum

def setupForcefields(system, space, iresnum, imolnum):

    # First we break the system into groups per residues
    # Except for solvent (single group) 
    
    all = system[ MGName("all") ]
    molecules = all.molecules()
    molnums = molecules.molNums()

    iresgrp = MoleculeGroup("ires")
    solvent = MoleculeGroup("solvent")

    samemol_other_residues = []
    othermols_residues = []

    # Dictionary to store interaction energies of ires with everything else
    energiesDict = {}
   
   # import pdb ; pdb.set_trace()

    for molnum in molnums:
        mol = molecules.at(molnum).molecule()
        residues = mol.residues()
        for x in range(0,residues.count()):
            residue = residues.at(x)

            if residue.number().value() == iresnum:
                # IRES
                #print " residue %s is IRES " % residue.toString()
                iresgrp.add( PartialMolecule(mol, residue.selection() ) )
            elif molnum.value() == imolnum:
                # Different residue in same molecule
                #print "residue %s is in same mol as IRES " % residue.toString()
                grp = MoleculeGroup("res%s" % residue.number().value())
                grp.add(  PartialMolecule(mol, residue.selection() ) )
                samemol_other_residues.append(grp) 
            elif ( residue.name().value() in SOLVENT_RESNAMES or 
                   residue.name().value() in IONS_RESNAMES):
                #print "residue %s is a water/ion molecule " % residue.toString()
                solvent.add(  PartialMolecule(mol, residue.selection() ) )
            else:
                #print "residue %s is a solute/protein residue in another molecule" % residue.toString()
                grp = MoleculeGroup("res%s" % residue.number().value())
                grp.add(  PartialMolecule(mol, residue.selection() ) )                
                othermols_residues.append(grp) 

    system.add(iresgrp)
    system.add(solvent)

    for grp in samemol_other_residues:
        system.add(grp)

    for grp in othermols_residues:
        system.add(grp)

    # Now setup the force fields

    # First the ires and solvent ffs

# LP - added empty list to each key to append energies at each timestep 
    # ires InternalFF
    ires_internalff = InternalFF("ires_internalff")
    ires_internalff.add( iresgrp )
    system.add( ires_internalff )
    energiesDict["E_{ires_internalff}^{internal}"] = [ 0.0, 0, [] ]

    # ires IntraMolecularCLJ
    ires_intraff = IntraCLJFF("ires_intraff")
    ires_intraff.add( iresgrp )
    system.add( ires_intraff )
    energiesDict["E_{ires_intraff}^{CLJ}"] = [ 0.0, 0, [] ]

    # ires - solvent
    ires_solventff = InterGroupCLJFF("ires:solvent")
    ires_solventff.add( iresgrp, MGIdx(0) )
    ires_solventff.add( solvent, MGIdx(1) )
    system.add( ires_solventff )


    energiesDict["E_{ires:solvent}^{coulomb}"] = [0.0, 0, []]
    energiesDict["E_{ires:solvent}^{LJ}"] = [0.0, 0, []]    
    energiesDict["E_{ires:solvent}^{CLJ}"] = [0.0, 0, []]

    total_nrg = ires_internalff.components().total() + ires_intraff.components().total() +\
        ires_solventff.components().total()

    energiesDict["intraresidues"] = []
    energiesDict["interresidues"] = []

    # Then the intra group forcefields
    for grp in samemol_other_residues:
        grpname =  grp.name().value()
        grpnum =  int( ("%s" % grpname).strip("res") )
        energiesDict["intraresidues"].append(grpnum)
        ires_jres_intraff = IntraGroupCLJFF("intra_ires:%s" % grpname )
        ires_jres_intraff.add( iresgrp, MGIdx(0) )
        ires_jres_intraff.add( grp, MGIdx(1) )
        system.add( ires_jres_intraff )
        energiesDict["E_{intra_ires:%s}^{LJ}" % grpname] = [ 0.0, 0, []]
        energiesDict["E_{intra_ires:%s}^{coulomb}" % grpname] = [0.0, 0, []]   
        energiesDict["E_{intra_ires:%s}^{CLJ}" % grpname] = [0.0, 0, []]   
        total_nrg += ires_jres_intraff.components().total()

    # Then the intergroup forcefields
    for grp in othermols_residues:
        grpname =  grp.name().value()
        grpnum = int( ("%s" % grpname).strip("res") )
        energiesDict["interresidues"].append(grpnum)
        ires_jres_interff = InterGroupCLJFF("inter_ires:%s" % grpname )
        ires_jres_interff.add( iresgrp, MGIdx(0) )
        ires_jres_interff.add( grp, MGIdx(1) )
        system.add( ires_jres_interff )
        energiesDict["E_{inter_ires:%s}^{LJ}" % grpname] = [0.0, 0, []]
        energiesDict["E_{inter_ires:%s}^{coulomb}" % grpname] = [ 0.0, 0, []]
        energiesDict["E_{inter_ires:%s}^{CLJ}" % grpname] = [ 0.0, 0, []]
        total_nrg += ires_jres_interff.components().total()

    system.setProperty( "space", space )
    system.setProperty( "switchingFunction", 
                        CHARMMSwitchingFunction(cutoff) )
    system.setProperty( "useReactionField", VariantProperty(True) )
    system.setProperty( "reactionFieldDielectric", VariantProperty(rfdielectric) )
    system.setProperty( "combiningRules", VariantProperty(combining_rules) )

    e_total = system.totalComponent()
    system.setComponent( e_total, total_nrg )

    return system, energiesDict

def updateSystemfromTraj(system, frame_xyz, cell_lengths, cell_angles):
    traj_coordinates = frame_xyz[0]

    traj_box_x = cell_lengths[0][0].tolist()
    traj_box_y = cell_lengths[0][1].tolist()
    traj_box_z = cell_lengths[0][2].tolist()

    traj_natoms = len(traj_coordinates)

    # Sire does not support non rectangular boxes
    newmols_coords = {}

    traj_index = 0
    mol_index = 0

    molnums = system.molNums()
    molnums.sort()

    for molnum in molnums:
        mol = system.molecule(molnum).molecule()
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()
        # Create an empty coord group using molecule so we get the correct layout
        newmol_coords = AtomCoords( mol.property("coordinates") )
        for x in range(0,molnatoms):
            tmparray = traj_coordinates[traj_index]
            atom_coord = Vector( tmparray[0].tolist() , tmparray[1].tolist() , tmparray[2].tolist() )
            atom = molatoms[x]
            cgatomidx = atom.cgAtomIdx()
            newmol_coords.set( cgatomidx, atom_coord)
            traj_index += 1
        newmols_coords[molnum] = newmol_coords
        mol_index += 1

    if traj_natoms != traj_index:
        print ("The number of atoms in the system is not equal to the number of atoms in the trajectory file ! Aborting.")
        sys.exit(-1)

    changedmols = MoleculeGroup("changedmols")
    mol_index = 0
    for molnum in molnums:
        mol = system.molecule(molnum).molecule()
        newmol_coords = newmols_coords[molnum]
        mol = mol.edit().setProperty("coordinates", newmol_coords).commit()
        changedmols.add(mol)
    system.update(changedmols)

    space = PeriodicBox(Vector( traj_box_x, traj_box_y, traj_box_z ) )
    system.setProperty("space",space)

    return system

if __name__ == "__main__":
    try:
        top_file = sys.argv[1]
        crd_file = sys.argv[2]
        dcd_file = sys.argv[3]
        residue_number = sys.argv[4]
        residue_number = int( residue_number )
    except IndexError:
        print ("Usage is %s top_file crd_file dcd_file residue_number [start_frame] [end_frame] [step_frame]" % ( sys.argv[0] ))
        sys.exit(-1)

    start_frame = -1
    end_frame = 1e10
    step_frame = 1 

    try:
        start_frame = sys.argv[5]
        end_frame = sys.argv[6]
        step_frame = sys.argv[7]
        start_frame = int(start_frame)
        end_frame = int(end_frame)
        step_frame = int(step_frame)
    except IndexError:
        pass

    # Create a Sire system
    amber = Amber()
    molecules, space = amber.readCrdTop(crd_file, top_file)
    system, molecule_number = createSystem(molecules, residue_number)        
    # Define forcefields
    system, energiesDict = setupForcefields(system, space, residue_number, molecule_number)
    
    energy_components = energiesDict.keys()

    # Create an output folder 
    outdir = "energies-%s"  % residue_number
    if not os.path.exists(outdir):
        cmd = "mkdir -p %s" % outdir
        os.system(cmd)

    # Load DCD trajectory
    mdtraj_trajfile = mdtraj.open(dcd_file,'r')
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    mdtraj_trajfile.seek(start_frame)
    current_frame = start_frame

    while (current_frame <= end_frame):
        print ("Processing frame %s " % current_frame)
        print ("CURRENT POSITION %s " % mdtraj_trajfile.tell() )
        frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
        # ...Update coordinates
        system = updateSystemfromTraj(system, frames_xyz, cell_lengths, cell_angles)
        # Compute energy
        system.energy()
        #print system.energy()

# LP: note: energy_components is the keys from energiesDict: energy_components = energiesDict.keys()
        
        # ...Accumulate results
        for component in energy_components:
        # above is "for key in keys"
            if component == "intraresidues" or component == "interresidues":
                continue
            nrg = system.energy( Symbol( component ) ).value()
            # LP append energy computed to the empty list which item [2] in each key
            energiesDict[component][2].append(nrg) 
            # print component, nrg
            energiesDict[component][1] += 1
            energiesDict[component][0] =  energiesDict[component][0] + \
                ( ( nrg -  energiesDict[component][0] )  /   energiesDict[component][1] )
        current_frame += step_frame
        mdtraj_trajfile.seek(current_frame)	

    # Dump energies in output folder
    outpath = os.path.join(outdir,"chunk-%s-to-%s" % (start_frame, end_frame) )
    outfile = open(outpath,'wb')
    pickle.dump(energiesDict, outfile)
    outfile.close()




