#! home/marie/Utilities/sire.app/bin/python
import os, sys, pickle,re#,argparse
import mdtraj
import math
import numpy
import cmath

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

from Sire.Tools import Parameter, resolveParameters, readParams



######  CALCULATION PARAMETERS ##
temperature = 298 * kelvin
## the dielectric for the reaction field
rfdielectric=78.3
kb = 0.0019872041 # Boltzmann constant kcal/mol K
#############################################
SOLVENT_RESNAMES = ["CYC","ZBT","WAT","T3P","HOH","T4P"]
IONS_RESNAMES = ["Na+","Cl-"]
################################################
shift_delta = Parameter("shift delta", 2.0,
                        """Value of the Lennard-Jones softcore parameter.""")

coulomb_power = Parameter("coulomb power", 0,
                          """Value of the Coulombic softcore parameter.""")

combining_rules = Parameter("combining rules", "arithmetic",
                            """Combining rules to use for the non-bonded interactions.""")

def createSystem(molecules):
    #print("Applying flexibility and zmatrix templates...")
    #print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system

def setupForcefields(system, space):

    #print("Creating force fields... ")

    all = system[MGName("all")]
    molecules = system[MGName("molecules")]
    ions = system[MGName("ions")]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedff = InterCLJFF("molecules:molecules")
    internonbondedff.add(molecules)

    #inter_ions_nonbondedff = InterCLJFF("ions:ions")
    #if (cutoff_type.val != "nocutoff"):
    #    inter_ions_nonbondedff.setUseReactionField(True)
    #    inter_ions_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)

#    inter_ions_nonbondedff.add(ions)

    #inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")
    #if (cutoff_type.val != "nocutoff"):
    #    inter_ions_molecules_nonbondedff.setUseReactionField(True)
    #    inter_ions_molecules_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)#

#    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0))
#    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1))

    # Now solute bond, angle, dihedral energy
    intrabondedff = InternalFF("molecules-intrabonded")
    intrabondedff.add(molecules)

    # Now solute intramolecular CLJ energy
    intranonbondedff = IntraCLJFF("molecules-intranonbonded")

    intranonbondedff.add(molecules)


    # Here is the list of all forcefields
    forcefields = [internonbondedff, intrabondedff, intranonbondedff,]
    #print(forcefields)
    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))

    total_nrg = internonbondedff.components().total() + \
                intranonbondedff.components().total() + intrabondedff.components().total()

    e_total = system.totalComponent()
    '''
    e_internonbondedff = system.totalComponent()
    e_intranonbondedff = system.totalComponent()
    e_intrabondedff = system.totalComponent()

    system.setComponent(e_internonbondedff,internonbondedff.components().total())
    system.setComponent(e_intranonbondedff,intranonbondedff.components().total())
    system.setComponent(e_intrabondedff,intrabondedff.components().total())
	'''
    system.setComponent(e_total, total_nrg)

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add("total_energy", MonitorComponent(e_total, Average()))
#    system.add("internonbondedff_energy", MonitorComponent(e_internonbondedff, Average()))
#    system.add("intranonbondedff_energy", MonitorComponent(e_intranonbondedff , Average()))
#    system.add("intrabondedff_energy", MonitorComponent(e_intrabondedff , Average()))
    return system






def singlepoint(topol, trajs):
 energies=[]
 for traj in trajs :
    top_file = topol
    #traj = sys.argv[2] #collection of mdcrd/coordiantes

    #load each frame and use it as a coordinate
    #mdtraj_top = mdtraj.load_prmtop(top_file)
    mdtraj_dcdfile = mdtraj.load_mdcrd(traj,top=top_file)
    nframes= len(mdtraj_dcdfile)
    #create a folder to store the crds
    if not os.path.exists("rst7_sire_files"):
        os.makedirs("rst7_sire_files")



    #print("reading the mdcrd file... %s frames" %(nframes))
    for framenumber in range(0, nframes):
        #create a Sire system
        rst_file = "rst7_sire_files/%i.rst7" % framenumber
        mdtraj_dcdfile[framenumber].save_amberrst7("rst7_sire_files/%i.rst7" %(framenumber))
        amber = Amber()
        molecules, space = amber.readCrdTop(rst_file, top_file)
        system=createSystem(molecules)
        # Define forcefields
        system = setupForcefields(system, space)
        #print(framenumber , system.energy())#.value())
        energies.append(system.energy().value())

    #print("removing folder...")
    cmd = "rm -r rst7_sire_files"
    os.system(cmd)

    minimum = min(energies)
    new_energies =[]
    for val in energies:
        print(val)
        new_val = val - minimum
        # if new_val <1000 : new_energies.append(new_val)
        new_energies.append(new_val)
    #print(len(new_energies))
    outputenergy=open('energySinglepoint.dat'  ,'w')
    for energy in new_energies :     outputenergy.write(str(energy) + '\n')
    return [energy for energy in new_energies ]

if __name__ == "__main__":

    trajs =  sys.argv[2:]
    topol = sys.argv[1] #topology file
    singlepoint(topol, trajs)
