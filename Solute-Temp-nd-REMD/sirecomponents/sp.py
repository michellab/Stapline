
import sys
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
import numpy

combining_rules = "arithmetic"

def createSystem(molecules):
    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")

    solvent = MoleculeGroup("solvent")
    solute = MoleculeGroup("solute")

    for molecule in moleculeList:
        resname =  molecule.residues()[0].name()
        #print(resname)
        if ( resname == ResName("WAT") ):
            solvent.add(molecule)
            #print('wat')
        elif ( resname == ResName("TFE") ):
            solvent.add(molecule)
            #print('tfe')
        elif ( resname == ResName('Cl-') ):
            #solvent.add(molecule)
            print('ion')
        elif ( resname ==  ResName('Na+') ):
            #solvent.add(molecule)
            print('ion')
        else :
            solute.add(molecule)
            #print(resname)
        #molecules.add(molecule)

    all = MoleculeGroup("all")
    #all.add(molecules)
    all.add(solvent)
    all.add(solute)

    # Add these groups to the System
    system = System()

    system.add(all)
    #system.add(molecules)
    system.add(solvent)
    system.add(solute)

    return system

def setupForcefields(system, space):

    print("Creating force fields... ")

    all = system[MGName("all")]
    #molecules = system[MGName("molecules")]
    solvent = system[MGName("solvent")]
    solute = system[MGName("solute")]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedffsolvent = InterCLJFF("solvent:solvent")
    internonbondedffsolvent.add(solvent)

    # Now solvent intramolecular CLJ energy
    intranonbondedffsolvent = IntraCLJFF("solvent-intranonbonded")
    intranonbondedffsolvent.add(solvent)


    # Now solvent bond, angle, dihedral energy
    intrabondedffsolvent = InternalFF("solvent-intrabonded")
    intrabondedffsolvent.add(solvent)


    # Now solute bond, angle, dihedral energy
    intrabondedffsolute = InternalFF("solute-intrabonded")
    intrabondedffsolute.add(solute)

    # Now solute intramolecular CLJ energy
    intranonbondedffsolute = IntraCLJFF("solute-intranonbonded")
    intranonbondedffsolute.add(solute)

    # Now solute - solvent intermolecular energy
    solutesolventinternonbondedff = InterGroupCLJFF("solute:solvent")
    solutesolventinternonbondedff.add(solute, MGIdx(0))
    solutesolventinternonbondedff.add(solvent, MGIdx(1))

    # Here is the list of all forcefields
    forcefields = [internonbondedffsolvent, intrabondedffsolvent, intrabondedffsolute, intranonbondedffsolute, solutesolventinternonbondedff, intranonbondedffsolvent]

    print(forcefields)
    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("combiningRules", VariantProperty(combining_rules))

    total_nrg = internonbondedffsolvent.components().total() +  intrabondedffsolvent.components().total() + \
                intranonbondedffsolute.components().total() + intrabondedffsolute.components().total() + \
                solutesolventinternonbondedff.components().total()

    e_total = system.totalComponent()
    system.setComponent(e_total, total_nrg)

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add("total_energy", MonitorComponent(e_total, Average()))
    return system



top_file = sys.argv[1]

rst_file = sys.argv[2]

amber = Amber()
print(rst_file ,top_file )
molecules, space = amber.readCrdTop(rst_file, top_file)

system = createSystem(molecules)
system = setupForcefields(system, space)

print (system.energy())

energysymbols = system.energySymbols()

energysymbols =numpy.array(energysymbols)
energysymbols = numpy.sort(energysymbols)
for symbol in energysymbols:
    print ( symbol,system.energy(symbol))
