#! $env python

d= \
"""
===============================================================================
Modification of Amber system topology and building peptide sequences with stapled residu to set up simulations for
Solute Tempering Replica Exchange Simulations method.

Modifies LJ pair potentials and solvent partial charges between F, CF and OF atoms (TFE),  OW (water) and heavy atomtypes
as well as  heavy atomtypes inside the protein


TODO

-scalling ladder : for the moment, takes a list as input or linear make a linear scale

-dihedral scalling for the solvent

-see whether the B factor for the lennard Jones has to be removed for controlling aggregation

------------------------
Examples:
./Solute-temp_v$VERSION.py -seq AAAAA -helical -nreps 6 -smax 0.5 -v
build a pentalaine peptide in an helical conformation and produce 6 replicas with
 REST2 scaling factor down to 0.5;

./Solute-temp_v$VERSION.py -f <prmtop_file> -nreps 6 -smax 0.5 -hmr -v
takes as input amber topology file <prmtop_file> to produce 6 replicas with
REST2 scaling factor down to 0.5; and default HMassRepartition by parmed.tools.

The script relies on ParmEd and Numpy tools.

Marie Bluntzer
s1772078@ed.ac.uk

version 0.2                                                   24/10/2019
================================================================================
"""

import os
import copy
import sys
import argparse
import numpy as np
import parmed
from parmed.constants import NTYPES
import time
import subprocess
import math

def do_amino_acid_dict():
    return {'A':{'code':'ALA','nat':True, 'desc': 'Alanine' } ,\
    'C':{'code':'CYS','nat':True, 'desc': 'Cysteine' } ,\
    'D':{'code':'ASP','nat':True, 'desc': 'Aspartic acid' } ,\
    'E':{'code':'GLU','nat':True, 'desc': 'Glutamic acid'},\
    'F':{'code':'PHE','nat':True, 'desc': 'Phenylalanine'},\
    'G':{'code':'GLY','nat':True, 'desc': 'Glycine'},\
    'H':{'code':'HIS','nat':True, 'desc': 'Histidine'},\
    'I':{'code':'ILE','nat':True, 'desc': 'Isoleucine'},\
    'K':{'code':'LYS','nat':True, 'desc': 'Lysine'},\
    'L':{'code':'LEU','nat':True, 'desc': 'Leucine'},\
    'M':{'code':'MET','nat':True, 'desc': 'Methionine'},\
    'N':{'code':'ASN','nat':True, 'desc': 'Asparagine'} ,\
    'P':{'code':'PRO','nat':True, 'desc': 'Proline'},\
    'Q':{'code':'GLN','nat':True, 'desc': 'Glutamine'},\
    'R':{'code':'ARG','nat':True, 'desc': 'Arginine'},\
    'S':{'code':'SER','nat':True, 'desc': 'Serine'}   ,\
    'T':{'code':'THR','nat':True, 'desc': 'Threonine'},\
    'V':{'code':'VAL','nat':True, 'desc': 'Valine'},\
    'W':{'code':'TRP','nat':True, 'desc': 'Tryptophan'},\
    'Y':{'code':'TYR','nat':True, 'desc': 'Tyrosine'},\
    'a':{'code':'AKS','nat':False, 'desc': '(S)-2-(4-pentenyl)glycine' } ,\
    'b':{'code':'AKR','nat':False, 'desc': '(R)-2-(4-pentenyl)glycine'}, \
    'm':{'code':'A5S','nat':False, 'desc': '(S)-2-(4-pentenyl)alanine' } ,\
    'n':{'code':'A5R','nat':False, 'desc': '(R)-2-(4-pentenyl)alanine' } ,\
    'c':{'code':'LNS','nat':False, 'desc': '(S)-Lysine-N3' }  ,\
    'd':{'code':'LNR','nat':False, 'desc': '(R)-Lysine-N3'}, \
    'e':{'code':'PGS','nat':False, 'desc': '(S)-Propargyglycine'},\
    'f':{'code':'PGR','nat':False, 'desc': '(R)-Propargyglycine'},\
    'g':{'code':'CPA','nat':False, 'desc': 'cystein crossed-linked Paramethylbenzene (part 1)'},\
    'h':{'code':'CPB','nat':False, 'desc': 'cystein crossed-linked Paramethylbenzene (part 2)'},\
    'i':{'code':'COA','nat':False, 'desc': 'cystein crossed-linked Orthomethylbenzene (part 1)'},\
    'j':{'code':'COB','nat':False, 'desc': 'cystein crossed-linked Orthomethylbenzene (part 2)'},\
    'k':{'code':'CMA','nat':False, 'desc': 'cystein crossed-linked Methamethylbenzene (part 1)'},\
	'l':{'code':'CMB','nat':False, 'desc': 'cystein crossed-linked Methamethylbenzene (part 2)'},\
	'w':{'code':'L11','nat':False, 'desc': 'lactam bridge (part 1)'},\
	'x':{'code':'L12','nat':False, 'desc': 'lactam bridge (part 2)'},\
	'o':{'code':'L21','nat':False, 'desc': 'lactam bridge (part 1)'},\
	'p':{'code':'L22','nat':False, 'desc': 'lactam bridge (part 2)'},\
	'q':{'code':'C11','nat':False, 'desc': 'lactam bridge (part 1)'},\
	'r':{'code':'C12','nat':False, 'desc': 'lactam bridge (part 2)'},\
	's':{'code':'C21','nat':False, 'desc': 'lactam bridge (part 1)'},\
	't':{'code':'C22','nat':False, 'desc': 'lactam bridge (part 2)'},\
    'u':{'code':'C31','nat':False, 'desc': 'lactam bridge (part 1)'},\
    'v':{'code':'C32','nat':False, 'desc': 'lactam bridge (part 2)'},\
    'y':{'code':'C91','nat':False, 'desc': 'lactam bridge (part 1)'},\
    'z':{'code':'C92','nat':False, 'desc': 'lactam bridge (part 2)'},\
    'B':{'code':'C81','nat':False, 'desc': 'lactam bridge (part 2)'}}


def make_3_letters_sequence(seq,stapled=False, one_letter=False , three_letters=False):
    '''
      Check whether a 1-letter or 3-letters sequence has been given
      if 1-letter convert it to  3-letters adding an Acetyl group and check for errors in inputs.
      if 1-letter the sequence 'AAAAA' is expected for pentaalanine
      and if 3-letters,  'ALA ALA ALA ALA ALA'  is expected
      howerer 'ALAALAALAALAALA'   will be understood is three_letters=True
      and  A A A A A will be understood if one_letter=True
      also output the xtra librairies to be loaded by tleap and amino acids to be bound

    '''
    aa=do_amino_acid_dict()

    indices_stapled_residues=[]
    sequence = ''
    xtralib=''
    if  one_letter==True and  three_letters==True:
        exit('Error : flag one_letter and three_letter cannot be used together!' )

    if  len(seq.split()) == 1 and three_letters==True : # 3-letters 'ALAALAALAALAALA'
         seq=' '.join([seq[i:i+3] for i in range(0, len(a), 3)])   #will fall in the next if

    if len(seq.split()) != 1 and one_letter==False : # 3-letters is guessed
        if seq[0:3] =='ACE': sequence = seq.upper()
        else :  sequence += seq.upper()
        for AAA in sequence.split() :  #spare the 'ACE'
            if AAA not in [   aa[key]['code'] for key in aa.keys()  ]  :  exit( 'Error: three_letters code for %s  was not recognised ' %AAA )
            elif False in [value['nat'] for key, value in aa.items() if value['code'] == AAA ] :
                if stapled :
                    #indices_stapled_residues.append(sequence.split().index(AAA)+1)
                    xtralib+= 'loadoff  %s/%s_stapled.lib \n loadamberparams  %s/%s_stapled.frcmod  \n'  %(FFdir, AAA , FFdir , AAA ) #Might load twice the same
                else : xtralib+= 'loadoff  %s/%s_open.lib \n loadamberparams  %s/%s_open.frcmod  \n'  %(FFdir, AAA , FFdir , AAA )



    if len(seq.split()) != 1 and one_letter==True :  #one_letter with spaces #will fall in the next if
        seq=''.join(seq.split(''))

    if len(seq.split()) ==1 and three_letters==False :

        for a in seq:
            if   aa.get(a, False) ==False : exit( 'Error: one_letter code for %s  was not recognised ' %a )
            else : sequence=' '.join([sequence , aa[a]['code'] ])
    sequence=' '.join([args.Nter ,  sequence,args.Cter])
# Now that the sequen ce has ben formated to tleap time to check ! :
    for position in range(1,len(sequence.split())) :  #spare the 'ACE'
            AAA = sequence.split()[position]

            if False in [value['nat'] for key, value in aa.items() if value['code'] == AAA ] :
                if stapled :
                    indices_stapled_residues.append( position-1)
                    print(AAA)
                    print(position)
                    xtralib+= 'loadoff  %s/%s_stapled.lib \n loadamberparams  %s/%s_stapled.frcmod  \n'  %(FFdir, AAA , FFdir , AAA ) #Might load twice the same
                else : xtralib+= 'loadoff  %s/%s_open.lib \n loadamberparams  %s/%s_open.frcmod  \n'  %(FFdir, AAA , FFdir , AAA )




#    Cterminal = 'C' + aa[args.seq[-1].upper()]
#    sequence=' '.join([sequence ,Cterminal ])

    if verbose: print(' input file will be generated for this sequence %s ' %(sequence ))
    #logfile.write(' input file have  been generated for this sequence %s '  %(sequence ) )
    return sequence , xtralib, indices_stapled_residues


def setupsimulation(seq, stapled=False, FractTFE=False , FFdir=False, one_letter=False , three_letters=False ) :
 '''
 Build an AMBER topology from an sequence
 '''

#shutil.copyfile(args.pdbfile,'inputs/peptide.pdb' )

 sequence , xtra_libs, indices_stapled_residues = make_3_letters_sequence(seq, stapled , one_letter , three_letters)
 print(indices_stapled_residues)

 #xtra_libs=''
 #indices_stapled_residues=''
 bond=''
 while indices_stapled_residues : bond+= 'bond m.%s.CY  m.%s.CY\n' %( indices_stapled_residues.pop()+2 ,indices_stapled_residues.pop() +2)
 print(bond)
 if FractTFE :
     if verbose: print('A mixture of %s/%s TFE/WATERTIP4EW will be used' %( FractTFE, 100 -FractTFE ))
     solvent = 'loadoff %sTFE%sBOX.lib \nloadamberparams %sTFE.frcmod \nsolvatebox m  TFE%sBOX 10.0 iso \n' %(FFdir , FractTFE , FFdir, FractTFE)
 else :
     solvent = 'solvatebox m  TIP4PEWBOX 8 iso \n'
     if verbose: print('WATERTIP4EW will be used' )

 #impose an helical start( up to 20 residues )
 helical=''
 if args.helical : helical = 'impose m { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 } { { "N" "CA" "C" "N" -40.0 } { "C" "N" "CA" "C" -60.0 } }'
 tleap=open('tleap.in' , 'w')
 tleap.writelines(['source leaprc.protein.ff14SB \n' , \
                  'loadamberparams frcmod.ionsjc_tip4pew \n' ,\
                  '%s\n'  %(xtra_libs) ,\
                  'set default PBradii mbondi2  \n'  ,\
                  'source leaprc.water.tip4pew \n' ,\
                  'm = sequence { %s } \n'%(sequence) ,\
                  '%s \n' %(helical) ,\
		          '%s \n' %(bond) ,\
                  'relax m  \n' , \
                  'savepdb m peptide.pdb \n',\
                  '%s '%(solvent)  , \
                  ### addions set to 0  --->  ions will be added to neutralise the system
                  'addIons m Na+ 0 \n' ,\
                  'addIons m Cl- 0  \n' ,\
                  'saveamberparm m system.wat.leap.prmtop system.wat.leap.rst7 \n' , \
                  'quit' ] )
 tleap.close()
 subprocess.call('tleap -f tleap.in' , shell=True)!=0


#----------------------------------------------
def get_protein_atomtypes(input_parm):
    """
    Read all residues and atoms and return a list of sidechain (+ C_alpha)
    C and S atomtypes to be used in STREMD scaling.
    """

    solventype= ['OW','HW', 'CS', 'F', 'OF', 'HO', 'Cl-','Na+', 'K+']
    solvent=['WAT', 'SOL' , 'TFE','Cl-','Na+', 'K+' ]

    # protein atomtypes
    proteic_atomtypes = []
    for resid in input_parm.residues:
        for atom in resid.atoms:

             # C, N , O and S
            if atom.type.isupper() and atom.type not in proteic_atomtypes and resid.name not in solvent:
                #if atom.atomic_number in [ 6, 7, 8, 16]:

                    proteic_atomtypes.append(atom.type)
    if verbose == True : print('protein atom types are : ' + str(proteic_atomtypes))
    return proteic_atomtypes

def get_protein_atomtypes_withoutHydrogen(input_parm, verbose=False):
    """
    Read all residues and atoms and return a list of sidechain (+ C_alpha)
    C and S atomtypes to be used in STREMD scaling.
    """

    solventype= ['OW','HW', 'CS', 'F', 'OF', 'HO', 'Cl-','Na+', 'K+']
    solvent=['WAT', 'SOL' , 'TFE','Cl-','Na+', 'K+' ]

    # protein atomtypes
    proteic_atomtypes = []
    for resid in input_parm.residues:
        for atom in resid.atoms:

             # C, N , O and S
            if atom.type.isupper() and atom.type not in proteic_atomtypes and resid.name not in solvent:
                if atom.atomic_number in [ 6, 7, 8, 16]:

                    proteic_atomtypes.append(atom.type)
    if verbose == True : print('protein atom types WITHOUT hydrogen are : ' + str(proteic_atomtypes))
    return proteic_atomtypes



def get_protein_atom(input_parm):
    """
    Read all residues and atoms and return a list of sidechain (+ C_alpha)
    C and S atomtypes to be used in STREMD scaling.
    """
    protein_atomtypes=get_protein_atomtypes(input_parm)

    # protein atomtypes
    proteic_atom = []
    for resid in input_parm.residues:
        for atom in resid.atoms:
            #if atom.atomic_number in [ 1, 6, 7, 8, 16]: # C, N , O and S
                #if atom.type.isupper() and atom.idx not in proteic_atom and atom.type not in solventype:
                if atom.type in protein_atomtypes and atom.idx not in proteic_atom:
                    proteic_atom.append(atom.idx)
    return proteic_atom

#----------------------------------------------
def do_HMR(input_parm, HM_Da):
    """
    Perform H-Mass Repartition.  ---- Probably ambitious but easy to implement ....
    HM_Da - hydrogen mass in Daltons (float).
    """
    act_HMR = parmed.tools.actions.HMassRepartition(input_parm, HM_Da,'dowater')
    log.write(str(act_HMR)+"\n")
    if verbose: print('Mass Repartioning was used')
    print(str(act_HMR))
    act_HMR.execute()

    #input_parm.parm_comments['USER_COMMENTS'] += [str(act_HMR)]

    return input_parm

#----------------------------------------------
def combineLJPair(rh_i, eps_i, rh_j, eps_j):
    """
    LJ pair interaction Lorentz/Berthelot mixing rules for sigma and epsilon
    """
    Rmin_ij = (rh_i + rh_j)
    eps_ij = (eps_i * eps_j)**(0.5)
    return Rmin_ij, eps_ij


#----------------------------------------------
def scale_solvent_protein_LJPairs(input_parm, atomtypes_changed, sfactor, verbose=False):
    """
    Not USED anymore
    """
    rh_OW = input_parm.LJ_radius[input_parm.LJ_types["OW"]-1] # Rmin/2
    eps_OW = input_parm.LJ_depth[input_parm.LJ_types["OW"]-1]

    # add COMMENTS
    input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % sfactor,\
                                                  'scaled OW LJ pair interactions with atomtypes',', '.join(atomtypes_changed)]
    # FOR GROMACS
    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor

    # select atomtnames for scaling atomtypes
    for atype in atomtypes_changed:
        rh_i = input_parm.LJ_radius[input_parm.LJ_types[atype]-1] # Rmin/2
        eps_i = input_parm.LJ_depth[input_parm.LJ_types[atype]-1] # epsilon
        Rmin_iOW, eps_iOW = combineLJPair(rh_i, eps_i, rh_OW, eps_OW)
        eps_iOW *= sfactor # apply the scaling

        # use atype.replace("*", "\*") in selection to avoid wild-card results!
        # strangely "\*" works instead of "'*'" here, unlike in the ParmEd manual
        # on selection masks - http://parmed.github.io/ParmEd/html/parmed.html#atom-selection-masks
        act_changeLJPair = parmed.tools.actions.changeLJPair(input_parm, "@%" + atype.replace("*", "\\\*"), "@%OW", Rmin_iOW, eps_iOW)
        # write action to log
        if list(act_changeLJPair.mask1.Selected()) == []:
            print("WARNING: empty selection %s !\n" % act_changeLJPair.mask1.mask)
            exit(3)

        act_changeLJPair.execute()

        # add equivalent line in GROMACS topology
        #    sigma = Rmin_iOW * 2**(1.0/6) / 10 # convert to A
        #    epsilon = eps_iOW * 4.184     #convert to kJ/mol
        #    gromacs_nonbond += "%-10s   %-10s  1   %8.6e   %8.6e\n" % ("OW", atype, sigma, epsilon)

    return input_parm, gromacs_nonbond

def scale_solvent_LJ(input_parm,  sfactor, verbose=False):
    """
    Multiplies water and TFE (atomtype = "OW", C O ) LJPair depth (epsilon) with solvent  by sfactor according to REST2 scheme.
    Returns an updated parmed parameter object.
    """


    # add COMMENTS
    input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % sfactor,\
                                                  'scaled OW , CS, F, OF LJ epsilon' ]
    # FOR GROMACS
    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor

    #Define solvent atomype : Water,  Water/TFE , TFE . Maybe have to come with something more elegant ?
    solventype=['OW','HW' , 'EP']


    if FractTFE not in [ 0 , False] :
            solventype= ['OW','HW' 'CS', 'F', 'OF']
    if FractTFE in [ 100] :
                solventype= ['CS', 'F', 'OF']

    for atype in solventype:
        try :
            rh_s = input_parm.LJ_radius[input_parm.LJ_types[atype]-1] # Rmin/2
            eps_s = input_parm.LJ_depth[input_parm.LJ_types[atype]-1]
        except : continue
        eps_s *= sfactor   #apply the factor
        act_changeLJ = parmed.tools.actions.changeLJSingleType(input_parm, "@%" + atype.replace("*", "\\\*"), rh_s, eps_s)
        act_changeLJ.execute()

    return input_parm
'''
            # add equivalent line in GROMACS topology
            sigma = Rmin_iS * 2**(1.0/6) / 10 # convert to A
            epsilon = eps_iS* 4.184     #convert to kJ/mol
            gromacs_nonbond += "%-10s   %-10s  1   %8.6e   %8.6e\n" % (atype, atype, sigma, epsilon)
'''

def scale_protein_LJ(input_parm,  sfactor, verbose=False):
    """
    MultipliesProtein atom type (atomtype = "OW", C O ) LJPair depth (epsilon)   by sfactor according to REST2 scheme.
    Returns an updated parmed parameter object.
    """


    # add COMMENTS

    # FOR GROMACS
    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor


    if verbose : print('Scalling LJ protein using scale_protein_LJ(%s, %s ,%s)' %(input_parm,  sfactor, verbose))
    atypes=get_protein_atomtypes_withoutHydrogen(input_parm)
    #for i  in range(len(atypes),-1, -1):
        #if atypes[i][-1].isdigit() ==True:
            #atypes.pop(i)

        #for atype in atypes :

    if verbose : print(input_parm.LJ_types)
    '''
    lj_types=input_parm.LJ_types
    print(lj_types['N'])
    lj_types_inv = {tuple(v): k for k, v in lj_types.items()}
    lj_types_unique_values={v: list(k) for k, v in lj_types_inv.items()}
    print(lj_types_unique_values)
    '''
    done=[]
    for atype in get_protein_atomtypes(input_parm ) :
        #print(atype)

        atype= formatatype(atype)
        rh_p = input_parm.LJ_radius[input_parm.LJ_types[  atype]-1] # Rmin/2
        eps_p = input_parm.LJ_depth[input_parm.LJ_types[  atype]-1]

        if verbose : print(atype ,  eps_p)
        eps_p *= sfactor   #apply the factor
        #if verbose : print(atype , eps_p)
        if input_parm.LJ_types[  atype] not in done:
                #print("@%" + atype )
                act_changeLJ = parmed.tools.actions.changeLJSingleType(input_parm, "@%" + atype , rh_p, eps_p)
                act_changeLJ.execute()
                done.append(input_parm.LJ_types[  atype])


        #input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % sfactor, 'scaled %s atomtypes epsilon' %get_protein_atomtypes_withoutHydrogen(input_parm) ]


    return input_parm
'''
            # add equivalent line in GROMACS topology
            sigma = Rmin_iS * 2**(1.0/6) / 10 # convert to A
            epsilon = eps_iS* 4.184     #convert to kJ/mol
            gromacs_nonbond += "%-10s   %-10s  1   %8.6e   %8.6e\n" % (atype, atype, sigma, epsilon)
'''

def formatatype(atype):
    j=0

    #for i in atype:
    '''
        if i.isdigit():
            #atype= atype[:j] + atype[j+1:]
            #atype= atype.pop(j)
            #j+=-1
            atype= atype[:j] + atype[j+1:]
            print(atype)
            if atype in get_protein_atomtypes_withoutHydrogen(input_parm, verbose=False):dupl =True
        j+=1
    '''
    #atype=atype.replace("*", "\\\*")
    #print(atype)
        #atype="\"" + atype + "\""
    return atype

#-------------------------------------------------------------------------------------------


def scale_solvent_charges(input_parm,  sfactor, verbose=False):
    """
    Multiplies water and TFE (atomtype = "OW", C O ) LJPair depth (epsilon) with solvent  by sfactor according to REST2 scheme.
    Returns an updated parmed parameter object.
    """


    # add COMMENTS
    #input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % math.sqrt(sfactor),\
                                            #      'scaled OW , CS, F, OF partial charges ' ]


    solventype=[]
    solventype= [ 'EP','HW' ,'OW' ]

    if FractTFE not in [ 0 , False] :
        solventype= ['EP','HW', 'CS', 'F', 'OF' ,'HF', 'HT']
    if FractTFE in [ 100] :
            solventype= [ 'CS', 'F', 'OF' ,'HF', 'HT']

    solventcharge=np.zeros(len(solventype))
    SOL_RES=[]
    for resid in input_parm.residues:
            if  0 not  in solventcharge  : break
            for atom in resid.atoms:
                if atom.type.isupper() and atom.type in solventype:
                    solventcharge [solventype.index(atom.type)] =atom.charge

    solventcharge *= math.sqrt(sfactor)   #apply the factor

    for i in solventype :
    #    try :
            act_changecharge=parmed.tools.actions.change(input_parm, input_parm.charge_flag , "@%" + i,solventcharge[solventype.index(i)] )
            act_changecharge.execute()
    #        except : continue

    solventcharge=np.zeros(len(solventype))
    #while 0 in solventcharge:
    for resid in input_parm.residues:
            for atom in resid.atoms:
                if atom.type.isupper() and atom.type in solventype:
                    solventcharge [solventype.index(atom.type)] =atom.charge

    return input_parm

#-------------------------------------------------------------------------------------------


def scale_protein_charges(input_parm,  sfactor, verbose=False):
    """
    Multiplies every not water or TFE (atomtype = "OW", C O ) LJPair depth (epsilon) with solvent  by sfactor according to REST2 scheme.
    Returns an updated parmed parameter object.
    """
    if verbose: print ( 'Scaling protein charges')

    # add COMMENTS
#    input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % math.sqrt(sfactor),\
                                                #  'scaled OW , CS, F, OF partial charges ' ]
    # FOR GROMACS
    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor


    protype= get_protein_atom(input_parm)
    #if verbose: print (protype)
    #protcharge=np.zeros(len(protype))
    protcharge=np.zeros(len(protype))
    for resid in input_parm.residues:

            if  0 not  in protcharge  : break # shorten the charge search : but not if one protein atom got charge of zero
            for atom in resid.atoms:
                #if verbose: print (protype.index(atom))
                if atom.type.isupper() and atom.idx in protype:

                    protcharge [protype.index(atom.idx)] =atom.charge


    protcharge *= math.sqrt(sfactor)   #apply the factor


    for i in protype :

        #try :
            act_changecharge=parmed.tools.actions.change(input_parm, input_parm.charge_flag , "@" + str(i+1),protcharge[protype.index(i)] )
            act_changecharge.execute()
        #except : continue

    protcharge=np.zeros(len(protype))
    #while 0 in protcharge:
    for resid in input_parm.residues:
            for atom in resid.atoms:
                if atom.type.isupper() and atom.idx  in protype:
                    protcharge [protype.index(atom.idx)] =atom.charge

    return input_parm

#-------------------------------------------------------------------------------------------

def scale_protein_dihedral(input_parm,  sfactor, verbose=False):
    """
    Multiplies every not water or TFE (atomtype = "OW", C O ) LJPair depth (epsilon) with solvent  by sfactor according to REST2 scheme.
    Returns an updated parmed parameter object.
    """

    protype= get_protein_atomtypes(input_parm)
    # add COMMENTS
    #input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % math.sqrt(sfactor),\
                                                  #'scaled dihedrals for atoms in ' +  str(protype)]
    # FOR GROMACS

    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor

    protcharge=np.zeros(len(protype))
    SOL_RES=[]

    protdihtype=[]
    for dih in input_parm.dihedrals:

            if (dih.atom2.type  in  protype)    and  dih.type not in protdihtype: # atom2 to exclude Hydrogens which are not in protype
                protdihtype.append( dih.type)
            else :
                if  dih.type not in protdihtype and   dih.atom2.type  not in ['CS','OF'] : print(dih)
    if verbose : print(str(len(protdihtype))+ ' protein dihedral types to scale on ' +str(len( input_parm.dihedral_types )) + ' factor ' +str(sfactor))
    for dih in protdihtype:
                if verbose : print(dih ,dih.phi_k *sfactor )
                dih.phi_k = dih.phi_k *sfactor


    return input_parm

#-------------------------------------------------------------------------------------------
def neutralize_protein_charges(input_parm,  verbose=True):
    '''
    Neutralize topology charge par repartitioning the protein charge over the couterions
    '''
    #poscounterions=['Na+','K+']
    #negcounterions=['Cl-']
    ## list auto-updated :avoid warning if mask is empty allow all kingd of couterions if they got + or - at the end of theyre atomtype name
    poscounterions=[]
    negcounterions=[]
    count=0
    charge=parmed.tools.actions.netCharge(input_parm).execute()
    if charge < 0 :
            for resid in input_parm.residues:
                for atom in resid.atoms:
                    if atom.type[-1] == '-':
                        count+=1
                        if atom.type not in negcounterions : negcounterions.append(atom.type)
            if verbose : print('%s neg counterions found' %(count))
            for i in negcounterions:
                act_changecharge=parmed.tools.actions.change(input_parm, input_parm.charge_flag , "@%" +i ,(1-1/count*charge) )
                act_changecharge.execute()
    else:
        for resid in input_parm.residues:
                for atom in resid.atoms:
                    if atom.type[-1] == '+':
                        count+=1
                        if atom.type not in poscounterions : poscounterions.append(atom.type)
        for i in poscounterions:
            act_changecharge=parmed.tools.actions.change(input_parm, input_parm.charge_flag , "@%" +i ,(1-1/count*charge) )
            act_changecharge.execute()
        if verbose :print('%s nposcounterions found' %(count)        )


    return input_parm

def scale_solvent_solvent_LJPairs(input_parm, atomtypes_changed, sfactor, verbose=False):
    """
    Multiplies water and TFE (atomtype = "OW", CS O ) LJPair depth (epsilon) with solvent  by sfactor according to REST2 scheme.
    Returns an updated parmed parameter object.   DEAD
    """


    # add COMMENTS
    input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % sfactor,\
                                                  'scaled solvent  LJ pair interactions with atomtypes',\
						  ', '.join(atomtypes_changed)]
    # FOR GROMACS
    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor

    # select atomtnames for scaling solvent = water and TFE
    for atype in ['OW','CS', 'F', 'OF']:
    	for btype in ['OW', 'CS', 'F', 'OF']:
            rh_S1 = input_parm.LJ_radius[input_parm.LJ_types[atype]-1] # Rmin/2
            eps_S1 = input_parm.LJ_depth[input_parm.LJ_types[atype]-1]
            rh_S2 = input_parm.LJ_radius[input_parm.LJ_types[btype]-1] # Rmin/2
            eps_S2 = input_parm.LJ_depth[input_parm.LJ_types[btype]-1] # epsilon
            Rmin_S, eps_S = combineLJPair(rh_S1, eps_S1, rh_S2, eps_S2)
            eps_iS *= sfactor # apply the scaling

            # use atype.replace("*", "\*") in selection to avoid wild-card results!
            # strangely "\*" works instead of "'*'" here, unlike in the ParmEd manual
            # on selection masks - http://parmed.github.io/ParmEd/html/parmed.html#atom-selection-masks
            act_changeLJPair = parmed.tools.actions.changeLJPair(input_parm, "@%" + atype.replace("*", "\\\*"), "@%" + atype.replace("*", "\\\*"), Rmin_S, eps_iS)
            # write action to log
            if list(act_changeLJPair.mask1.Selected()) == []:
                print("WARNING: empty selection %s !\n" % act_changeLJPair.mask1.mask)
            exit(3)

            act_changeLJPair.execute()

            # add equivalent line in GROMACS topology
            sigma = Rmin_iS * 2**(1.0/6) / 10 # convert to A
            epsilon = eps_iS* 4.184     #convert to kJ/mol
            gromacs_nonbond += "%-10s   %-10s  1   %8.6e   %8.6e\n" % (atype, btype, sigma, epsilon)

    return input_parm, gromacs_nonbond


#----------------------------------------------
def scale_protein_protein_LJPairs(input_parm, atomtypes_changed, sfactor, verbose=False):  #REST2
    """
  DEAD
    """

    # add COMMENTS
    input_parm.parm_comments['USER_COMMENTS'] += ['REST2 METHOD (scaling factor = %.3f)' % sfactor,\
                                                  'scaled P-P LJ pair interactions with atomtypes',\
						  ', '.join(atomtypes_changed)]
    # FOR GROMACS
    gromacs_nonbond = "; REST2 scaling - %.2f\n" % sfactor

    # select atomtnames for scaling atomtypes
    for k in range(len(atomtypes_changed)):
        atype=atomtypes_changed[k]
        for j in range(k,len(atomtypes_changed)):
            btype=atomtypes_changed[j]

            rh_k = input_parm.LJ_radius[input_parm.LJ_types[atype]-1] # Rmin/2
            eps_k = input_parm.LJ_depth[input_parm.LJ_types[atype]-1] # epsilon
            rh_j = input_parm.LJ_radius[input_parm.LJ_types[atype]-1] # Rmin/2
            eps_j = input_parm.LJ_depth[input_parm.LJ_types[atype]-1] # epsilon
            Rmin_s, eps_s = combineLJPair(rh_k, eps_k, rh_j, eps_j)
            eps_j *= sfactor # apply the scaling

            act_changeLJPair = parmed.tools.actions.changeLJPair(input_parm, "@%" + atype.replace("*", "\\\*"), "@%"  + atype.replace("*", "\\\*"), Rmin_s, eps_s)
            # write action to log
            if list(act_changeLJPair.mask1.Selected()) == []:
                print("WARNING: empty selection %s !\n" % act_changeLJPair.mask1.mask)
                exit(3)

            act_changeLJPair.execute()

            # add equivalent line in GROMACS topology
            sigma = Rmin_s* 2**(1.0/6) / 10 # convert to A
            epsilon = eps_j * 4.184     #convert to kJ/mol
            gromacs_nonbond += "%-10s   %-10s  1   %8.6e   %8.6e\n" % ("OW", atype, sigma, epsilon)

    return input_parm, gromacs_nonbond


#----------------------------------------------
def get_NONBONDED_PARM_INDEX(input_parm, atype_i, atype_j):
    """
    Get NONBONDED_PARM_INDEX for atomtype_i and atomtype_j LJ interaction.   called by remove_BCOEF
    Returns and index using Amber prmtop file convention (start index = 1).
    If parm_data is modified as Python list - the index needs to be reduced by 1.
    """

    # get ATOM_TYPE_INDEX for both i and j
    i, j = input_parm.LJ_types[atype_i], input_parm.LJ_types[atype_j]
    # max is equal to total value of types (indexed starting with 1)
    #NTYPES = np.max(parm.LJ_types.values())
    ntypes = ntypes = parm.pointers['NTYPES']

    # NONBONDED_PARM_INDEX = NTYPES * [ ATOM_TYPE_INDEX(i) - 1 ] +  ATOM_TYPE_INDEX(j)
    #ij = parm.pointers['NTYPES'] * (i - 1) + j
    ij = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*(i-1)+j-1] - 1

    return ij


#----------------------------------------------
def remove_BCOEF(input_parm, atomtypes):
    """
    Remove BCOEF (C6) for the attractive part of LJ.   #### Might be useful ??? see artcicle from
    Usage eg. to prevent ligand aggregation.
    ----
    input:
        input_parm - default system parameters
        atomtypes - list of atomtypes

    """
    for i in np.arange(len(atomtypes)):
        atype_i = atomtypes[i]
        for j in np.arange(i, len(atomtypes)):
            atype_j = atomtypes[j]

            # get cross interaction coefficient index
            ij = get_NONBONDED_PARM_INDEX(input_parm, atype_i, atype_j)

            # set BCOEF to zero!
            # adjust to pythonic start indexing at 0
            input_parm.parm_data["LENNARD_JONES_BCOEF"][ij - 1] = 0.0

            str_action = "Removing LJ attractive @%"+atype_i+"-@%"+atype_j+" pairwise interaction. (C6, BCOEF = 0.0)"
            log.write(str_action)

    return




#####################################################
###################### MAIN #########################
#####################################################

# keep the time!



if __name__ == "__main__":
    #start_time = time.clock()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
    parser.add_argument("-f", type=str, default='system.wat.leap.prmtop', help='input AMBER topology file name (default: %(default)s)')
    parser.add_argument("-o", type=str, default=None, help='output AMBER topology file name (default: system[HMR]_REST2.XXX.prmtop)')
    parser.add_argument("-nreps", type=int, default=8, help='number of replicas (default: %(default)s)')
    parser.add_argument("-smax", type=float, default=1.0, help="maximum scaling factor (default: %(default)s)")
    parser.add_argument("-smin", type=float, default=0.4, help="minimum scaling factor (default: %(default)s)")
    parser.add_argument("-scale", type=float, nargs='+', default=[], help="list of scaling factor(s); overrides -nreps, -smax, -smin (default: %(default)s)")
    parser.add_argument("-sf", type=str,  default='lin', help="function for scaling : lin/exp (linear/exponetial) (default: %(default)s)")
    parser.add_argument("-ignore", default=False, action='store_true', help="ignore warnings (default: %(default)s)")
    parser.add_argument("-hmr", default=False, action='store_true', help="HMassRepartition (default: %(default)s)")
    parser.add_argument("-hm_Da", type=float, default=None, help="H mass in Daltons used in HMR (default: use parmed default, 3.024)")
    parser.add_argument("-v", "--verbose", action='store_true',default=False,  help="be verbose")
    parser.add_argument("-xyz", type=str, default='system.wat.leap.rst7', help='input AMBER coordinate file name (default: %(default)s)')
    parser.add_argument("-gmx", default=False, action='store_true', help="save GROMACS topology, too (not functional !!!) (default: %(default)s)")
    parser.add_argument('-seq'  ,  type=str , help='input sequence')
    parser.add_argument('-Nter'  ,  type=str ,default='ACE',  help='Nterminal  (default: %(default)s) for none enter  \'\' ')
    parser.add_argument('-Cter'  ,  type=str  ,default='NHE' , help='Nterminal  (default: %(default)s) (amine) for carboxylic enter \'\' ' )

    parser.add_argument('-helical' , action='store_true' ,default='True' , help='whether the peptide start in an helical conformation ' )
    parser.add_argument('-stapled' , action='store_true',default='True' ,  help='whether the peptide is stapled' )
    parser.add_argument('-FracTFE'  ,  default=False,  type=int, help='input sequence')
    parser.add_argument('-FFdir' ,default=os.path.dirname(os.path.realpath(sys.argv[0]))+'/ForceFieldFiles/' , action='store_true',  help='Additional forcefield files for custom residus and solvent  %(default)s' )
    ########################################################
    # PARSING INPUTS
    ########################################################

    args = parser.parse_args()

    inputfile = args.f
    outputname = args.o
    nreps = args.nreps
    smax = args.smax
    smin = args.smin
    scale = [ sf for sf in args.scale ]
    HMR = args.hmr
    HM_Da = args.hm_Da
    verbose = args.verbose

    ignore = args.ignore
    inputcoords = args.xyz
    gmx = args.gmx
    sequence=args.seq
    FractTFE=args.FracTFE
    #FFdir=args.FFdir
    FFdir='/home/marie/Desktop/newparam/libs/'
    stapled=args.stapled
    helical=args.helical
    scale_func=args.sf

    if sequence :
        setupsimulation(sequence,stapled=stapled, FractTFE=FractTFE,  FFdir=FFdir)

    #else :
        #if FractTFE or stapled or helical : exit('Error : input  FractTFE or stapled or helical only valid if paramater files are generated from the file  please provide a sequence')

    #######################
    ##### check the input
    if not os.path.isfile(inputfile):
        print("ERROR: file %s not found or no sequence provided, check your inputs" % inputfile)
        exit(1)
    else:
        input_parm = parmed.amber.readparm.AmberParm(inputfile, xyz=inputcoords)

    print("%s\n SOLUTE TEMPERING SETUP\n%s\n\n" % ("#"*42, "#"*42) +\
          "Input files read!\nStarting the script. This should take about a minute or two...\n")

    if scale == [] and scale_func =='exp'  :
        if ( smax <2 and smax > 0) and ( smin <2 and smin > 0):

            sfactors = np.linspace(smin, smax, nreps)

            for i in range(len(sfactors)):
                sfactors[i] = (1/(math.exp(1)-math.exp(2)))*(math.exp(sfactors[i])-math.exp(2))

        else :
            sys.exit('smax and smin sould be between 1 and 2')
    elif scale == [] and scale_func =='lin'  :
        #smax= 1/smax
        #smin= 1/smin
        sfactors = np.linspace( smax, smin, nreps)
    else:
        sfactors = scale
        if verbose : print(sfactors)
        nreps = len(scale)


    if 1.0 not in sfactors:
        print("WARNING: scaling factors does not include 1.0 ")
        if not ignore:
            print("If above list looks OK, and you know what you are doing: restart with a flag '-ignore' to proceed.")
            exit(2)

    ######################
    ### begin the script

    # define output file root!
    if outputname == None:
        # use the same root as the input file with modifications to be added to the name!
        outputname = ".".join(inputfile.split(".")[:-1])

    dirpath="./"

    logfile = "REST2_%s.log" % time.strftime("%Y-%m-%d", time.gmtime())
    print("Writing all the modifications to the log file: %s\n" % logfile)
    log = open(logfile, 'w')
    log.write("# This log was created on %s \n# by executing this command line:\n# %s\n###\n" %\
              (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), ' '.join(sys.argv)))


    # add USER_COMMENTS flag to .prmtop file!
    #input_parm.add_flag("USER_COMMENTS", "(20a4)", data=None, num_items=0)

    ###########
    #### HMR
    if HMR:
        # perform HMR
        input_parm = do_HMR(input_parm, HM_Da)
        outputname += "+HMR"

    ####Solute Tempering !!

    # use all heavy atomtypes for scaling

    for rep in np.arange(nreps):
        sfactor = sfactors[rep]

        print_replica = '### Replica %03d - scaling factor %.3f\n' % (rep, sfactor)
        log.write(print_replica)
        if verbose: print(print_replica)

        # apply REST2-type scaling!
        input_parm_copy = copy.copy(input_parm)
        #new_parm, REST2_gmx_nonbond = scale_solvent_protein_LJPairs(input_parm_copy, atomtypes_REST2, np.sqrt(sfactor), verbose=True)
        #new_parm, REST2_gmx_nonbond = scale_protein_protein_LJPairs(input_parm_copy, atomtypes_REST2, sfactor, verbose=True)
        if verbose: print('scaling LJ interactions')
        #scale_solvent_LJ  : lorentz berthelot rules : e12 = sqrt(e1 * e2)
        new_parm_LJ = scale_protein_LJ(input_parm_copy, sfactor, verbose)
        if verbose: print('scaling charges interactions')
        #scale_solvent_charges
        new_parm_LJ_charge = scale_protein_charges(new_parm_LJ, sfactor, verbose)
        if parmed.tools.actions.netCharge(new_parm_LJ_charge).execute() !=0 :

            new_parm_charge_reweigthed =neutralize_protein_charges(new_parm_LJ_charge,  verbose)
        else : new_parm_charge_reweigthed = new_parm_LJ_charge
        # check if off-diagonal terms have been changed
        if verbose: print("Off-diagonal terms changed: %s\n" % str(new_parm_LJ_charge.has_NBFIX()))
        if verbose: print('scaling protein dihedral')
        new_parm = scale_protein_dihedral(new_parm_charge_reweigthed, sfactor, verbose)
        # SAVE AMBER .PRMTOP
        new_parm.save(dirpath+'%s+REST2.%03d.prmtop' % (outputname, rep), format="amber", overwrite=True)
        print_saved = '### New topology saved as %s+REST2.%03d.prmtop\n' % (outputname, rep)
        log.write(print_saved)
        if verbose: print(print_saved)

        # SAVE GMX .TOP
        if gmx:
            # GROMACS - modify nonbond
            gmx_nonbond = lig_gmx_nonbond + REST2_gmx_nonbond
            # save
            add_nonbond_params(gmx_output+".top", gmx_nonbond, dirpath+outputname+"+REST2%i.top" % rep)
            print_saved_gmx = '### New GMX topology saved as '+ dirpath+outputname+'+REST2%i.top\n' % rep
            log.write(print_saved_gmx)
            if verbose: print(print_saved_gmx)


    # finish writing the log
    log.close()

    #print("DONE!\nSee logfile for details: %s\n\nRunning time: %.2f s" % (logfile, time.clock() - start_time))
