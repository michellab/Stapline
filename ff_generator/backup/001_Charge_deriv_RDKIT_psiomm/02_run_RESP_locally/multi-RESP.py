from multiprocessing import Pool
import os
import parmed
import numpy as np
import sys
import subprocess
import string
from multiprocessing import Process
import time
import copy
import argparse


def run_resp(input):
    #prepare gaussian files run gaussian jobs and read gaussian files


    #Uncomment the next if statement to rerun gaussian jobs
    if os.path.isfile('./%s.log'%input) == False : # or  'Normal termination' not in  subprocess.getoutput('tail -1 ./%s.log'%input)  :

        	print('running RESP for conformation  %s '%input)
        	os.system('g09 %s.gau' %input  )
        	print('RESP calculation over for %s '%input)
    #time.wait(1)

def sdf2pdb(input):
    # read the sdf (input) and convert it to multiple pdb and return the number of conformations
    tobeparsed= subprocess.getoutput('babel -m -isdf %s -O conf.pdb --split' %(input) )
    n_conf=tobeparsed.split('\n')[-2].split()[0]
    #print(tobeparsed.split('\n')[-2])
    os.system('sed -i  s/HETATM/ATOM\ \ /g  conf*pdb' )
    #os.system('grep -v ATOM  1.pdb > model.pdb')
    return n_conf


def run_allconf(n_conf, n_cpu=10) :
    #call run_resp for all conformation using multiprocessing
    # !!!!! For a RE-RUN  DO NOT USE multiprocessing :
    # 10 instances of antechamber running at the same time seems to break it ! use the run_resp in the 'for' loop

    args=[]

    for i in range(1,n_conf+1):
      args.append('%s' %(i))
      if 'Fatal Error!' in subprocess.getoutput('antechamber -i conf%s.pdb -fi pdb -o %s.gau  -fo gcrt -gv 2 -ge %s.gesp'  %(i, i, i)) : sys.exit('Antechamber couldn\'t convert your file, check your input ')
      os.system('sed -i  s/molecule/conf%s/g  %s.gau' %(i, i) )
      #run_resp('%s' %(i))

    p =Pool(processes=n_cpu)
    r = p.map_async(run_resp,args)
    #print(args)

    r.wait()

    for input in range(1,n_conf+1):
        subprocess.getoutput('antechamber -i  %s.log -fi gout -o %s.mol2 -fo mol2 -c resp -at amber   -eq 2  -pf   -dr n' %(input,input) )


def return_charges(parm, i=0):
    #remove unwanted fragment (any heavy atom fwhich name finishes with a number and bounded Hydrogen)

    frag=[]
    charges=[]

    for atom in  parm.atoms:

      if atom.name[-1].isdigit() ==False and atom.element not in [1] :

        frag.append(atom.idx)
        charges.append(atom.charge)
        for atm in atom.bond_partners:

            if atm.element == 1 :

                frag.append(atm.idx)
                charges.append(atm.charge)
    ####TO DO : here check if similar H had same charge independant
    totcharge= np.sum(charges)
    #print(charges)
    newcharges=np.add( charges ,-totcharge/len(charges) )
    new_parm = recharge(newcharges,frag , parm)
    # could save here for each using
    # new_parm.save(conf%i.mol2 %i)
    return frag, newcharges


def return_charges_full(parm, i=0):
    #os.system('grep -v ATOM  1.pdb > model.pdb')
    #remove unwanted fragment (any heavy atom fwhich name finishes with a number and bounded Hydrogen)

    index=[]
    charges=[]

    for atom in  parm.atoms:

                charges.append(atom.charge)
                index.append(atom.idx)
    ####TO DO : here check if similar H had same charge independant
    totcharge= np.sum(charges)
    #print(charges)
    newcharges=np.add( charges ,-totcharge/len(charges) )
    new_parm = recharge(newcharges,index , parm)
    # could save here for each using
    # new_parm.save(conf%i.mol2 %i)
    return index, newcharges


def multiple(n_conf, to_removefromlib ,resname='MOL'):
    # run return_charges and average all charge for a same atom
    I=[]

    C=[]
    for file in range(1,n_conf+1):
        #print(file)
        ##### TO DO :check index is conserved here :
        parm = parmed.load_file('%s-reo.mol2'%file)
        index, charges =return_charges(parm)
        I.append(index)
        #print(charges)
        if C!=[]  : C= np.concatenate( (C, [charges]) , axis = 0 )
        else :  C=[copy.deepcopy(charges)]

    avg_charge= np.average(C, axis=0)
    #print(avg_charge , index)
    new_parm = recharge(avg_charge,index ,  parm)
    for atm in parm.atoms:
        print(resname)
        print(atm.residue.name)
        atm.residue.name = resname
        print(atm.residue.name)
    parm.save('newcharges.mol2')
    #os.system('sed -i s/MOL/%s/g newcharges.mol2' %(resname))

    makeamberlibfile(index, parm, to_removefromlib ,resname)

def makeamberlibfile(frag, parm ,to_removefromlib, resname='MOL'):
    tleapinput=open('tleap.in' ,'w')
    os.system( 'parmchk2 -i %s  -o  parmchk2_all.frcmod -f  mol2 ' %( 'newcharges.mol2'))
    os.system( "sed '/^DIHE/,/^NONBON/{/^NONBON/!d}'  parmchk2_all.frcmod > parmchk2.frcmod  ")
    tleapinput.write('loadAmberParams parmchk2.frcmod  \n' )
    tleapinput.write('source leaprc.protein.ff14SB \n' )
    tleapinput.write('%s = loadmol2 newcharges.mol2 \n' %(resname) )
    for atom in  parm.atoms:
        if atom.idx not in frag or atom.idx in to_removefromlib :
            tleapinput.write('remove %s %s.1.%s\n' %(resname, resname, atom.idx +1) )
    #for atom  in to_removefromlib :
        #tleapinput.write('remove %s %s.1.%s\n' %(resname, resname, atom) )
    tleapinput.write('set %s head %s.1.N \n' %(resname,resname) )
    tleapinput.write('set %s tail %s.1.C \n' %(resname,resname) )
    tleapinput.write('impose %s { 1 2 3 } { { "N" "CA" "C" "O" 0.0 } { "H" "N" "CA" "C"  0.0 } { "CA" "CB" "CC" "CD" 180 } { "CB" "CC" "CD" "CE" 180 } } \n'%(resname) )
    tleapinput.write('saveoff  %s %s.lib \n'  %(resname,resname) )
    #os.system('grep -v ATOM  1.pdb > model.pdb')

    tleapinput.write('seq = sequence { ACE %s NME }\n' %(resname))
    tleapinput.write('savemol2 seq ACE-%s-NME.mol2 1 \n'%(resname))
    tleapinput.write('loadoff %s/ethyl.lib  \n' %(os.path.dirname(os.path.abspath( __file__ ))))
    tleapinput.write('seq = combine { ETH  seq } \n')
    tleapinput.write('bond seq.1.CY seq.3.CY \n')
    tleapinput.write('select seq.1 \n')
    tleapinput.write('relax seq.1 \n')
    tleapinput.write('select seq \n')
    tleapinput.write('relax seq\n')
    tleapinput.write('savemol2 seq ACE-%s-NME.mol2 1 \n'%(resname))
    tleapinput.write('quit')
    tleapinput.close()
    os.system('tleap -f tleap.in')

def multiple_full(n_conf):
    # run return_charges and average all charge for a same atom
    I=[]

    C=[]
    for file in range(1,n_conf+1):
        #print(file)
        ##### TO DO :check index is conserved here :
        parm = parmed.load_file('%s-reo.mol2'%file)
        index, charges =return_charges_full(parm)
        #print(index)
        #print(charges)
        I.append(index)
        #print(charges)
        if C!=[]  : C= np.concatenate( (C, [charges]) , axis = 0 )
        else :  C=[copy.deepcopy(charges)]

    avg_charge= np.average(C, axis=0)
    #print(avg_charge , index)
    new_parm = recharge(avg_charge,index , parm)
    parm.save('newcharges_full.mol2')

def recharge(newcharges,frag , parm):
    for atom in  parm.atoms:
      if atom.idx in frag :
         #print(atom.name,atom.idx, frag.index(atom.idx),np.where(frag ==atom.idx ) )
         #atom.charge=newcharges[np.where(frag ==atom.idx )[0][0]]
         atom.charge=newcharges[frag.index(atom.idx)]
      else:
        atom.charge=0

        #print(atom.name, atom.charge )
    return parm

def rename_atoms(n_conf, Calpha, LastAtomInSidechain=[] , extin='.mol2', extout='.mol2' ):

    ABC= list(string.ascii_uppercase)
    for i in range(1,n_conf+1):
        #os.system('antechamber -i %s%s -fi pdb -o %sout%s  -fo pdb' %(i,extin,i,extin))
        print('loading : %s%s'%(i,extin))
        parm = parmed.load_file('%s%s'%(i,extin))
        AmberOrder=['N', 'H','C', 'O', 'CA']
        k=1
        for atom in  parm.atoms:
            if atom.name in AmberOrder :
                atom.name = 'R' +str(k)
                atom.number=55
            k+=1
        prev_order=[-1,-1,-1,-1,-1]
        current=[]
        index=1
        for atom in  parm.atoms:


            if atom.name=='CA' or atom.idx==Calpha :
                atom.name=AmberOrder[4]
                #atom.type='CX'
                prev_order[4]= atom.idx
                atom.number=4
                atom.type='CX'
                for atm in atom.bond_partners:
                    if atm.element==7:
                      prev_order[0]= atm.idx
                      atm.name=AmberOrder[0]
                      atom.number=0
                      for at in atm.bond_partners:
                        if at.element==1 :
                          prev_order[1]= at.idx
                          at.name=AmberOrder[1]
                          at.number=1
                    elif atm.element==6:

                        if 8 in [at.element for at in atm.bond_partners] :
                            prev_order[2]= atm.idx
                            atm.name=AmberOrder[2]
                            atm.number=2
                            for at in atm.bond_partners:
                                if at.element==8 :
                                  prev_order[3]= at.idx
                                  at.name=AmberOrder[3]
                                  at.number=3
                        else:
                            current.append(atm)
                            t=0

                            for at in current :
                               #print(at.name , at.idx , t,  len(current))
                               if len(current) == 1 : t=''
                               else:  t+=1
                               at.name='C'+str(t)+'B'
                               prev_order.append(at.idx)
                               at.number=6
                               #print(at.name , at.idx)


                    elif atm.element==1 :
                         atm.number=4
                         atm.name='HA'


        if -1 in prev_order:
            parm.save( '%s-debug.mol2'%i )
            sys.exit('debug: missing atoms in backbone conformation : %s ' %i)
        if len(current) ==0 : sys.exit('debug: no next atom')

        to_removefromlib=[]
        index+=1
        while len(current) >0:
                Hbound=[]
                #parm.save( '%s-reo.mol2'%index )
                #print('atom : '+str(index))
                for jojo in current :
                    if jojo.type=='CT' : jojo.type ='CU'
                current, index , parm, prev_order,to_removefromlib=find_next(current, index, parm, prev_order, to_removefromlib)
                for jojo in current :
                    if  jojo.idx == atomtocut  :      atom.name='CY'

                #print('after atom : '+str(index))

                #print(to_removefromlib)
                for atom in current :

                    if  atom.idx == atomtocut or len(to_removefromlib) > 0 : to_removefromlib.append(atom.idx)

                    #print(LastAtomInSidechain,atom.idx, len(current))
                    current=np.array(current)  # ###   HERE a bug was found using the pop function! the remaining handle in the array were changed ?!?
                    if  atom.idx in LastAtomInSidechain   :
                        #if len(to_removefromlib) > 0 : to_removefromlib.append(atom.idx+1)
                        last=False
                        #last = current.pop(current.index(atom))
                        #last = current[current.index(atom)]
                        #last=current[np.where(current==atom)]
                        current= np.delete(current,np.where(current==atom))
                        last=atom
                        #print('LAST!!', len(current),last.element, last.name )
                        for at in last.bond_partners:

                            if  at.element==1 :
                                if len(to_removefromlib) > 0 : to_removefromlib.append(at.idx)
                                Hbound.append(at.idx)
                            if  at.element==8 :
                                if len(to_removefromlib) > 0 : to_removefromlib.append(at.idx)
                                at.name= 'OF'


                        if len(Hbound)==1 :j=''
                        else: j=1

                for atm in  parm.atoms:

                            if atm.idx in Hbound :
                                if len(to_removefromlib) > 0 : to_removefromlib.append(atm.idx)
                                atm.name= 'H' +ABC[index-1] +str(j)
                                #print(atm.name)
                                if j!='' :   j+=1
                                atm.number=index+4

                prev_order.append(index)
                #print('after atom end of loop:' +str(index) )
                #print('next :' +str(len(current)))

                Hbound=[]
        new=len(prev_order)+1
        for atom in  parm.atoms:
            if atom.idx not in prev_order :
                atom.number=new
                new+=1

                #for  i in Hbound : prev_order.append(i)
        #parm.atoms.sort(key=lambda x: x.number)
        #for atom in  parm.atoms:
            #print(atom.number)

        #print(prev_order)
        parm.save( '%s-reo.mol2'%i )
    return to_removefromlib

def find_atom_by_name(parm, name):
    for atom in  parm.atoms:
        if atom.name==name:
            return atom
def find_atom_by_idx(parm, idx):
    for atom in  parm.atoms:
        if atom.idx==idx:
            return atom

def find_next(current_atoms, index, parm, previous, to_removefromlib):
    ABC= list(string.ascii_uppercase)
    atom_element= dict({6:'C' , 16:'S' , 8:'O', 7:'N'  })
    Hbound=[]
    Heavybound=[]
    next_atoms=[]
    for current in current_atoms:
        for at in current.bond_partners:
            #print(at.element, at.number)
            if  at.element==1 :
                Hbound.append(at.idx)
                #print(Hbound)
            elif at.element!=1 and at.number==-1 or at.idx not in previous  :
                Heavybound.append(at.idx)
                #print(Heavybound)
                next_atoms.append(at)

        if len(Heavybound) == 1 : s=''
        else : s= 1 #str(len(Heavybound))
        if len(Hbound) == 1 : i=''
        else : i=1
    for atm in parm.atoms:
            if  atm.idx in Heavybound :
                #print(in1.mol2 dex)
                atm.name=atom_element[atm.element] + str(s) +ABC[index]

                #print(atom_element[atm.element] + str(s) +ABC[index])
                if s != '' : s+=1
                #print(atm.idx)
                atm.number=index+5
            elif atm.idx in Hbound :
                if len(to_removefromlib)>0 : to_removefromlib.append(atm.idx)
                atm.name= 'H' +ABC[index-1] +str(i)
                if i!='' : i+=1
                atm.number=index+4
    prev_order= previous+Hbound+Heavybound
    #if index==6 :sys.exit()
    return  next_atoms, index+1 , parm, prev_order, to_removefromlib


def run_rdkit(smiles, numConfs=20) :
    Chem.MolFromSmiles(smiles)
    m2=Chem.AddHs(m)

    cids = AllChem.EmbedMultipleConfs(m2, numConfs=numConfs)
    res = AllChem.MMFFOptimizeMoleculeConfs(m2)
    rmslist = []
    AllChem.AlignMolConformers(m2, RMSlist=rmslist)
    w = Chem.SDWriter('all_conf.sdf')
    for m in cids : w.write(m)




parser = argparse.ArgumentParser(description='Run multi RESP ')
parser.add_argument('--input',  type=str, help='an sdf file containning few conformations of the molecule')
parser.add_argument('--Calpha',  type=int,  help='Calpha of the molecule, necessary ! ')
parser.add_argument('--CY',  type=int,  default=0,  help='last atom in side chain for the librairy  ')
parser.add_argument('--last',  nargs='+', default=[], help='last atom in side chain to take in account for partial charge calculation, default:none ')
parser.add_argument('--name_res', default='MOL' , help='name of your residu , default MOL' )
args = parser.parse_args()




input=args.input
Calpha=args.Calpha
LastAtomInSidechain=args.last
resname=args.name_res
atomtocut=args.CY




#print(resname)
### parmed numerotation starts at 0 like python ###
Calpha += -1
atomtocut += -1
for i in range(len(LastAtomInSidechain)) : LastAtomInSidechain[i] =int(LastAtomInSidechain[i]) -1
#print(LastAtomInSidechain)
##################################################

if input[-3:]=='.smi':
    file=open(input,'r')
    smile=file.readlines().replace(' ','')
    smile=smile.replace('\n','')
    run_rdkit(smile, numConfs=20)

n_conf= int(sdf2pdb(input))-1
print('There is %s conformations in input file'%n_conf)
if  n_conf == -1 : sys.exit('Error, check your sdf file!')
run_allconf(int(n_conf))

rename_atoms(int(n_conf), Calpha,LastAtomInSidechain, extin='.mol2')
to_removefromlib=rename_atoms(int(n_conf), Calpha,LastAtomInSidechain, '.mol2')
multiple(int(n_conf), to_removefromlib,resname)
multiple_full(int(n_conf))
