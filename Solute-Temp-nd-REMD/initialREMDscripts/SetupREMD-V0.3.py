#--------
# This script setup and makes a 100 ns test-run Gromacs simulation
# using AMBER99SB-ILDN forcefield
# and TIP3P water molecule type
#------------------------------
import numpy as np
import sys
import os
import argparse
import shutil
from shutil import copyfile
import fileinput
import re
import multiprocessing
from collections import Counter
import string
from math import *

templdir='/home/marie/Templates/'

def check_arguments() :
 if os.path.isfile(args.pdbfile)==False:
   sys.exit('specify a valid input file')
 if os.path.exists(templdir)==False:
   sys.exit('check your templates folder')


def get_nTFE(conc):
    # nTFE=C*Vfree*avogadro
    # hyp1 : linear small peptide --> volume occupied max 2 % Vfree=Vbox*0.98--> wont work for a normal protein ..
    volfile=open('volume','r')
    Vfree=float(volfile.readline().split()[4])*0.98
    return ceil(conc*Vfree*6.002)


def setupsimulation(conc=0) :
 n=13
 if not os.path.exists('inputs'):
  os.mkdir('inputs')
 shutil.copyfile(args.pdbfile,'inputs/peptide.pdb' )
 if not os.path.exists('simprep'):
  os.mkdir('simprep')
 shutil.copyfile('%s/mini-peptide.mdp' %(templdir), 'simprep/prep.mdp')
 os.chdir('simprep')
 cmd='pdb2gmx  -quiet -f ../inputs/peptide.pdb -ignh -o peptide.gro << EOF \n 6 \n 1 \n EOF  &> /dev/null'
 if os.system(cmd)<>0 : sys.exit( ' error : try to run /n cd %s /n %s %(os.getcwd(), cmd ) ')
 namegro='peptide'
 cmd='gmx editconf  -quiet -f %s.gro -o %s+box.gro -bt octahedron  -d .3 | grep "new box volume" > volume &> /dev/null ' %(namegro,namegro)
 if os.system(cmd)<>0 : sys.exit( ' error : try to run /n cd %s /n %s %(os.getcwd(), cmd ) ')
 namegro+='+box'

 if conc!=0:
    n=15
    #os.system('genbox   -cp %s.gro -ci %s/TFE.acpype/TFE_GMX.gro -nmol %s -o %s+TFE%s.gro >> /dev/null' %(namegro,templdir,conc,namegro,str(conc)))
    nTFE=int(get_nTFE(conc))
    cmd='genbox  -quiet  -cp %s.gro -ci %s/TFE_lib/TFE_GMX.gro -nmol %s -o %s+TFE%s.gro &> /dev/null' %(namegro,templdir,nTFE,namegro,str(conc))
    if os.system(cmd)<>0 : sys.exit( ' error : try to run /n cd %s /n %s %(os.getcwd(), cmd ) ')
    namegro+='+TFE%s' %(str(conc))
    os.system('mv topol.top topol_old.top')
    topolfile=open('topol_old.top','r' )
    lines=topolfile.readlines()
    topolfile.close()
    topolfile=open('topol.top','a' )
    i=0
    print((len(lines)))
    while i<len(lines) and 'amber99sb-ildn.ff/forcefield.itp' not in lines[i]:
          topolfile.write(lines[i])
          i+=1
    topolfile.writelines(lines[i])
    topolfile.write('; Include TFE parameters \n')
    topolfile.write('#include "%sTFE_lib/TFE_GMX_B15.itp" \n' %(templdir) )
    # topolfile.write('#include "%sTFE.acpype/TFE_GMX.itp" \n' %(templdir) )
    topolfile.writelines(lines[i+1:])
    topolfile.flush()
    topolfile.write('TFE	           %s\n' %(str(nTFE)))
    topolfile.flush()
    topolfile.close
    '''
     cnt=Counter()
     for l in open('%s.gro' %(namegro) ): cnt['TFE'] +=1
     print cnt.get('TFE')/9
    '''
 os.system('cp topol.top topolchk2.top')
 cmd='gmx solvate  -quiet -cp %s.gro -cs spc216.gro -p topol.top -o %s+solv.gro &> /dev/null'%(namegro,namegro)
 if os.system(cmd)<>0 : sys.exit( ' error : try to run /n cd %s /n %s %(os.getcwd(), cmd ) ')
 namegro+='+solv'


 #print 'gmx grompp -f prep.mdp -c %s.gro -p topol.top -o topol.tpr'%(namegro)
 cmd='gmx grompp  -quiet -f prep.mdp -c %s.gro -p topol.top -o topol.tpr -maxwarn 5  &> /dev/null'%(namegro)
 if os.system(cmd) )<>0 : sys.exit( ' error : try to run /n cd %s /n %s %(os.getcwd(), cmd ) ')
 namegro+='+ions'
 cmd='echo %s | gmx genion  -quiet -s topol.tpr -o %s.gro -p topol.top -pname NA -nname CL -neutral  &> /dev/null'%(n,namegro)
 if os.system(cmd) )<>0 : sys.exit( ' error : try to run /n cd %s /n %s %(os.getcwd(), cmd ) ')
 shutil.copyfile('%s.gro' %(namegro), 'confin.gro')
 os.chdir(workdir)


def minimization():
  if not os.path.exists('minimization'):
   os.mkdir('minimization')
  shutil.copyfile('%s/mini-peptide.mdp' %(templdir), 'minimization/min.mdp')
  os.chdir('minimization')
  print()
  os.system('gmx grompp -quiet  -f min.mdp -c ../simprep/confin.gro -p ../simprep/topol.top -o topol.tpr -maxwarn 5')
  os.system('gmx mdrun')
  os.chdir(workdir)


def equilibration():
  if not os.path.exists('equilibration'):
   os.mkdir('equilibration')
  shutil.copyfile('%s/equi-peptide.mdp' %(templdir), 'equilibration/equi.mdp')
  os.chdir('equilibration')
  os.system('gmx grompp  -quiet -f equi.mdp -c ../minimization/confout.gro -p ../simprep/topol.top -o topol.tpr -maxwarn 5')
  os.system('gmx mdrun')
  os.chdir(workdir)

def equilibration2(T=300):
 direquil2='equilibration2/equilibration2-'+str(T)
 if not os.path.exists('equilibration2'):
   os.mkdir('equilibration2')
 if not os.path.exists(direquil2):
  os.mkdir(direquil2)
 conffile=open('%s/equi2-peptide.mdp' %(templdir),'r')
 lines=conffile.readlines()
 if len(lines)==0:
     sys.exit()
 conffile.close()
 conffile_new=open('%s/equi2.mdp' %(direquil2),'w')
 for line in lines:
   conffile_new.write(string.replace(line,'$TEMP',str(T)))
 conffile_new.close()
 os.chdir(direquil2)
 os.system('gmx grompp  -quiet -f equi2.mdp -c ../../equilibration/confout.gro -p ../../simprep/topol.top -o topol.tpr -maxwarn 5')
 os.system('gmx mdrun')
 os.chdir(workdir)

def prodtest():
 if not os.path.exists('prodtest'):
  os.mkdir('prodtest')
 shutil.copyfile('%s/prodtest-peptide.mdp' %(templdir), 'prodtest/prodtest.mdp')
 os.chdir('prodtest')
 os.system('gmx grompp  -quiet -f jprodtest.mdp -c ../equilibration2/equilibration2-300/confout.gro -p ../simprep/topol.top -o topol.tpr -maxwarn 5')
 os.system('gmx mdrun  -quiet')
 os.chdir(workdir)

def REMDtempcalc(n_replicas=8,Tmin=300,Tmax=400):
### Setup a linear temperature distribution according to how many cores are available
### Thi has to be changed
    temperatures=np.arange(Tmin,Tmax,(Tmax-Tmin)/(n_replicas-1))
    print(temperatures)
    return temperatures




def prodREMDsetup(T):
 proddir='prodREMD/prod-'+str(T)
 filesetup='prod-'+str(T)
 if not os.path.exists('prodREMD'):
  os.mkdir('prodREMD')
 if not os.path.exists(proddir):
  os.mkdir(proddir)
 filemdp=open('%s/prodREMD.mdp' %(proddir),'w')
 filemdptemplate=open('%sprodTEMP-peptide.mdp' %(templdir),'r')

 lines=filemdptemplate.readlines()
 filemdptemplate.close()
 for line in lines:

  filemdp.write(string.replace(line,'$TEMP',str(T)))
 filemdp.close()
 os.chdir(proddir)
 os.system('gmx grompp  -quiet -f prodREMD.mdp -c %s/equilibration2/equilibration2-%s/confout.gro -p %s/simprep/topol.top -o topol.tpr -maxwarn 5' %(workdir,T,workdir ))
 os.chdir(workdir)

def prodREMDrun(temperatures):
    proddir='prodREMD/prod-'+str(T)
    os.chdir(proddir)
    os.system('mpirun -np %s  mdrun_mpi -v -multidir prod-%s replex 100' %( n_replicas, temperatures))
    os.chdir(workdir)



############################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='# This script setup an REMD Gromacs simulation \
    /n  using a concentration of TFE (default pure water) \
    /n# using AMBER99SB-ILDN forcefield \
    /n # and TIP3P water molecule type')
    parser.add_argument('--TFAconcentration', default=0 , type=float,   help=' concentration of TFA in M (Seems that maximum is a little less than 1 M, more investigation to come...)')
    parser.add_argument('--pdbfile', default='peptide.pdb' , type=str , help='input PDB file ')
    parser.add_argument('--start_from' ,  default='prep' , help='step to restart from (default system preparation values prep, min, equ1,  equ2, prodtest  )' )
    #parser.add_argument('--prodtest' , default=True )
    parser.add_argument('--n_replicas' , default=1 , type=int)
    parser.add_argument('--Tmin' , default=300 , type=float  )
    parser.add_argument('--Tmax' , default=400 , type=float)
    #parser.add_argument('--n_cores', type=int , default=multiprocessing.cpu_count() , help='number of cores to use, default maximum available'  )
    args = parser.parse_args()

    steps={'prep':1, 'min':2, 'equ1':3,  'equ2':4, 'prod':5}

    workdir=os.getcwd()

    if steps.get(args.start_from)<2:
        check_arguments()
        setupsimulation(args.TFAconcentration)

    if steps.get(args.start_from)<3:
        minimization()

    if steps.get(args.start_from)<4:
        equilibration()
    if steps.get(args.start_from)<5 and args.n_replicas==1:
        equilibration2(args.Tmin)



    if steps.get(args.start_from)<6  and args.n_replicas==1:
        prodtest()

    if steps.get(args.start_from)<6 and args.n_replicas>1:
     temperatures=REMDtempcalc(args.n_replicas,args.Tmin,args.Tmax)

     #run the equilibration in parrallele
     p =Process(target=equilibration2 , args=temperatures)
     p.start()
     p.join()

     for T in temperatures:
     # equilibration2(T)
      prodREMDsetup(T)
     prodREMDrun(temperatures)
