#!/usr/bin/python
############
# Name:
# SetupREMD-VERSION.py
#
###########
############
# Author:
# Marie BLUNTZER
#
###########
# Version :
# 0.4 (developpment)
#
###########
# Description:
# This script setup and makes a 100 ns test-run Gromacs simulation
# using AMBER99SB-ILDN forcefield
# parametrized TFE molecule solvent if required
# and TIP4P-Ewald water molecule type
#
#
###########


import numpy as np
import sys
import subprocess
import os
import argparse
import shutil
from shutil import copyfile
import fileinput
import re
import multiprocessing
from multiprocessing import Process ,Queue ,Pool
from collections import Counter
import string
from math import *


templdir='/home/marie/Stapled_peptide_git/'

def check_arguments() :
 if args.pdbfile==True and args.pdbfile[-3:]!='pdb' and os.path.isfile(args.pdbfile)!=True:
   sys.exit('specify a valid input file')
 if os.path.exists(templdir)==False:
   sys.exit('check your templates folder')


def get_nTFE():
    volfile=open('volume','r')
    Vfree=float(volfile.readline().split(':')[1].split()[0])*0.98
    if args.conc==True:
        print args.conc
        # nTFE=C*Vfree*avogadro
        # hyp1 : linear small peptide --> volume occupied max 2 % Vfree=Vbox*0.98--> wont work for a normal protein ..

        return ceil(args.conc*Vfree*6.002)
    if args.frac!=0 :
        return  ceil(args.frac*Vfree*600.2*1.39/100.02)
    else : return False


def setupsimulation() :
 if not os.path.exists('inputs'):
   os.mkdir('inputs')
#shutil.copyfile(args.pdbfile,'inputs/peptide.pdb' )
 os.chdir('inputs')
 log=open('gmxcmd.bsh', 'w')
 aa={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN', \
     'P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR','s':'AUS',  'r' : 'AUR'  }
 sequence = 'ACE'

 indices_stapled_residues = [i for i, x in enumerate(args.seq) if x == "s" or x == "r" ]
 bond=''
 if args.stapled :
    lib= 'loadoff   %sRES_lib/AKS/AKS.lib \n loadoff  %sRES_lib/AKS/AKS.lib  \n'  %(templdir ,templdir)
    while indices_stapled_residues : bond= 'bond m.%s.CE  m.%s.CE' %( indices_stapled_residues.pop()+2 ,indices_stapled_residues.pop() +2)
 else  : lib= 'loadoff   %sRES_lib/AUS/AUS.lib \n loadoff  %sRES_lib/AUS/AUS.lib  \n'  %(templdir ,templdir)



 for a in args.seq[:-1]:
  sequence=' '.join([sequence , aa[a] ])
 Cterminal = 'C' + aa[args.seq[-1].upper()]
 sequence=' '.join([sequence ,Cterminal ])
 #impose an helical start( up to 20 residues )
 helical=''
 if args.helical : helical = 'impose m { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 } { { "N" "CA" "C" "N" -40.0 } { "C" "N" "CA" "C" -60.0 } }'

 tleap=open('tleap.in' , 'w')
 tleap.writelines(['source leaprc.protein.ff14SB \n' , \
                  'loadamberparams frcmod.ionsjc_tip4pew \n' ,\
                  'loadamberparams %sRES_lib/parameters-atoms-bonds-angles.frcmod \n' %(templdir), \
                  'loadamberparams %sRES_lib/parameters-dihedrals-gaff2.frcmod \n' %(templdir), \
                  '%s\n'  %(lib) ,\
                  'set default PBradii mbondi2  \n'  ,\
                  'source leaprc.water.tip4pew \n' ,\
                  'm = sequence { %s } \n'%(sequence) ,\
                  '%s \n' %(helical) ,\
		  '%s \n' %(bond) ,\
                  'relax m  \n' , \
                  'savepdb m peptide.pdb \n',\
                  'charge m \n' ,\
                  '# solvatebox m prot TIP4PBOX 10.0 iso \n'  , \
                  '# addIonsRand prot Na+ 0 Cl- 0 \n' ,\
                  'saveamberparm m peptide.prmtop peptide.inpcrd \n' , \
                  'quit' ] )
 tleap.close()
 subprocess.call('tleap -f tleap.in' , shell=True)!=0
 subprocess.call('acpype -x peptide.inpcrd -p peptide.prmtop '  , shell=True )

 os.chdir('..')

 if not os.path.exists('simprep'):
  os.mkdir('simprep')
 shutil.copyfile('%s/GROMACS_SIM/mini-peptide.mdp' %(templdir), 'simprep/prep.mdp')

 os.chdir('simprep')
 shutil.copyfile('../inputs/peptide_GMX.top', 'topol.top')
 ## 1 --> TIP3P, 2 --> TIP4P, 3 --> TIP4P-Ewald
 #cmd='echo 3 | pdb2gmx  -quiet -f ../inputs/peptide.pdb -ignh -o peptide.gro -ff amber99sb-ildn '
 #log.write(cmd+'\n')
 #if subprocess.check_output(cmd)!=0 : sys.exit( ' ######## \n !!!  ERROR !!! \n ########## \n DEBUG : try to run : \n cd %s \n %s' %(os.getcwd(), cmd ))
 namegro='peptide_GMX'
 cmd='gmx editconf  -quiet -f ../inputs/%s.gro -o %s+box.gro -bt octahedron  -d 1.5 | grep "new box volume" > volume' %(namegro , namegro)
 log.write(cmd+'\n')
 subprocess.check_output(cmd  , shell=True)
 namegro+='+box'


 subprocess.check_output('mv topol.top topol_old.top' , shell=True  )
 topolfile=open('topol_old.top','r' )
 lines=topolfile.readlines()
 topolfile.close()
 topolfile=open('topol.top','a' )
 i=0
 print((len(lines)))

 while i<len(lines) and 'moleculetype' not in lines[i]:
           topolfile.write(lines[i])
           i+=1
 topolfile.write(';   Include forcefield parameters \n')
 topolfile.write('#include "/home/marie/Stapled_peptide_git/TFE_lib/fffield.itp" \n \n' )

 if args.conc!=False or args.frac!=0:
    #subprocess.check_output('genbox   -cp %s.gro -ci %s/TFE.acpype/TFE_GMX.gro -nmol %s -o %s+TFE%s.gro >> /dev\null' %(namegro,templdir,conc,namegro,str(conc)))
    nTFE=int(get_nTFE())

    cmd='gmx insert-molecules  -quiet  -f %s.gro -ci %s/TFE_lib/TFE_GMX.gro -nmol %s -o %s+TFE.gro ' %(namegro,templdir,nTFE,namegro)
    log.write(cmd+'\n')
    subprocess.check_output(cmd + ' &> /dev/null'  , shell=True )
    namegro+='+TFE'
#    topolfile.write('; Include TFE parameters \n')
#    topolfile.write('#include "%sTFE_lib/TFE_GMX_B15.itp" \n \n' %(templdir) )


# while i<len(lines) and 'system' not in lines[i]:
#         topolfile.write(lines[i])
#         i+=1
# topolfile.write('; Include Water Parameters \n')
# topolfile.write('#include "/usr/local/gromacs/share/gromacs/top/amber99sb-ildn.ff/tip4pew.itp" \n \n')
#    topolfile.write('; Include Ions Parameters \n')
#    topolfile.write('#include "/usr/local/gromacs/share/gromacs/top/amber99sb-ildn.ff/ions.itp" \n \n')

    # topolfile.write('#include "%sTFE.acpype/TFE_GMX.itp" \n' %(templdir) )
 topolfile.writelines(lines[i:])
 topolfile.flush()
 if args.conc!=False or args.frac!=0: topolfile.write('TFE	           %s\n' %(str(nTFE)))
 topolfile.flush()
 topolfile.close

 cmd='gmx solvate  -quiet -cp %s.gro -cs tip4p.gro -p topol.top -o %s+solv.gro'%(namegro,namegro)
 log.write(cmd+'\n')

 subprocess.check_output(cmd + ' &> /dev/null' , shell=True )
 namegro+='+solv'


 #print 'gmx grompp -f prep.mdp -c %s.gro -p topol.top -o topol.tpr'%(namegro)
 cmd='gmx grompp  -quiet -f prep.mdp -c %s.gro -p topol.top -o topol.tpr -maxwarn 10 ' %(namegro)
 log.write(cmd + '\n')
 subprocess.check_output(cmd  , shell=True)!=0
 namegro+='+ions'
 cmd='echo sol | gmx genion  -quiet -s topol.tpr  -o %s.gro -p topol.top -pname NA -nname CL -neutral '%(namegro)
 log.write(cmd+'\n')
 subprocess.check_output(cmd + ' &> /dev/null'   , shell=True)!=0
 shutil.copyfile('%s.gro' %(namegro), 'conf-in.gro')
 os.chdir(workdir)


def minimization():
  if not os.path.exists('minimization'):
   os.mkdir('minimization')
  shutil.copyfile('%s/GROMACS_SIM/mini-peptide.mdp' %(templdir), 'minimization/min.mdp')
  os.chdir('minimization')
  cmd = 'gmx grompp  -f min.mdp -c ../simprep/conf-in.gro -p ../simprep/topol.top -o topol.tpr -maxwarn 5 '
  subprocess.check_output(cmd + ' &> /dev/null'  , shell=True )
  subprocess.check_output('gmx mdrun -gpu_id 0 -nt 3', shell=True)
  os.chdir(workdir)



def equilibration(T=300):
     if not os.path.exists('equilibration'):
          os.mkdir('equilibration')
     direquil2='equilibration/equilibration-'+str(T)
     if not os.path.exists(direquil2):
           os.mkdir(direquil2)

     conffile=open('%s//GROMACS_SIM/equi-peptide.mdp' %(templdir),'r')
     lines=conffile.readlines()
     if len(lines)==0:
             sys.exit()
     conffile.close()
     conffile_new=open('%s/equi.mdp' %(direquil2),'w')
     for line in lines:
           conffile_new.write(line.replace('$TEMP',str(T)))
     conffile_new.close()
     os.chdir(direquil2)
     cmd='gmx grompp -maxwarn 1 -quiet -f equi.mdp -c ../../minimization/confout.gro -p ../../simprep/topol.top -o topol.tpr'
     subprocess.check_output(cmd + ' &> /dev/null ', shell=True)
     subprocess.check_output('gmx mdrun -quiet -gpu_id 0 -nt 3  ', shell=True )
     os.chdir(workdir)
     conffile.close()

def equilibration2(T=300):
 direquil2='equilibration2/equilibration2-'+str(T)
 if not os.path.exists('equilibration2'):
   os.mkdir('equilibration2')
 if not os.path.exists(direquil2):
  os.mkdir(direquil2)
 conffile=open('%s//GROMACS_SIM/equi2-peptide.mdp' %(templdir),'r')
 lines=conffile.readlines()
 if len(lines)==0:
     sys.exit()
 conffile.close()
 conffile_new=open('%s/equi2.mdp' %(direquil2),'w')
 for line in lines:
   conffile_new.write(line.replace('$TEMP',str(T)))
 conffile_new.close()
 os.chdir(direquil2)
 subprocess.check_output('gmx grompp  -quiet -f equi2.mdp -c ../../equilibration/equilibration-%s/confout.gro -p ../../simprep/topol.top -o topol.tpr -maxwarn 5' %(str(T)), shell=True)
 subprocess.check_output('gmx mdrun  -gpu_id 1 -nt 3 ' , shell=True)
 os.chdir(workdir)
 conffile.close()


def prodtest():
 if not os.path.exists('prodtest'):
  os.mkdir('prodtest')
 shutil.copyfile('%s/GROMACS_SIM/prodtest-peptide.mdp' %(templdir), 'prodtest/prodtest.mdp')
 os.chdir('prodtest')
 subprocess.check_output('gmx grompp  -quiet -f prodtest.mdp -c ../equilibration2/equilibration2-300/confout.gro -p ../simprep/topol.top -o topol.tpr -maxwarn 5')
 subprocess.check_output('gmx mdrun  -quiet')
 os.chdir(workdir)

def REMDtempcalc():
### Setup a linear temperature distribution according to how many cores are available
### Thi has to be changed
#temperatures=np.arange(Tmin,Tmax,int(ceil(((Tmax-Tmin)/(n_replicas-1)))))

    #New Set of Temperaratures generated by http://folding.bmc.uu.se/remd/tgenerator.php
    # 'aranged' (by changing the number of atoms) to fit the 12 cores required by ARCHER
    '''temperatures=[300.00, 302.14, 304.30, 306.48, 308.67, 310.88, 313.11, 315.35, 317.60, 319.88, 322.17, 324.48, 326.80, 329.14, 331.50, \
    333.88, 336.27, 338.68, 341.11, 343.56, 346.02, 348.50, 351.00, 353.52, 356.06, 358.62, 361.19, 363.77, 366.38, 369.01, 371.67, 374.34,\
     377.03, 379.74, 382.47, 385.22, 387.98, 390.78, 393.59, 396.42, 399.27, 402.14, 405.03, 407.94, 410.88, 413.83, 416.81, 419.81, 422.83 \
     , 425.88, 428.94, 432.03, 435.15, 438.28, 441.43, 444.61, 447.81, 451.04, 454.28, 457.56, 460.85, 464.17, 467.50, 470.87, 474.26, \
     477.68, 481.12, 484.58, 488.07, 491.59, 495.13, 498.70, 502.28, 505.90, 509.54, 513.22, 516.92, 520.65, 524.40, 528.18, 531.99, 535.82\
     , 539.69, 543.57, 547.48, 551.43, 555.41, 559.41, 563.44, 567.51, 571.60, 575.72, 579.87, 584.05, 588.26, 592.50, 596.77, 601.08, \
     605.41, 609.77, 614.17, 618.60, 623.06, 627.55, 632.07, 636.62, 641.21, 645.83, 650.49, 655.21, 659.93, 664.69, 669.48, 674.30, 679.17\
     , 684.07, 689.00, 693.97, 699.97, 705.01]'''
    '''
    temperatures=[290.00, 290.90, 291.80, 292.70, 293.61, 294.51, 295.42, 296.34, 297.25, 298.17, 299.09, 300.01, 300.93, 301.86, 302.78, \
     303.72, 304.65, 305.58, 306.52, 307.46, 308.41, 309.36, 310.30, 311.25, 312.21, 313.16, 314.12, 315.08, 316.04, 317.01, 317.97, 318.94, \
     319.92, 320.89, 321.87, 322.85, 323.83, 324.82, 325.81, 326.80, 327.79, 328.79, 329.79, 330.79, 331.79, 332.80, 333.80, 334.82, 335.83,\
     336.85, 337.87, 338.89, 339.91, 340.94, 341.97, 343.00, 344.03, 345.07, 346.11, 347.16, 348.20, 349.25, 350.30, 351.36, 352.41, 353.47, \
     354.54, 355.60, 356.67, 357.74, 358.81, 359.89, 360.97, 362.05, 363.13, 364.22, 365.31, 366.40, 367.50, 368.60, 369.70, 370.80, 371.91,\
     373.02, 374.13, 375.25, 376.36, 377.49, 378.61, 379.74, 380.87, 382.00, 383.14, 384.28, 385.42, 386.56, 387.71, 388.86, 390.02, 391.17, \
     392.33, 393.49, 394.66, 395.83, 397.01, 398.18, 399.36, 400.54, 401.73, 402.88, 404.07, 405.27, 406.46, 407.66, 408.87, 410.07, 411.28, \
     412.49, 413.71, 414.93, 416.15, 417.37, 418.60, 419.83, 421.06, 422.30, 423.54, 424.78, 426.03, 427.28, 428.53, 429.79, 431.05, 432.31, 433.58, 434.85, 436.12, \
     437.40, 438.68, 439.96, 441.25, 442.54, 443.83, 445.13, 446.43, 447.73, 449.04, 450.35, 451.66, 452.98, 454.30, 455.63, 456.95, 458.28, 459.62, 460.96, 462.30,\
     463.64, 464.99, 466.34, 467.69, 469.05, 470.41, 471.78, 473.15, 474.52, 475.90, 477.28, 478.66, 480.05, 481.43, 482.83, 484.22, 485.63, 487.03, 488.44, 489.85, \
     491.27, 492.69, 494.11, 495.54, 496.97, 498.40, 499.84, 501.28, 502.73, 504.18, 505.64, 507.10, 508.56, 510.02, 511.49, 512.96, 514.44, 515.92, 517.40, 518.89, 520.38,\
     521.87, 523.37, 524.88, 526.38, 527.90, 529.41, 530.93, 532.45, 533.98, 535.51, 537.04, 538.58, 540.12, 541.66, 543.21, 544.77, 546.33, 547.89, 549.46, 551.03, 552.61, \
     554.19, 555.77, 557.36, 558.95, 560.54, 562.14, 563.75, 565.36, 566.97, 568.58, 570.21, 571.77, 573.40, 575.03, 576.67, 578.31, 579.96, 581.61, 583.26, 584.92, 586.58, \
     588.25, 589.92, 591.59, 593.27, 594.96, 596.65, 598.34, 600.04, 601.74, 603.44, 605.16, 606.87]

    '''
    temperatures=[290.00, 291.88, 293.78, 295.68, 297.59, 299.51, 301.45, 303.39, 305.35, 307.31, 309.29, 311.28, 313.28, 315.29, 317.31, 319.35, 321.40, 323.45, 325.51, 327.59, 329.68, 331.77, 333.88, 336.01, 338.14, 340.29, 342.43, 344.59, 346.77, 348.96, 351.17, 353.39,\
     355.61, 357.85, 360.11, 362.37, 364.65, 366.93, 369.24, 371.55, 373.88, 376.21, 378.56, 380.93, 383.30, 385.69, 388.10, 390.51, 392.94, 395.38, 397.84, 400.32, 402.81, 405.30, 407.81, 410.33, 412.86, 415.41, 417.98, 420.56, 423.16, 425.76, 428.38, 431.02, 433.67, \
     436.33, 439.01, 441.70, 444.41, 447.13, 449.86, 452.62, 455.38, 458.17, 460.97, 463.81, 466.63, 469.46, 472.32, 475.20, 478.09, 480.99, 483.91, 486.81, 489.77, 492.74, 495.73, 498.73, 501.75, 504.78, 507.84, 510.90, 513.99, 517.09, 520.20, 523.34, 526.49, 529.66, 532.85,\
     536.05, 539.27, 542.51, 545.77, 549.04, 552.33, 555.64, 558.98, 562.32, 565.69, 569.08, 572.48, 575.90, 579.34, 582.80]

    return temperatures[:108]


def prodREMDsetup(i,T):

     proddir='prodREMD/'
     filesetup='prod-'+str(i)
     if not os.path.exists('prodREMD'):
      os.mkdir('prodREMD')

     filemdptemplate=open('%s/GROMACS_SIM/prodTEMP-peptide.mdp' %(templdir),'r')
     filemdp=open('prodREMD/prod_%s.mdp' %('{:03d}'.format(i)), 'w')
     lines=filemdptemplate.readlines()
     filemdptemplate.close()
     for line in lines:

      filemdp.write(line.replace('$TEMP',str(T)))
     filemdp.close()

     shutil.copyfile('./equilibration2/equilibration2-%s/confout.gro' %(T) ,  './prodREMD/conf_%s.gro' %('{:03d}'.format(i)))




def prodREMDrun(temperatures):
    proddir='prodREMD/prod-'+str(T)
    os.chdir(proddir)
    subprocess.check_output('mpirun -np %s  mdrun_mpi -v -multidir prod-%s replex 100' %( args.n_replicas, temperatures))
    os.chdir(workdir)



############################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='# This script setup an REMD Gromacs simulation \
    \n  using a concentration of TFE (default pure water) \
    \n# using AMBER99SB-ILDN forcefield \
    \n # and TIP4P-Ewald water molecule type')
    parser.add_argument('--conc', type=float,  default=False,  help=' concentration of TFA in M (Seems that maximum is a little less than 1 M, more investigation to come...)')
    parser.add_argument('--frac', default=0.25 , type=float,   help=' concentration of TFA in M (Seems that maximum is a little less than 1 M, more investigation to come...)')
    parser.add_argument('--pdbfile', default='peptide.pdb' , type=str , help='input PDB file ')
    parser.add_argument('--seq'  , required=True, type=str , help='input sequence')
    parser.add_argument('--start_from' ,  default='prep' , help='step to restart from (default system preparation values prep, min, equ1,  equ2, prodtest  )' )
    #parser.add_argument('--prodtest' , default=True )
    parser.add_argument('--helical' , action='store_true' )
    parser.add_argument('--noREMD' , action='store_true'  )
    parser.add_argument('--stapled' , action='store_true'  )
    #parser.add_argument('--n_cores', type=int , default=multiprocessing.cpu_count() , help='number of cores to use, default maximum available'  )
    args = parser.parse_args()

    steps={'prep':1, 'min':2, 'equ1':3,  'equ2':4, 'prod':5}

    workdir=os.getcwd()
    if args.noREMD==True : REMD = False
    else : REMD=True
    if steps.get(args.start_from)<2:
        check_arguments()
        setupsimulation()
        print('run minimization')
    if steps.get(args.start_from)<3:
        minimization()

    if steps.get(args.start_from)<4:


        temperatures=REMDtempcalc()
        print(temperatures)

        for T in temperatures:
              print( 'equilibration2 at', T)
              equilibration(T)

    if steps.get(args.start_from)<5 and REMD==0:
        equilibration2()



    if steps.get(args.start_from)<6 and REMD==1:

     temperatures=REMDtempcalc()
     print(temperatures)

     for T in temperatures:
      print( 'equilibration2 at', T)

      equilibration2(T)
     '''
     ### paralezitian 'works', but also makes gromacs crashe even with 2 processes...'
     pool = Pool(processes=2)
     print (pool.map(equilibration2,temperatures) )


     #### Other way same issue
     processes = [multiprocessing.Process(target=equilibration2, args=temperatures)]
     # Run processes
     for p in processes:
         p.start()
     # Exit the completed processes
     for p in processes:
         p.join()
     '''
     for i in range(len(temperatures)):
        prodREMDsetup(i,temperatures[i])
     shutil.copyfile( 'simprep/topol.top' , 'prodREMD/topol.top' )

    # prodREMDrun(temperatures)
