import mdtraj as md
import sklearn.metrics as sk1
import parmed
#import 02_prepare_paramfit_inputs
from  get_all_dihedral import find_dihedrals
from  energyDecomposition import singlepoint
from prepare_paramfit_inputs import *

import sys
import os
import subprocess
import numpy as np
import shutil


import pylab as plt

import seaborn as sns



def add_new_parms(inputfile,newfrcmod):
    tleapinput=open('temp.in','w')
    tleapinput.writelines(['loadAmberParams best.frcmod \n',  'loadAmberParams %s \n' %newfrcmod, \
    'source leaprc.protein.ff14SB \n'  ,  'new = loadmol2 %s \n '  %(inputfile) ,\
    'saveAmberParm  new  newtopol.prmtop new.inpcrd\n ', 'quit'])
    tleapinput.close()
    os.system( 'tleap -f temp.in >tleap' )

def get_converted_QM_data(prefix):
    if not os.path.isfile(prefix+'-energySP.dat')  :sys.exit(prefix+'-energySP.dat   energy file not found'  )
    qm = open(prefix+'-energySP.dat','r')
    #1 Hartree	627.509 kcal mol-1
    QM_data=[ float(i)*627.509 for i in qm.readlines() ]
    minimum = min(QM_data)
    converted_QM= [val - minimum for val in QM_data]
    return converted_QM

def find_outiers (data,threshold = 50 ):

    datarollingmedian=[]
    for i in range(len(data))  :
        # Periodic data : easy to get rolling median :-)
        datarollingmedian.append(np.median(data[i-5:i+5]))

    data=np.array(data)
    datarollingmedian= np.array (datarollingmedian)
    difference= np.subtract(data, datarollingmedian)

    outlier_idx = difference > threshold

    return outlier_idx

def add_order():
    F= open('amber.frcmod', 'r')
    N=open('oneorder.frcmod','w')
    line=F.readline()
    while line[0:4] !='DIHE':
        N.writelines(line)
        line=F.readline()

    N.writelines(line)
    previous=F.readline()
    line=F.readline()
    next=F.readline()

    while line[:6] not in ['IMPROP','NONBON' ]:
        if line[0:11]!=previous[0:11]:
                N.writelines( previous[0:11]+ '   1    0.000       180.000          -6.000     \n')
        N.writelines(previous)
        previous=line
        line=next
        next=F.readline()


    N.writelines(previous)
    N.writelines(line)
    N.writelines(F.readlines())
    N.close()

def single_point_calculations(topol, prefix):
    #print(topol, prefix+'-traj.mdcrd')
    singlepointenergy=singlepoint(topol,[ prefix+'-traj.mdcrd'])
    return singlepointenergy

	#return K
def paramfit_fit_K(prefix,topol):
    cmd='paramfit -i %s-K-fit.in -p  %s   -c %s-traj.mdcrd  -q %s-energySP.dat   > paramrun.out' %(prefix, topol ,prefix,prefix)
    os.system(cmd)
    file=open('paramrun.out','r')
    for line in file.readlines():
        if '*K' in line :
            K=float(line.split()[2] )
    file.close()
    return K

def run_paramfit(topolin, prefix,K,mode,optimizations, max_gen, gen_to_conv , gen_simplex, gen_without_simplex, mutation_rate , parent_percent, search_space , dihedralstorsions ):
    # K value of K
    # mode = BOTH , GENETIC, SIMPLEX
    # mutation_rate

    old_fitP= open(prefix + '-P-fit.in','r')
    new_fitP= open(prefix + '-P-fit-run.in','w')
    for line in old_fitP.readlines():
        newline=line.replace('$K', str(K)).replace('$ALGORITHM', str(mode))
        newline=newline.replace('$SSPACE', str(search_space)).replace('$PARENT', str(parent_percent) ).replace('$MUT_RATE',str(mutation_rate))
        newline=newline.replace('$GEN_WT_SIMPLEX', str(gen_without_simplex)).replace('$GEN_SIMPLEX',str(gen_simplex ))
        newline=newline.replace('$MAX_GEN',str(max_gen)).replace('$GEN_CONV',str(gen_to_conv)).replace('$OPT', str(optimizations))
        new_fitP.writelines(newline )
    new_fitP.close()
    cmd='paramfit -i %s-P-fit-run.in -p %s -c %s-traj.mdcrd -q %s-energySP.dat >jobparamfit ' %(prefix, topol ,prefix,prefix)
    #cmd='paramfit -i %s-P-fit-run.in -p %s -c %s-traj.mdcrd -q %s-energySP.dat' %(prefix, topol, prefix ,prefix)
    #os.system(cmd)
    subprocess.check_call(cmd, shell=True)

    file=open(prefix+'paramfit.frcmod','r')
    top = parmed.load_file(topol)
    previousline=[0,0,0,0]
    new=[]
    for line in file.readlines():

            s=line.split('-')

            if len( s) ==4:
                a=[s[0].replace(' ',''),s[1].replace(' ',''),s[2].replace(' ',''),s[3].split()[0] ]
                phi=float(s[3].split()[1])
                phase=float(s[3].split()[2])
                if phi < 0 :
                    phi=-phi
                    phase= phase + 180 %180    #phase == 0 or 180 ::::   (phi*cos(angle +180)= -phi*cos(angle)
                    '''
                if a  != previousline :
                    parmed.tools.deleteDihedral(top, "@%" +a[0] ,"@%"+a[1], "@%" +a[2] ,"@%"+a[3]).execute()  #delete all multiplicity for one dihedral

                parmed.tools.addDihedral(top, "@%"+a[0],"@%"+a[1],"@%"+a[2],"@%"+a[3], phi , s[3].split()[3] , phase ).execute()
                previousline=[a[0],a[1],a[2],a[3]]
                    '''
                new.append([ [ a[0],a[1],a[2],a[3]], [phi , s[3].split()[3] , phase]])
    '''
    for dihe in dihedralstorsions:
        print(dihe)
        ai = [ str(dihe.atom1.idx), str(dihe.atom2.idx), str(dihe.atom3.idx), str(dihe.atom4.idx)]
        at = [ dihe.atom1.type, dihe.atom2.type,dihe.atom3.type,dihe.atom4.type]
        parmed.tools.deleteDihedral(top, "@" +ai[0] ,"@"+ai[1], "@" +ai[2] ,"@"+ai[3]).execute()

        phi=[a[1][0] for a  in new if a[0] == at]
        per=[a[1][1] for a  in new if a[0] == at][0]
        phase=[a[1][2] for a  in new if a[0] == at][0]

        parmed.tools.addDihedral(top, "@" +ai[0] ,"@"+ai[1], "@" +ai[2] ,"@"+ai[3], phi, per ,phase  )
    if os.path.isfile('newtopol.prmtop')==True:os.remove('newtopol.prmtop')
    top.save('newtopol.prmtop')
    '''
    add_new_parms(inputfile, prefix+'paramfit.frcmod' )
    new=parmed.load_file('newtopol.prmtop')
    parmed.tools.writeFrcmod(new,'new.frcmod').execute()

if __name__=='__main__' :


    if len(sys.argv) != 2:
        sys.exit('need the name of inputfile')
    inputfile=sys.argv[1]


    all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name ,torsion_names , dihedrals_heavy_index= find_dihedrals(inputfile)
    add_order()
    to_be_removed=[]
    for p in torsion_names:
        prefix = p + '_'
        extract_stationary_structures(prefix, Opt=False)

        nframes= make_traj(prefix)
        if nframes != False :
            prepare_paramfit_job_files (prefix, nframes)
            prepare_paramfit_param_files (prefix,torsion_names, all_dihedrals_type)

            prepare_mdgx_job_files (prefix,torsion_names, all_dihedrals_type)
        else :
            to_be_removed.append(torsion_names.index(p)-len(to_be_removed))
            print(p)
    print(torsion_names)
    for i in to_be_removed:

        torsion_names.pop(i)
        print(torsion_names)
###########################################################





    all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names, dihedrals_heavy_index=find_dihedrals(inputfile)
    topol='input.prmtop'
    #print (sys.argv[1])
    #print(topol)
    besttopol='input.prmtop'
    add_order()
    shutil.copy('oneorder.frcmod' ,'best.frcmod')

    for prefix in  torsion_names:
        #print(prefix)
    #for prefix in  ['PSI', 'PHI','CHI1' , 'CHI2', 'CHI3', 'CHI4']:
        K=paramfit_fit_K(prefix+'_', topol)
        #print(K
        plt.figure()
        prefix += '_'
        mode='BOTH'
        '''
        optimizations=5     #30
        max_gen=5         #50
        gen_to_conv=3     #50
        gen_simplex=5    #20
        gen_without_simplex=0
        mutation_rate=0.1
        parent_percent=0.80
        search_space=0.01
        conv_limit= 1.0E-15

        '''
        ## DEBUG:
        conv_limit= 1.0E-15
        optimizations=10     #30
        max_gen=50         #50
        gen_to_conv=10      #50
        gen_simplex=5    #20
        gen_without_simplex=5
        mutation_rate=0.1
        parent_percent=0.80
        search_space=0.10
        #'''

#        solutions=[1,2,3,4,5]
        solutions = [1,2,1,2]
        sol=0
        singlepointenergy=single_point_calculations(topol, prefix)
        QM_data=get_converted_QM_data(prefix)
        plt.errorbar( range(len(singlepointenergy)),singlepointenergy,label='gaff2',linewidth=1)
        #print(QM_data)
        plt.errorbar( range(len(QM_data)),QM_data,label='QM', linewidth=1)

        count=0
        while len(solutions) != 0 :
            all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name ,torsion_names , dihedrals_heavy_index= find_dihedrals(inputfile,'best.frcmod')
            prepare_paramfit_param_files(prefix , torsion_names, all_dihedrals_type)

            run_paramfit(besttopol, prefix,K,mode,optimizations, max_gen, gen_to_conv , gen_simplex,\
            gen_without_simplex, mutation_rate , parent_percent, search_space, all_dihedrals[torsion_names.index(prefix[0:-1])])


            singlepointenergy=single_point_calculations(besttopol, prefix)
            singlepointenergyparamfit=single_point_calculations('newtopol.prmtop', prefix)
            score_before= sk1.r2_score(QM_data, singlepointenergy)
            score_after =  sk1.r2_score(QM_data, singlepointenergyparamfit )
            #print(score_before , score_after)
            if score_after >  score_before :
                shutil.copy('newtopol.prmtop' ,'best.prmtop')
                besttopol='best.prmtop'
                shutil.copy('new.frcmod' ,'best.frcmod')
                plt.errorbar( range(len(singlepointenergyparamfit)),singlepointenergyparamfit,label=prefix,linewidth=1)
                #legend.append(first_line[col])
            if score_after > 0.985 or score_before > 0.999  : solutions=[]
            if score_after < 0 and count<3 :  # Add the oreder in otpimisation
                solutions.append(3)
                if count >2   :
                    all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name ,torsion_names , dihedrals_heavy_index= find_dihedrals(inputmol2,'best.frcmod')
                    prepare_paramfit_param_files(prefix , torsion_names, all_dihedrals_type)
            #else :
                #sol = solutions.pop(0)  #got to the next step
                #print( prefix )
            if sol==0:

                pass
                #Default values
            if sol==1 :
                # INCREASE convergence limit on the simplex ALGORITHM
                conv_limit *= 5E-1
            if sol==2 and parent_percent > 0.5 :
                #Try Genetic algorithm
                conv_limit = 1.0E-15
                parent_percent+=-0.1
                gen_simplex=20
            if sol==3  and parent_percent > 0.5 :
                conv_limit = 1.0E-8
                parent_percent+=-0.1
                count+=1
                mutation_rate=0.4
                optimizations=10     #30
                max_gen=20         #50
                gen_to_conv=10     #50
                gen_simplex=20  #20
                gen_without_simplex=10
            else : sol = solutions.pop()
        plt.legend()
        plt.savefig('%s.png' %(prefix), bbox_inches='tight')
        plt.close()

'''
            if sol==3:
                gen_without_simplex+=10
                gen_simplex=30
                mutation_rate+=0.1
           if sol==4:
                gen_without_simplex+=10
                search_space+=0.1
            if sol==5:
                optimizations+=20
'''
