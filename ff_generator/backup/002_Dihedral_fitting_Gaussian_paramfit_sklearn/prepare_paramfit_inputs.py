d=' \
############################ \
Prepare Paramfit inputs input from gaussians input files and the mol2 file used to produce those \
!!!! It HAS TO BE THE SAME MOL2 FILE OR AT LEAST ATOMS HAVE TO BE IN THE SAME ORDER! !!!!! \
Call process_gaussians_files which produces  \
        the stationary structures from gaussians,\
        -energy files \
        - listof files (succesful gaussians calculations) \
this script produce: \
    - the initial prmtop with gaff2 parameters, \
            KNOWN BUG :  you will have to remove  atoms types in the initial mol2 file...  \
    - mdcrd files \
    - corresponding formated paramfit input files \
############################\
Author : Marie Bluntzer \
############################\
'

from  process_gaussians_files  import extract_stationary_structures
from  get_all_dihedral import find_dihedrals
import mdtraj as md
import shutil
import sys
import os.path

def get_gaff2_parmcheck(input):

     os.system('parmchk2 -f mol2 -i %s -o gaff2.frcmod  -s  gaff2' %(input))

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



def make_traj(prefix):
    list_pdb=open(prefix +'-listfiles')
    lines=list_pdb.readlines()
    if len(lines) > 0 : first=lines.pop(0).replace(' \n','')
    else :
        #sys.exit('BOUH!' + prefix +   str(lines) )
        return False
    t = md.load(first, top='input.prmtop')

    for i in range(len(lines)):
        filename=lines.pop(0).replace(' \n','')
        print(filename)
        tnext= md.load(filename.replace(' \n',''), top='input.prmtop')
        t= md.join([t , tnext], check_topology=True, discard_overlapping_frames=False)
    t.save_mdcrd(prefix+ '-traj.mdcrd' )
    t.save_netcdf(prefix+ '-traj.nc' )
    return t.n_frames

def prepare_paramfit_job_files (prefix, nframes) :
    ## This fonction prepares
    # - fit_K.in for paramfit : find constant to fit Kcal to Hartree
    # fit_P.in for torsion

    temp_fitK= open(os.path.dirname(os.path.realpath(sys.argv[0]))+'/fit_K_A.in' , 'r')
    temp_fitP= open(os.path.dirname(os.path.realpath(sys.argv[0]))+'/job_fit_1.in' , 'r')
    new_fitK= open(prefix + '-K-fit.in','w')
    new_fitP= open(prefix + '-P-fit.in','w')

    for line in temp_fitK.readlines():
        new_fitK.writelines(line.replace('$PREFIX', prefix).replace('$NSTRUCT',  str(nframes)))
    temp_fitK.close()
    for line in temp_fitP.readlines():
        new_fitP.writelines(line.replace('$PREFIX', prefix).replace('$NSTRUCT',  str(nframes)).replace('$PARAMFILE', prefix + '-P-fit.in' ))
    temp_fitP.close()


def prepare_paramfit_param_files (prefix,torsion_names, all_dihedrals_type) :
    template= os.path.dirname(os.path.realpath(sys.argv[0])) +'/param.in'
    jobfile='parameter_'+prefix +'.in'
    shutil.copy(template ,jobfile )
    temp_fit=open(jobfile, 'a')
    index=torsion_names.index(prefix[0:-1])
    print(prefix , all_dihedrals_type)
    added = []
    count=[]
    for types in all_dihedrals_type[index]:
        #print(types)
        #    for mult in range(types[4]):
        if types[0] + types[1] +types[2]+types[3] not in added :
        # order != mult   if order 1 3 6 / mult = 3
            count.append(0)
            added.append(types[0] + types[1] +types[2]+types[3])
            mult=0
        else :
            count[added.index(types[0] + types[1] +types[2]+types[3] )]+=1
            mult=count[added.index(types[0] + types[1] +types[2]+types[3] )]
        temp_fit.writelines('	%s %s %s %s	%s	1	0	0\n' %(types[0], types[1], types[2],types[3] ,mult) )
    temp_fit.close()

def prepare_mdgx_job_files (prefix,torsion_names, all_dihedrals_type) :
    ## This fonction prepares
    # - fit_K.in for paramfit : find constant to fit Kcal to Hartree
    # fit_P.in for torsion

    new_fit= open(prefix + '-mdgx-fit.in','w')

    new_fit.writelines('&files\n')
    new_fit.writelines(' -parm   /home/marie/Utilities/amber16/dat/leap/parm/frcmod.ff14SB\n')
    new_fit.writelines(' -o      fit.out\n')
    new_fit.writelines(' -fmod amber.frcmod\n')
    new_fit.writelines('&end\n\n')

    new_fit.writelines('&param\n')
    new_fit.writelines(' System  input.prmtop %s-traj.nc %s-energySP.dat\n' %(prefix, prefix)   )
    new_fit.writelines(' ParmOutput frcmod\n')
    new_fit.writelines(' eunits    hartree,\n')
    new_fit.writelines(' accrep    report.m\n')
    new_fit.writelines(' verbose    1,\n\n')
    new_fit.writelines(' % Angle fitting input\n')
    angle_added=[]
    dih_added=[]
    index=torsion_names.index(prefix[0:-1])
    for types in all_dihedrals_type[index]:

        if sorted([types[0] , types[1] ,types[2] ])  not in angle_added :
         angle_added.append(sorted([types[0], types[1] ,types[2]]))
         new_fit.writelines('  fita       %s %s %s\n'%(types[0], types[1] ,types[2]))
        if sorted([types[1] ,types[2], types[3] ])  not in angle_added :
         angle_added.append(sorted([types[1] ,types[2],types[3]]))
         new_fit.writelines('  fita       %s %s %s\n'%(types[0], types[1] ,types[2]))

    new_fit.writelines('  FitAnglEq  1,\n  arst       0.0002,\n  arstcpl    1.0,\n\n')
    new_fit.writelines(' % Torsion fitting input\n')

    for types in all_dihedrals_type[index]:

        if [types[0] , types[1] ,types[2],types[4] ]  not in dih_added :
         dih_added.append([types[0], types[1] ,types[2],types[3] ])
         new_fit.writelines('  fith       %s %s %s %s\n'%(types[0], types[1] ,types[2],types[3]))
    new_fit.writelines('  hrst       0.0002,\n&end')
    new_fit.close()
    print(angle_added)




if __name__=='__main__':
    if len(sys.argv) != 2:
        sys.exit('need the name of inputfile')
    inputmol2=sys.argv[1]

    all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name ,torsion_names , dihedrals_heavy_index= find_dihedrals(inputmol2)
    add_order()

    for p in torsion_names:
        prefix = p + '_'
        extract_stationary_structures(prefix, Opt=False)
        nframes= make_traj(prefix)
    if  nframes :
        prepare_paramfit_job_files (prefix, nframes)
        prepare_paramfit_param_files (prefix,torsion_names, all_dihedrals_type)
        prepare_mdgx_job_files (prefix,torsion_names, all_dihedrals_type)
