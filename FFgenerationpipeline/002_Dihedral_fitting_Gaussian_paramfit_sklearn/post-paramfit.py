#!/usr/bin/env python
# coding: utf-8

# In[19]:


import numpy as np
from lmfit import Minimizer, Parameters, report_fit
import math
import mdtraj


import numpy
#import MDAnalysis as MDA
import parmed
import os
import sys
import itertools

# In[74]:

def make_missing_parms_parmck2(inputfile, frcmod='amber.frcmod'):
    os.system( 'parmchk2 -i  %s  -o frcmod -f mol2 -s 2 -o amber.frcmod -a "Y" ' %(inputfile))
    tleapinput=open('temp.in','w')
    tleapinput.writelines([  'source leaprc.protein.ff14SB \n'  , 'loadAmberParams %s \n'%(frcmod), 'input = loadmol2 %s \n '  %(inputfile ) ,    'saveAmberParm  input  input.prmtop input.inpcrd\n ', 'quit'])
    tleapinput.close()
    os.system( 'tleap -f temp.in >tleap' )
inputfile='newcharges.mol2'
make_missing_parms_parmck2(inputfile)



def find_dihedrals(inputfile, frcmod='amber.frcmod'):
    make_missing_parms_parmck2(inputfile, frcmod=frcmod)
    all_dihedrals=[[]]
    dihedrals_heavy=[[]]
    dihedrals_heavy_index=[]
    dihe=[]
    heavy_bonds=[]
    dihedrals_heavy_centralbond_name=[]
    dihedrals_heavy_name=[]
    topol = parmed.load_file( 'input.prmtop')
    for resid in topol.residues:
             for atom in resid.atoms:

                 if atom.atomic_number in [6, 7, 8, 16]: # C, N , O and S

                    for bond in atom.bonds :
                        b = [bond.atom1.idx , bond.atom2.idx]
                        b.sort()
                        if b not in heavy_bonds and bond.atom1.atomic_number in [6, 7, 8, 16] and   bond.atom2.atomic_number in[6, 7, 8, 16]                         and sorted([bond.atom1.type ,  bond.atom2.type]) not in [sorted([e[0], e[1]]) for e in  itertools.combinations_with_replacement(('CA',  'C' ,'CV' , 'N*' , 'NB','O' ), 2) ]                        and bond.atom1.name[-1].isalpha() and bond.atom2.name[-1].isalpha() :
                            heavy_bonds.append(b)




    dihe= topol.dihedrals

    for bond in heavy_bonds:
        if all_dihedrals[-1] != [] :
            all_dihedrals.append([])
            dihedrals_heavy.append([])
        for i  in range (len(dihe)) :
             d = [dihe[i].atom2.idx,dihe[i].atom3.idx]
             d.sort()
             if bond ==d :
                 #all_dihedrals[-1].append(dihe[i])
                 if dihe[i] not in  all_dihedrals[-1]:

                     ans = [' "%s" "%s" "%s" "%s" ' %( dihe[i].atom1.name, dihe[i].atom2.name,dihe[i].atom3.name,dihe[i].atom4.name)]
                     an = [ dihe[i].atom1.name, dihe[i].atom2.name,dihe[i].atom3.name,dihe[i].atom4.name]

                     ai = [ dihe[i].atom1.idx, dihe[i].atom2.idx,dihe[i].atom3.idx,dihe[i].atom4.idx]
                     #ai = [ ' "%s" "%s" "%s" "%s" ' %(dihe[i].atom1.idx, dihe[i].atom2.idx,dihe[i].atom3.idx,dihe[i].atom4.idx)]
                     all_dihedrals[-1].append(dihe[i])
                     if dihe[i].atom1.atomic_number in [6, 7, 8, 16] and  dihe[i].atom4.atomic_number in [6, 7, 8, 16] and ai not in  dihedrals_heavy_index:
                            dihedrals_heavy[-1].append([dihe[i]])
#                            if dihe[i].atom2.type != 'CA' and dihe[i].atom3.type != 'CA'  :  # no aromatic or it will be a mess!
                            d=[an[1], an[2]]
                            d.sort()

                            if d not in dihedrals_heavy_centralbond_name   :

                             dihedrals_heavy_centralbond_name.append(d)
                             dihedrals_heavy_index.append(ai)
                             dihedrals_heavy_name.append( ans)
    torsion_names=[]
    j=1
    all_dihedrals_type=[]
    all_dihedrals_type_index=[]
#    print(dihedrals_heavy_centralbond_name)
    for bond in dihedrals_heavy_centralbond_name:
        bond.sort()
        all_dihedrals_type.append([])
        all_dihedrals_type_index.append([])
        for i  in range (len(dihe)) :
            if  dihe[i].improper == False :
                d = [dihe[i].atom2.name,dihe[i].atom3.name]
                d.sort()
                if bond ==d :

                    #at = [dihe[i].atom1.type,dihe[i].atom2.type,dihe[i].atom3.type,dihe[i].atom4.type , dihe[i].type.per]

                    at = [dihe[i].atom1.type,dihe[i].atom2.type,dihe[i].atom3.type,dihe[i].atom4.type ]
                    if at not in all_dihedrals_type[-1]  : #and ai not in  dihedrals_heavy_index:
                        all_dihedrals_type[-1].append(at)
                    at_index = [dihe[i].atom1.idx,dihe[i].atom2.idx,dihe[i].atom3.idx,dihe[i].atom4.idx , all_dihedrals_type[-1].index(at) ]
                    #print(all_dihedrals_type_index[-1], )
                    if at_index not in all_dihedrals_type_index[-1]  :
                        all_dihedrals_type_index[-1].append(at_index)
                    ### ---> return all indexes of atoms around the same bond and last column is the type of dihedral by reference to all_dihedrals_type

        if ['CA' ,'N'] == bond:
            torsion_names.append('PHI')
        elif [ "C" , "CA" ] == bond :
            torsion_names.append('PSI')
        else :
            torsion_names.append('CHI'+str(j) )
            j+=1

    return all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names, dihedrals_heavy_index ,all_dihedrals_type_index











# In[429]:


def plot(data, final,x):
    try:
        import matplotlib.pyplot as plt
        plt.plot(x.T[0], data, 'k+')
        plt.plot(x.T[0], final, 'r')
        #plt.plot( data2, 'b')
        plt.show()
    except ImportError:
        pass


# In[466]:


def write_frcmod(result, torsion_name , all_dihedrals_type,all_dihedrals_type_index):
    frcmod_file=open('%s.frcmod' %(torsion_name),'w')
    frcmod_file.writelines('\nDIHE\n')
    new_frcmod=[]
    threshold=0.1

    for parm  in result.params:

        if  parm[0:3]== 'amp' and result.params[parm].value != 0 and parm[-3:] != 'bis' : # dependency.keys():
            #print( parm[3] , parm[-1]  ,  result.params['order'+parm[3:]].value, result.params[parm].value)
            if math.fabs(result.params[parm].value) > threshold :
                #new_frcmod.append([all_dihedrals_type[torsion][all_dihedrals_type_index[torsion][int(parm[-1])][-1]],  result.params[parm].value ,  result.params['order'+parm[3:]].value])
                print( [i.ljust(2) for i in parm.split('_')[-4:]])
                new_frcmod.append(['-'.join([i.ljust(2) for i in parm.split('_')[-4:]]),  result.params[parm].value ,  result.params['order'+parm[3:]].value])
            #else :
                #print(result.params[parm].value )

    new_frcmod.sort(key=lambda x: x[0])

    for l in  range(len(new_frcmod)) :
        line=new_frcmod[l]

        if line[1] <0 :  phase = '180'
        else : phase = '  0'
        if l != len(new_frcmod)-1 : nextdih=new_frcmod[l+1][-1]
        else : nextdih=False
        if line[-1] == nextdih  : s= -1
        else : s= 1
        frcmod_file.writelines('%s  1  %s   %s   %s \n' %(line[0] , str(math.fabs(line[1]))[0:6], phase , s*line[2]  ))

        #print( '%s   1  %s   %s   %s' %(line[0][0].ljust(2), line[0][1].ljust(2) ,line[0][2].ljust(2) ,line[0][3].ljust(2) , str(math.fabs(line[1]))[0:6], phase , s*line[2]  ))
#write_frcmod(result,'test', all_dihedrals_type=all_dihedrals_type,all_dihedrals_type_index=all_dihedrals_type_index)


# In[490]:


import parmed
def add_new_parms(inputfile,newfrcmod):
    tleapinput=open('temp.in','w')
    tleapinput.writelines(['source leaprc.protein.ff14SB \n' , 'source leaprc.gaff2 \n' ,'loadAmberParams gaff2.frcmod',  ])
    for name in newfrcmod:
         tleapinput.writelines( 'loadAmberParams %s.frcmod \n' % name  )
    tleapinput.writelines( [ 'new = loadmol2 %s \n '  %(inputfile) ,    'saveAmberParm  new  newtopol.prmtop new.inpcrd\n ', 'quit'])
    tleapinput.close()
    os.system( 'tleap -f temp.in >tleap' )
    new=parmed.load_file('newtopol.prmtop')
    parmed.tools.writeFrcmod(new,'new.frcmod').execute()


# In[491]:


# In[476]:


def leastsquare(x,qm,dih_idx, dih_type,dependency, maxorder=6):
    n_dih=len(dih_idx)

    while maxorder*n_dih > 35 :
        maxorder-=1
        print('warning maxorder has been reduced to : %s ' %maxorder)
    # create a set of Parameters
    params = Parameters()
    names=[]
    for i in range(1, maxorder+1):
        for j in range (len(dih_idx)):
            #print(i ,j)
            if dependency.get(j) :
                expr='amp%s_%s'%(i,names[dependency.get(j)])
                #print(expr)
                name= '_'.join(dih_type[dih_idx[j][-1]])+ '_bis'
                #print(name)
                names.append(name)
                params.add('amp%s_%s'%(i,name), value=0 , expr=expr)
            else :
                name='_'.join(dih_type[dih_idx[j][-1]])
                names.append(name)
                params.add('amp%s_%s'%(i,name), value=0)
            # print(name)
            print('order%s_%s'%(i,name))
            params.add('order%s_%s'%(i,name), value=i,vary=False)
    params.add('K', value=0)
    #print(params)
    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(n_dih,maxorder,x, data,names))
    result = minner.minimize()
    return result
    # calculate final result





make_missing_parms_parmck2(sys.argv[1])
all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names, dihedrals_heavy_index , all_dihedrals_type_index = find_dihedrals('input.prmtop')



def make_traj(prefix, index ):
    list_pdb=open(prefix +'-listfiles')
    lines=list_pdb.readlines()
    if len(lines) > 0 : first=lines.pop(0).replace(' \n','')
    else :
        return False
    t = mdtraj.load(first, top='input.prmtop')

    table, bonds = t.topology.to_dataframe()

    for i in range(len(lines)):
        filename=lines.pop(0).replace(' \n','')
        #print(filename)
        tnext= mdtraj.load(filename.replace(' \n',''), top='input.prmtop')
        t= mdtraj.join([t , tnext], check_topology=True, discard_overlapping_frames=False)
    t.save_mdcrd(prefix+ '-traj.mdcrd' )
    t = mdtraj.load(prefix+ '-traj.mdcrd', top='input.prmtop')
    table, bonds = t.topology.to_dataframe()
    #print(table)
    dih= mdtraj.compute_dihedrals(t,index)
    t.save_netcdf(prefix+ '-traj.nc' )
    return t.n_frames ,dih


# In[344]:


def fcn2min(params, n_dih,maxorder,  x, data ,names):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = 0
    for i in range(1,maxorder+1):
        for j in range (0,n_dih):
            amp = params['amp%s_%s'%(i,names[j]) ]
            #shift = params['shift%s'%i]
            order = params['order%s_%s'%(i,names[j]) ]
        #decay = params['decay']

            model +=  amp * np.sin(x[j]*order )
    model+=params['K']
    return model - data





torsion_idx= 0


def make_dependency(torsion_idx,all_dihedrals_type=all_dihedrals_type, all_dihedrals_type_index=all_dihedrals_type_index):
    dependency={}
    count=np.zeros(len( all_dihedrals_type[torsion_idx] ))
    for j in all_dihedrals_type_index[torsion_idx]:
        if count[j[-1]]!= 0 :
            dependency.update({all_dihedrals_type_index[torsion_idx].index(j): int(count[j[-1]] )})
        else :
            count[j[-1]]=all_dihedrals_type_index[torsion_idx].index(j)
    return  dependency




# In[477]:

#x = np.linspace(0, 2*math.pi, 36)
# not true if missing some points... SOLVED!
for torsion in range(len(torsion_names)):
#for torsion in range(4,5):
    print(torsion_names[torsion])
    nframes , x = make_traj(torsion_names[torsion]+'_',[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion] ])

    n_dih=len(all_dihedrals_type_index[torsion])
    print(n_dih)
    data=np.loadtxt(torsion_names[torsion]+'_-energySP.dat')
    k=np.average(data)
    data=data-k
    dependency=make_dependency(torsion)
    #print(all_dihedrals_type[torsion],all_dihedrals_type_index[torsion])
    #print(dependency)
    result=leastsquare(x.T,data,all_dihedrals_type_index[torsion],all_dihedrals_type[torsion],dependency)

    final = data + result.residual

    # write error report
    #report_fit(result)
    plot(data, final,x)
    write_frcmod(result,torsion_names[torsion],  all_dihedrals_type,all_dihedrals_type_index)
add_new_parms(inputfile,torsion_names)


new=parmed.load_file('newtopol.prmtop')

print(parmed.tools.writeFrcmod(new,'new.frcmod').execute())
