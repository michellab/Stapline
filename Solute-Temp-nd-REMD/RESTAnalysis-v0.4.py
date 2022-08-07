#!/usr/bin/env python3


d= '##########################################  \
# This Script mean to group all kind of analysis that can be made on the REST simulations  \
# RMSF on Calpha, MORE OPTIONS to come \
##################'
import os
import sys
import argparse
import mdtraj as md
import copy
import MDAnalysis
import MDAnalysis.analysis.encore as encore
import argparse
from MDAnalysis.analysis import rms
import MDAnalysis.coordinates.TRJ
import numpy as np
from  matplotlib import cm

from sklearn.cluster import KMeans
from matplotlib import dates
#from msmbuilder.cluster import KMeans
import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import matplotlib.pyplot as plt
import math
from sklearn.neighbors.kde import KernelDensity
import random
from multiprocessing import Pool
from mdtraj.formats import PDBTrajectoryFile
import msmbuilder as msm
import msmbuilder.featurizer
import msmbuilder.cluster
import msmbuilder.msm
from itertools import groupby
from pyemma.coordinates import clustering
import mdshare
ipython
import mdshare


def count_event(topfile,trajfile, rep):
    Binary=[0]
    for n in range(rep):

		#trajfile="traj%s.xtc" %(str(n))
		#topfile = 'conf%s.gro'%(str(n))

        traj = md.load(str(n)+'.'+trajfile, top=topfile)

        idxs=traj.topology.select('protein')
        #traj = mdtraj.load(trajfile , atom_indices=idxs,  top=topfile )
        Hy=[]
        Hx=[]
        Hya=[]
        #dssp2 = mdtraj.compute_dssp(traj, simplified=False)
        dssp = md.compute_dssp(traj,simplified=True)
        #addavgHya=[]

        print(str(n)+'.'+trajfile)

        #numpy.set_printoptions(threshold=sys.maxsize)
        #print(dssp)

        count=0

        Hyrunavg=[]


        for i in range(int(len(dssp))):

            unique, counts = np.unique(dssp[i : i+1], return_counts=True)
            #print(unique, counts)
            dictss = { key:val for key, val in dict(zip(unique, counts)).items() if key != 'NA' }
            print(dictss)

            lengthseq= len([j for j in dssp[0]  if j  != 'NA' ])-2
            Hyrunavg.append(dictss.get('H',0)/(lengthseq)*100)
            for val in dictss.values() :
                if np.average(Hyrunavg[-40:]) > 75 and Binary[-1]==0 :
                    Binary.append(1)
                    print( dictss)
                elif np.average(Hyrunavg[-40:])  < 10 and Binary[-1]==1 :
                    Binary.append(0)
                    print( dictss)

                #print( Binary)
            #Hyrunavg.append(dictss.get('E',0)/(lengthseq)*100)
            #print(Binary)
            #Binary = [x[0] for x in groupby(Binary)]

    transition = len(Binary)-1
    count+=transition
    print(count)
    return count

def helicity_per_residus(dssp,n_frames ):
    Hy=np.zeros(len(dssp[0]))
    Ey=np.zeros(len(dssp[0]))
    print(Hy)
    for j in range(len(dssp[0])):
            unique, counts = np.unique(dssp[:,j], return_counts=True)
            #print(unique, counts)
            dictss = { key:val for key, val in dict(zip(unique, counts)).items() if key != 'NA' }
            print(dictss)

            lengthtraj= n_frames
            Hy[j]=dictss.get('H',0)/(lengthtraj)*100
            Ey[j]=dictss.get('E',0)/(lengthtraj)*100
            x=np.arange(len(dssp[0]))
    plt.bar(x, Hy , align='center')
    plt.bar(x, Ey ,align='center', bottom=Hy )

def kde2D(x, y, bandwidth=0.15, xbins=100j, ybins=100j):
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[-math.pi:math.pi:bandwidth , -math.pi:math.pi:bandwidth]


    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)


def Convergence(topolPrefix,  reftopol =False ,reftraj =False ):
     # make a plot computing the distance of the surface of 2 ramachadan plots
     # if ref is given compare it to this simulation , otherwise  compare to the whole simulation

    traj = md.load( '0.striped_traj.nc' ,  top=topolPrefix )
    bandwidth=0.15
    phia=md.compute_phi(traj)[1]
    psia=md.compute_psi(traj)[1]

    if reftopol != 0 :
        print(reftraj, reftopol)
        ref = md.load(  reftraj ,  top=reftopol)
        print(ref)
        phiref = np.concatenate(md.compute_phi(ref)[1][0:])
        psiref = np.concatenate(md.compute_psi(ref)[1][0:])
        xx, yy , zzref =  kde2D(phiref, psiref, xbins=100j, ybins=100j)

    else :
        phi = np.concatenate( phia[0:])
        psi = np.concatenate(psia[0:])

        xx, yy , zzref =  kde2D(phi, psi,  xbins=100j, ybins=100j)
    stride =200

    z = []
    t = []
    '''

    p = Pool(10)
    t = np.arange(50,traj.n_frames  , 120 ) #plot every 2 ns
    band=[bandwidth for i in t]

    args= [(np.concatenate( phia[0:i]),np.concatenate(psia[0:i]) ) for  i in t]

    results = p.starmap(kde2D, args)
    t= 8*t/1000         #  from frames number to ns

    zz=[results[i][2] for i in range(len(t))]

    for zi in zz:
        z.append(np.sqrt(np.sum(np.sum( np.power( zi- zzref,2)) )))
    print(t )
    print(z )
    '''
    started = False
    while z[-1]> 0.05 or started == False :
        started= True
    for i in range(50,traj.n_frames,stride):
        phi = np.concatenate( phia[0:i])
        psi = np.concatenate(psia[0:i])
        xx, yy , zz =  kde2D(phi, psi, bandwidth, xbins=100j, ybins=100j)

        t.append( i)
        z.append(np.sqrt(np.sum(np.sum( np.power( zz- zzref,2)) ))) # to the power to get a sum of positives would work with #np.sum( np.absolute( zz- zzref))
        sys.stdout.write( '\r%s %s' %(t[-1],z[-1]))

    plt.plot(t,z , linewidth=2, markersize=6)
    plt.xlabel('time in ns')
    plt.ylabel('distance in probality (no unit)')
    plt.ylim(top=1)
    plt.savefig('convergence0.png')
    plt.close()
    dz = np.diff(z)/np.diff(t)

    plt.plot(t[:-1],dz , linewidth=2, markersize=6)
    plt.xlabel('time in ns')
    plt.ylabel('d(distance in probality (no unit))/dt')

    plt.savefig('derivconv0.png')
    plt.close()


def runRama(topolPrefix, trajPrefix  , nreps) :
    locator = dates.HourLocator(interval=1)
    locator.MAXTICKS = 1000
    #sns.set()
    mcmap =  cm.get_cmap('gist_ncar_r')
    # https://matplotlib.org/tutorials/colors/colormaps.html  https://python-graph-gallery.com/100-calling-a-color-with-seaborn/
    #gist_heat_r   gist_earth_r   gist_stern    inferno_r')   # gnuplot2')# afmhot_r') #nipy_spectral_r')CMRmap_r  magma_r inferno_r
    mcmap.set_under('w')
    sns.set( style='white' ,font_scale=0.4 , rc={'figure.facecolor':'gray'})  #'darkgrey'
    trajs = []
    #indicesphi = np.array( [[C[residue-1], N[residue-1],CA[residue-1], C[residue]]])
    #indicespsi = np.array( [[ N[residue-1],CA[residue-1], C[residue],  N[residue]] ])
    phis=[]
    psis=[]


    for replica in range(nreps):
        trajs.append(md.load( str(replica ) + '.demux.nc' , stride = 1 ,  top=topolPrefix ))
        traj = trajs[replica]
        phis.append(  md.compute_phi(traj))
        psis.append( md.compute_psi(traj))

    #print(phis)
        #z = np.array(t.n_frames)
        #t = np.array(t.n_frames)
    ramatime=[10, 20, 30, 40 ,  50, 60 , 70 , 80 ,90,  100 , 500 , 600, 700, 800, 900, 1000, 5000 , 10000]
    #ramatime =[15000 , 20000, 25000,  30000,  35000,  40000,  45000 , 50000]
    for laps in ramatime:

        fig, axes = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10), sharex=True, sharey=True)
        cbar_ax = fig.add_axes([0.91, 0.3, 0.03, 0.4])
        #figsnap, axessnap = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10))
        #axes.xaxis.set_minor_locator(locator)
        for repl in range(nreps):

            stride =  1  #laps// 100

                            #CA = t.topology.select('name CA')
            #C = [atom.index for atom in t.topology.atoms if atom.name == 'C']
            #N = [atom.index for atom in t.topology.atoms if atom.name == 'N']
            #CA = [atom.index for atom in t.topology.atoms if atom.name == 'CA']
            #print(N, C , CA)

            #residues = range(2,7)
            residues= range(0,5)
            #print(residues)
            phiall =[]
            psiall =[]
            if phiall : del phiall
            phiall = np.array([])
            if psiall : del psiall
            psiall = np.array([])
            for residue in residues:
                sys.stdout.write('\rPrinting density plot between 0 and %s for replicas %s'%(laps,repl))
                phi=phis[repl][1].T[residue][0:laps]
                psi=psis[repl][1].T[residue][0:laps]




                df = pd.DataFrame( { 'phi' : phi , 'psi' :psi} )
                #print(df.phi, df.psi)
                #sns.kdeplot(df.phi, df.psi, cmap='nipy_spectral_r',  n_levels=1000 , shade=True , ax=axes[replica, residue-2 ]) # norm=LogNorm())

                axes[repl][residue].set(xlim=(-math.pi, math.pi), ylim=(-math.pi,math.pi  ) )
                #ax.xaxis.set_minor_locator(locator)


                #sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[replica][residue-1 ])

                try : sns.kdeplot(df.phi, df.psi,n_levels=1000 ,cmap=mcmap,  shade=True ,  ax=axes[repl][residue] , vmin=0, vmax=1)
                except (RuntimeError, TypeError, Exception) :  sns.kdeplot(df.phi, df.psi,n_levels=100 ,cmap=mcmap, shade=True , ax=axes[repl][residue],  vmin=0, vmax=1)
                phiall=np.concatenate((phiall , phi))
                psiall=np.concatenate((psiall , psi))
                #phiall.append(phi.T[0][:])
                #psiall.append(phi.T[0][:])
                #print(psiall)
            #Overall phi/psi values
            #phi =  md.compute_phi(t , periodic=True, opt=True  )
            #phi =  md.compute_psi(t , periodic=True, opt=True  )

            #phiall = np.array(phiall )
            #psiall= np.array(psiall )
            df = pd.DataFrame( { 'phi' : phiall, 'psi' :psiall  } )
            axes[repl][-1 ].set(xlim=(-math.pi, math.pi), ylim=(-math.pi, math.pi) )
            #print(axes[replica][-1 ])
            im = sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[repl][-1 ] , vmin=0, vmax=1)
        plt.title(laps)

        #fig.tight_layout(rect=[0, 0, .9, 1])


        #fig.colorbar(im, cax=cbar_ax)
        #fig.colorbar()
        #plt.colorbar(fig)
        plt.savefig(str(laps*8)+'ps-rama.jpg')
        plt.close()

def runInstantRama(topolPrefix, trajPrefix , nreps) :
    locator = dates.HourLocator(interval=1)
    locator.MAXTICKS = 1000
    #sns.set()
    mcmap =  cm.get_cmap('gist_ncar_r')
    # https://matplotlib.org/tutorials/colors/colormaps.html  https://python-graph-gallery.com/100-calling-a-color-with-seaborn/
    #gist_heat_r   gist_earth_r   gist_stern    inferno_r')   # gnuplot2')# afmhot_r') #nipy_spectral_r')CMRmap_r  magma_r inferno_r
    mcmap.set_under('w')
    sns.set( style='white' ,font_scale=0.4 , rc={'figure.facecolor':'gray'})  #'darkgrey'
    trajs = []
    #indicesphi = np.array( [[C[residue-1], N[residue-1],CA[residue-1], C[residue]]])
    #indicespsi = np.array( [[ N[residue-1],CA[residue-1], C[residue],  N[residue]] ])
    phis=[]
    psis=[]


    for replica in range(nreps):
        trajs.append(md.load( str(replica ) + '.striped_traj.nc' , stride = 1 ,  top=topolPrefix ))
        traj = trajs[replica]
        phis.append(  md.compute_phi(traj))
        psis.append( md.compute_psi(traj))
        print(len(phis[-1][1]))
    #print(phis)
        #z = np.array(t.n_frames)
        #t = np.array(t.n_frames)

    for laps in range(20,100, 50):

        fig, axes = plt.subplots(ncols=5, nrows=nreps, figsize=(20,10))
        #figsnap, axessnap = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10))
        #axes.xaxis.set_minor_locator(locator)
        for repl in range(nreps):

            stride =  1  #laps// 100
            residues= range(0,5)
            #print(residues)
            phiall =[]
            psiall =[]
            if phiall : del phiall
            phiall = np.array([])
            if psiall : del psiall
            psiall = np.array([])
            for residue in residues:
                sys.stdout.write('\rPrinting instant plot between 0 and %s for replicas %s'%(laps,repl))
                snapphi=phis[repl][1].T[residue][laps-20:laps+20]
                snappsi=psis[repl][1].T[residue][laps-20:laps+20]

                df = pd.DataFrame( { 'phi' :snapphi , 'psi' :snappsi} )
                #print(df.phi, df.psi)
                #sns.kdeplot(df.phi, df.psi, cmap='nipy_spectral_r',  n_levels=1000 , shade=True , ax=axes[replica, residue-2 ]) # norm=LogNorm())

                axes[repl][residue].set(xlim=(-math.pi, math.pi), ylim=(-math.pi,math.pi  ) )

                #ax.xaxis.set_minor_locator(locator)


                #sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[replica][residue-1 ])

                axes[repl][residue].plot(df.phi, df.psi, color='green', marker='o', linestyle='dashed', linewidth=1, markersize=4)

                #phiall=np.concatenate((phiall , snapphi))
                #psiall=np.concatenate((psiall ,snappsi))

            #df = pd.DataFrame( { 'phi' : phiall, 'psi' :psiall  } )
            #axes[repl][-1 ].set(xlim=(-math.pi, math.pi), ylim=(-math.pi, math.pi) )
            #print(axes[replica][-1 ])
            #axes[repl][-1 ].plot(df.phi, df.psi, 'r+')

        plt.title(laps)
        #plt.colorbar(fig)ashrc
        plt.savefig(str(laps)+'-instrama.jpg')
        plt.close()



def runCluster_pyemma(t) :
    feat = pyemma.coordinates.featurizer('system.wat.leap.strip.prmtop')
    feat.add_backbone_torsions()
    data = pyemma.coordinates.load('0.striped_traj.nc' ,features=feat)

    nb_cluster =math.floor(min(math.sqrt(math.sqrt(t.n_frames)), 50 ))
    clu_method=clustering.KmeansClustering(nb_cluster, max_iter=5, metric='euclidean', tolerance=1e-05, init_strategy='kmeans++', fixed_seed=False, oom_strategy='memmap', stride=1, n_jobs=None, skip=0, clustercenters=None, keep_data=False)
    clusters=clu_method.fit(data)
    for i in range(nb_cluster):
        sliced_all=t.slice([ j for j in clusters[i]])
        sliced_all.superpose(sliced_all,frame=0)
        filename='clusterpyemma%s.nc' %(i)
        sliced_all.save_netcdf(filename)
    return nb_clusters, clusters


def runCluster(t ):

        #universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
        #cluster_collection = encore.cluster(universe , method= KMeans(100,n_jobs=5 ))

        dataset = []
        noTFE=t.topology.select("not resname TFE")


        # compute over psi phi and chi1
        phi = md.compute_phi(t[::1])
        psi = md.compute_psi(t[::1])
        #chi1 = md.compute_chi1(t[::1])

        ## remove the first 2 and 2 last residus
        metric = np.concatenate(( phi[1], psi[1] ) , axis = 1 )
        #list to keep initial frame numbers as metric values will be deleted as they are used
        # Only usefull if the snapshot are output
        indexes= np.arange(0,t.n_frames)


        cutoff = 5

        nb_cluster=0
        list_old=[]
        while  metric.shape[0] !=0 :

            nb_cluster+=0
            #medoid = find_medoid(metric)
            #Not really a medoid here but a random frame optfully representative of the most common conformation
            # ! the medoid will not be reajusted
            medoid = metric[random.randint(0,len(metric)-1)]

            #for structure in metric.shape[0] :
            diff = np.mod( metric - medoid +math.pi , 2* math.pi )  -math.pi
            #dist = np.sum(metric[])
        #for i in range(len(metric2)): metric2[i] = np.reshape(metric[i] , (5, 2), order='F')


            dist = np.sqrt(np.sum(np.power(diff,2) , axis= 1 ))
            '''
            weight=np.concatenate([np.linspace(0, 5, num=4),np.linspace( 5,0, num=4)])

            weight=np.ones(diff.shape[1])
            print(weight)
            print(np.average(np.power(diff,2) , axis= 1 , weights=weight, returned=True  ))
            dist = np.sqrt(np.average(np.power(diff,2) , axis= 1 , weights=weight, returned=True  )[0])
            '''
            list_in= np.array([])

            for i in range (len(dist)):
                #print(dist[i])
                if dist[i] < cutoff:
                    #metric = np.delete(metric, i, 0)

                    list_in = np.append(list_in , [int(i)] )
            #print(list_in)
            if nb_cluster == 0 :
                #clusters = np.array(list , dtype=object)
                clusters = []



            nb_cluster += 1

            i=0
            sliced_all = False

            #print(t.n_frames, len(indexes), indexes[-1])

            sliced_all=t.slice([int(indexes[int(j)]) for j in list_in])
            sliced_all.superpose(sliced_all,frame=0)
            filename='cluster%s.nc' %(nb_cluster-1)
            filenamepdb='cluster%s.pdb' %(nb_cluster-1)
            sliced_all.save_netcdf(filename)
            #sliced_all.save_pdb(filenamepdb)
            #del sliced_all


            clusters.append([int(indexes[int(j)]) for j in list_in])
            for j in list_in[::-1]:
                j=int(j)
                metric = np.delete(metric, j, 0)
                indexes = np.delete(indexes, j, 0)
                #filename='cluster%s_snap%s.pdb' %(nb_cluster, indexes[j])

                '''
                #sliced.save_pdb(filename)
                if sliced_all == False :
                    sliced_all=sliced
                else :  sliced_all.join( [sliced_all,sliced])
                print(sliced_all)
                '''

        np.set_printoptions(threshold=sys.maxsize)
        trajCLUSTERED=np.empty((t.n_frames),dtype=int)
        for i in range(len(clusters)):
            for j in clusters[i]:
                trajCLUSTERED[j]=i
        #print(trajCLUSTERED)



        lmax=10
        markov=np.zeros((lmax,nb_cluster,nb_cluster),dtype=int)


        for l in range(1,lmax):
            for i in range(len(clusters)):
                for j in  range(len(clusters[i])):

                    #if len([f  for f in  range(len(clusters)) if  clusters[i][j]-1  in  clusters[f]) >1 :
                            #sys.exit('bug')

                    if   clusters[i][j]-l >= 0 and len([f  for f in  range(len(clusters)) if  clusters[i][j]-l !=0 and clusters[i][j]-1  in  clusters[f]])==1 :
                        index=  [f  for f in  range(len(clusters)) if  clusters[i][j]-l !=0 and clusters[i][j]-l  in  clusters[f]][0]
                        markov[l,i,index]+=1
                    elif len([f  for f in  range(len(clusters)) if  clusters[i][j]-l !=0 and clusters[i][j]-l  in  clusters[f]]) >1 or  clusters[i][j]-l !=0 :
                        #print(clusters[i][j] )
                        index=0
        #keep biggest 10 clusters

        #print(markov.shape)
        for i in range(min(markov.shape[1],10) ) :
            for j in range(min(markov.shape[1],10) )  :
                plt.plot(markov[:,i,j])
        plt.savefig('markov.png' )
        plt.close()
        print(markov)

        return  nb_cluster,  clusters

def find_correlation(nb_cluster, rep=1):
    j=1

    plt.figure(figsize=(14,14))
    nb_cluster=nb_cluster[:15]

    for i in range(len(nb_cluster)):

        traj=md.load('cluster%s_%s.nc' %(rep,i) , top='system.wat.leap.strip.prmtop' )
        dssp = md.compute_dssp(traj)
        #stapled_residues= traj.topology.select("resname  AS1 or resname  AS2 or resname AR1 or resname AR2  ")
        #stapled_residues= traj.atom_slice('resname in AS1 AS2 AR1 AR2')
        #stapled_residues= traj.atom_slice(traj.topology.select('resname  AS* or resname AR*'))
        stapled_residues= traj.atom_slice([atom.index for atom in traj.topology.atoms if atom.residue.name in  ( 'AS1' ,'AR1', 'AS2' ,'AR2' ,'AKS' , 'AKR' ,'A1R' ,'A1S' ,'A2R' ,'A2S'  )])
        #print(stapled_residues)
        #stapled_residues=[residue.index for residue in traj.topology.chain(0).residues if residue.name in  ( 'AS1' ,'AR1', 'AS2' ,'AR2'   )]
        #print([atom.index for atom in stapled_residues.topology.atoms if atom.name = 'CA' ])

        atoms1 = ['CA', 'CB', 'CC', 'CD' , 'CE', 'CF','CY' ]
        atoms2 =  ['CA', 'CB', 'CY' ]
        atoms2.reverse()
        atoms=copy.deepcopy( atoms1 )
        atoms+= atoms2
        stapled_residue_index = np.zeros(len(atoms), dtype=int)
        #print(stapled_residues.topology.atoms)

        for atom in stapled_residues.topology.atoms :

            if atom.name in atoms1:

                index=atoms1.index(atom.name)

                if stapled_residue_index[index]==0 :
                    stapled_residue_index[index]=atom.residue.index

                elif atom.name in atoms2:
                    index=atoms2.index(atom.name)

                    stapled_residue_index[-index-1]=atom.residue.index

        if 0 in  stapled_residue_index  : sys.exit('check your residu name and atoms names')
        #'CY', 'CD',   'CC','CB', 'CA'

        chi1_atoms= np.array([[atom.index for atom in stapled_residues.topology.atoms if atom.name == 'N' and atom.residue.index==stapled_residue_index[0] or  atom.name == 'N' and atom.residue.index==stapled_residue_index[-1]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[1] and atom.residue.index==stapled_residue_index[0] or  atom.name == atoms[-1] and atom.residue.index==stapled_residue_index[-1]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[2] and atom.residue.index==stapled_residue_index[1] or  atom.name == atoms[-2] and atom.residue.index==stapled_residue_index[-2]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[3] and atom.residue.index==stapled_residue_index[2] or  atom.name == atoms[-3] and atom.residue.index==stapled_residue_index[-3]] ]).T


        chi2_atoms= np.array([[atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[0] and atom.residue.index==stapled_residue_index[0] or  atom.name == atoms[-1] and atom.residue.index==stapled_residue_index[-1]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[1] and atom.residue.index==stapled_residue_index[1] or  atom.name == atoms[-2] and atom.residue.index==stapled_residue_index[-2]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[2] and atom.residue.index==stapled_residue_index[2] or  atom.name == atoms[-3] and atom.residue.index==stapled_residue_index[-3]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[3] and atom.residue.index==stapled_residue_index[3] or  atom.name == atoms[-4] and atom.residue.index==stapled_residue_index[-4]] ]).T

        chi3_atoms= np.array([[atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[1] and atom.residue.index==stapled_residue_index[1] or  atom.name == atoms[-2] and atom.residue.index==stapled_residue_index[-2]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[2] and atom.residue.index==stapled_residue_index[2] or  atom.name == atoms[-3] and atom.residue.index==stapled_residue_index[-3]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[3] and atom.residue.index==stapled_residue_index[3] or  atom.name == atoms[-4] and atom.residue.index==stapled_residue_index[-4]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[4] and atom.residue.index==stapled_residue_index[4] or  atom.name == atoms[-5] and atom.residue.index==stapled_residue_index[-5]] ]).T

        chi4_atoms= np.array([[atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[2] and atom.residue.index==stapled_residue_index[2] or  atom.name == atoms[-3] and atom.residue.index==stapled_residue_index[-3]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[3] and atom.residue.index==stapled_residue_index[3] or  atom.name == atoms[-4] and atom.residue.index==stapled_residue_index[-4]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[4] and atom.residue.index==stapled_residue_index[4] or  atom.name == atoms[-5] and atom.residue.index==stapled_residue_index[-5]] , \
                              [atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[5] and atom.residue.index==stapled_residue_index[5] or  atom.name == atoms[-6] and atom.residue.index==stapled_residue_index[-6]] ]).T



        chi1 =  md.compute_dihedrals(stapled_residues, chi1_atoms)
        chi2 =  md.compute_dihedrals(stapled_residues, chi2_atoms)
        chi3 =  md.compute_dihedrals(stapled_residues, chi3_atoms)
        chi4 =  md.compute_dihedrals(stapled_residues, chi4_atoms)

        CatoCA= md.compute_distances(stapled_residues, np.array([[atom.index for atom in stapled_residues.topology.atoms if atom.name == 'CA']]), periodic=True, opt=True)
        print(CatoCA)
        plt.subplot2grid( (len(nb_cluster) , 5), (i,2)  )
        plt.ylim(0.55,0.85)
        plt.plot(CatoCA ,label = dssp[-1])
        plt.subplot2grid( (len(nb_cluster) , 5), (i,3)  )
        helicity_per_residus(dssp,traj.n_frames)
        plt.ylim(0,50)



        vectors_index =np.array([[[  atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[0] and atom.residue.index==stapled_residue_index[0] ][0] , \
                            [ atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[1] and atom.residue.index==stapled_residue_index[1]][0] ], \
                           [[ atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[0] and atom.residue.index==stapled_residue_index[4]][0] , \
                            [ atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[4] and atom.residue.index==stapled_residue_index[5] ][0]], \
                           [[ atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[-1] and atom.residue.index==stapled_residue_index[-1] ][0] , \
                            [  atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[-2] and atom.residue.index==stapled_residue_index[-2] ][0]] ])
        angle_index =np.array([[[  atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[0] and atom.residue.index==stapled_residue_index[0] ][0] , \
                            [ atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[1] and atom.residue.index==stapled_residue_index[1]][0] ], \
                            [ atom.index for atom in stapled_residues.topology.atoms if atom.name == atoms[4] and atom.residue.index==stapled_residue_index[5] ][0] ])


        vectors=md.compute_displacements(stapled_residues ,vectors_index)
        angle=md.compute_angles(stapled_residues ,angle_index)
        """ Returns the angle in radians between vectors   """
        cosang = [ np.dot(vectors[i,0], vectors[i,0])  for i in range(traj.n_frames)]
        sinang =  [ np.linalg.norm(  np.cross(vectors[i,0], vectors[i,1]) )   for i in range(traj.n_frames)]
        direction_staple  = np.arctan2(sinang, cosang)
        plt.subplot2grid( (len(nb_cluster) , 5), (i,4)  )
        plt.hist(direction_staple, bins=100,    label = 'direction of the staple')
        plt.hist(angle[i,0], bins=100,    label = 'direction of the staple')
        plt.xlim(0,math.pi)

        for h in range(2) :
            plt.subplot2grid( (len(nb_cluster) , 5), (i,h)  )
            plt.hist(chi1[:,h], bins=100,    label = 'chi1%s'  %(h))
            plt.hist(chi2[:,h], bins=100,  label = 'chi2%s'%(h))
            plt.hist(chi3[:,h], bins=100,  label = 'chi3%s'%(h))
            #plt.suptitle('cluster%s   ' %(i))
            #plt.title('dssp ; %s' %(dssp[-1]), fontsize=5)
            j+=1


    plt.legend()
    plt.savefig('rep%s.png' %(rep))
    plt.close()




def find_medoid(metric) :
    #for the moment that is not a metroid but a sort of centroid
    avgs =[]
    length = metric.shape[0]
    #for i in range(10):
            #avgs.append( np.sum( metric2[:, i])/length)
    avgs=np.median(metric, axis=0)
    return avgs





    '''
        for frame in range(len(metric)):
            for res in range(len(metric[frame]):
                 metric[frame][res] = [ phi[1][frame][res] ,  psi[1][frame][res] ]

    '''

        #kmeans = KMeans(n_clusters=10, random_state=0)


        #kmeans.cluster_centers_

        #cluster = KMeans(n_clusters=10)
        #cluster.fit(dataset)

        #print( kmeans)

def markovbuilder( topolPrefix, trajPrefix , i,relevantcluster):
        dataset = []

        t = md.load( str(i) + '.striped_traj.nc' , top='system.wat.leap.strip.prmtop' )



        # compute over psi phi and chi1
        #phi = md.compute_phi(t[::1])
        #psi = md.compute_psi(t[::1])
        #chi1 = md.compute_chi1(t[::1])

        ## remove the first 2 and 2 last residus
        #metric = np.concatenate(( phi[1], psi[1] ) , axis = 1 )
        feat = msm.featurizer.DihedralFeaturizer(types=['phi', 'psi'], sincos=True)
        descriptor = feat.partial_transform(t)
        #list to keep initial frame numbers as metric values will be deleted as they are used
        # Only usefull if the snapshot are output
        indexes= np.arange(0,t.n_frames)
        algo = msm.cluster.KMedoids(n_clusters=relevantcluster, metric='euclidean')
        new=algo.transform(descriptor)[0]
        markov = np.zeros((relevantcluster,relevantcluster,10), dtype=int)
        for i in range(10) :
            for j in range(i, t.n_frames):
                print(new[j], new[j-i])
                markov[new[j], new[j-i], i]+=1





def runCluster2(topolPrefix , trajPrefix, nreps, ncluster):
    for i in range(nreps) :
        #universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
        #cluster_collection = encore.cluster(universe , method= KMeans(100,n_jobs=5 ))
        dataset = []

        t = md.load( str(i) + '.demux.nc' , top=topolPrefix )
        # pick  nclusters random structures
        centroids=np.zeros( ncluster)
        clusterframes=[]
        rmsd2centroid = []
        for i in range(ncluster):
            centroids[i] = (t.n_frames/  ncluster )*i
            rmsd2centroid.append([])
            clusterframes.append([])
        # sart :
        for conv in range(10) :   # ten for now but
            # will have to find a way to quantify the convergenge something like :
            # sum ((new centroid - prev. centroid)^2 ) < some_value
            for i in range(ncluster):

                rmsd2centroid [i] = md.rmsd(t, t , frame=centroids[i],  parallel=True, precentered=False).T
            j = range(t.n_frames)

            for i in range(ncluster) :
                del clusterframes[i]

            while j != 0:
                for i in range(ncluster) :
                    if j != 0 :
                        idx = rmsd+centroid[i].index( min(rmsd+centroid[i] ))
                        for k in range(ncluster) :
                            rmsd+centroid[k].pop(index)
                        j.pop(index)
                        clusterframes[i].append(index)

            #find the new centroid of each cluster
            for i in range(ncluster) :
                distances = np.empty((len(clusterframes[i]), clusterframes[i]))
                for j in range(clusterframes[i]):
                    distances[i] = md.rmsd(traj, traj, j, atom_indices=atom_indices)

        print( kmeans)




def diffusionAnalysis(nreps):
    f = open('rem.log','r')
    started =0
    out=[0]
    prev=[0]
    caption = '#exchange, '
    for i in range(1,nreps+1):
        out.append(i)
        prev.append(i)
        caption += 'traj'+str(i)
    # position of ham0
    index0 = 1
    # output file
    outfile= open ('outdemux.dat' ,'w')

    outfile.writelines( caption +  'position of ham0')
    # variable if the replica goes from 0 to the higest up=1 otherwise for down up =0
    Up=1
    # count the number of up and down of repl 1
    compteur = 0
    for line in f.readlines():
        #print(started , prev ,out, int(line.split()[0]) , out[int(line.split()[0])], prev[0])
        if line[:10] != '# exchange':
            if started ==0 :
                continue
            else:
                if line.split()[7]== 'F':

                    out[int(line.split()[0])]=prev[int(line.split()[0])]
                else :
                    out[int(line.split()[0])]=prev[int(line.split()[1])]
        else  :
                if started ==1 :
                    index0=out.index(1)
                    if Up== 1 and index0==nreps:
                        Up=0
                    elif Up==0 and index0==1 :
                        Up=1
                        compteur += 1
                        if compteur ==1 : print('it took %s exchanges attempts for the replica 1 to climb the ladder up and down the first time ' %out[0])

                    lineout=''
                    for j in out :
                         lineout+=  str(j) + ' '
                    lineout+= '            ' + str(index0)
                    outfile.writelines( lineout)
                    prev = copy.deepcopy(out)

                    out[0] =f'{int(line.split()[2]):05}'#print("{:02d}".format(1))
                else :
                    started =1
                    out[0] = '00000'

    f.close()
    outfile.close()
    print('Overall, %s exchanges were attempted and the replicas 1 climbed the hamiltonian ladder up and down %s times ' %(str(out[0]), str(compteur)))



def diffusionAnalysis_bck(nreps,logfile,topolPrefix, trajPrefix,relevantcluster ):
    f = open('rem.log','r')
    started =0
    out=[0]
    prev=[0]
    caption = '#exchange, '
    for i in range(1,nreps+1):
        out.append(i)
        prev.append(i)
        caption += 'traj'+str(i) + ', '
    # position of ham0
    index0 = 1
    # output file
    outfile= open ('outdemux.dat' ,'w')

    outfile.writelines( caption +  'position of ham0')
    # variable if the replica goes from 0 to the higest up=1 otherwise for down up =0
    Up=1
    # count the number of up and down of repl 1
    compteur = 0
    for line in f.readlines():

        #print(started , prev ,out, int(line.split()[0]) , out[int(line.split()[0])], prev[0])
        if line[:10] != '# exchange':

            if started ==0 :

                continue
            else:
                if line.split()[7]== 'F':

                    out[int(line.split()[0])]=prev[int(line.split()[0])]
                    print(line)
                else :
                    out[int(line.split()[0])]=prev[int(line.split()[1])]

        else  :

                if started ==1 :
                    index0=out.index(1)
                    if Up== 1 and index0==nreps:
                        Up=0
                    elif Up==0 and index0==1 :
                        #Up=markovbuilder( topolPrefix, trajPrefix , i-1,relevantcluster)
                        compteur += 1
                        if compteur ==1 :
                            comment='it took %s exchanges attempts for the replica 1 to climb the ladder up and down the first time \n' %out[0]
                            print(comment)
                            logfile.write(comment)

                    lineout=''
                    for j in out :
                         print(lineout)
                         lineout+=  str(j) + ' '
                    lineout+= '            ' + str(index0)
                    outfile.writelines( lineout + '\n')
                    prev = copy.deepcopy(out)

                    out[0] =f'{int(line.split()[2]):05}'#print("{:02d}".format(1))
                else :
                    #started =markovbuilder( topolPrefix, trajPrefix , i,relevantcluster)
                    out[0] = '00000'

    f.close()
    outfile.close()
    comment='Overall, %s exchanges were attempted and the replicas 1 climbed the hamiltonian ladder up and down %s times \n ' %(str(out[0]), str(compteur))
    print(comment)
    logfile.write(comment)


def demuxTraj(trajPrefix, topolPrefix , nreps):
 cpptrajInputFile1=open('cpptrajDemux.in' ,'w')
 cpptrajInputFile2=open('cpptrajStrip.in' ,'w')
 cpptrajInput='parm %s \nensemble %s_0.nc trajnames ' %(topolPrefix , trajPrefix)
 for i in range(1 ,nreps)  :
  cpptrajInput+='%s_%s.nc,' %(trajPrefix ,i)
 cpptrajInput1 = cpptrajInput[:-1] +' remlog rem.log nstlim 25000 ntwx 5000\nautoimage\nstrip :WAT,TFE \ntrajout  demux.nc \nparm %s name top \nparmstrip :WAT,TFE parm top \nparmwrite parm top out system.wat.leap.strip.prmtop\nrun\n'%topolPrefix
 #cpptrajInput1 = cpptrajInput[:-1] +' remlog rem.log\nautoimage\nstrip :WAT,TFE \ntrajout  demux.nc \nrun\nparm %s name top \nparmstrip :WAT,TFE parm top \nparmwrite parm top out system.wat.leap.strip.prmtop\nrun\n'%topolPrefix

 cpptrajInput2 = cpptrajInput[:-1] +' nosort remlog rem.log \nautoimage\nstrip :WAT,TFE \ntrajout  striped_traj.nc \nparm %s name top \nparmstrip :WAT,TFE parm top \nparmwrite parm top out system.wat.leap.strip.prmtop\nrun\n'%topolPrefix
 cpptrajInputFile1.writelines(cpptrajInput2)
 cpptrajInputFile2.writelines(cpptrajInput1)
 cpptrajInputFile1.close()
 cpptrajInputFile2.close()
 os.system('cpptraj -i cpptrajDemux.in' )
 os.system('cpptraj -i cpptrajStrip.in' )
 for i in range(nreps)  :
  try : os.rename('demux.nc.' +str(i) , str(i)+'.demux.nc')
  except :pass
  os.rename('striped_traj.nc.' +str(i) , str(i)+'.striped_traj.nc')
def runRMSF(topolPrefix, trajPrefix, nreps, atommask, logfile):
 for i in range(nreps) :
  repl= md.load(str(i)+'.demux.nc',top=topolPrefix )
  atom_indices= repl.topology.select('name  %s' %(atommask))   # [atom.index for atom.name== 'CA' in traj.topology.repl.
  #rmsf = md.rmsf(repl,   frame=0, atom_indices=atom_indices, parallel=True, precentered=False)
  #print(rmsf)

  repl.superpose(repl, atom_indices=atom_indices, ref_atom_indices=atom_indices)
  avg_xyz = np.mean(repl.xyz[:, atom_indices, :], axis=0)
  rmsf = np.sqrt(3*np.mean((repl.xyz[:, atom_indices, :] - avg_xyz)**2, axis=(0,2)))


  '''
  print(topolPrefix , str(i) + '.demux.nc')
  universe = MDAnalysis.Universlogfile.write(comment)e(topolPrefix , str(i) + '.demux.nc' )
  universe = MDAnalysis.Universe()
  proteine=universe.select_atoms("protein and name CA")
  R = rms.RMSF(protein)
  R.run()
  ax = plt.subplot(111)
  '''
  plt.plot( atom_indices, rmsf,  linewidth=1 ,label='repl %s'%(i))
  #plt.fill_between(  atom_indices,  rmsf, color="red", alpha=0.1)
  #sns.despine(ax=ax, offset=10)
  #plt.set_ylabel(r"RMSF ($\AA$)")
 #plt.set_xlabel("atom number")
  comment='for %s the rmsf average for %s is %s and the standard deviation is %s'  %(i ,atommask, np.average(rmsf) ,  np.std(rmsf) )
  print(comment)
  logfile.writelines(comment +'\n')
 plt.legend()
 plt.savefig('rmsf%s.png' %(atommask) )
 plt.close()




def main() :
 print('starting script')
 parser = argparse.ArgumentParser(description=d)

 parser.add_argument('--nreps', type=int,  help='number of replicas', required=True)
 parser.add_argument('-t', '--topol',  default='system.wat.leap.prmtop', help='topology')
 parser.add_argument('-xyz', '--trajPrefix',  default='traj', help='traj_prfix')
 parser.add_argument('--rewrite',  action='store_true',help='Force the script ro rewrite the trajectories if those have alredy been found')
 parser.add_argument('--rmsf', action='store_true',  help='rmsf')
 parser.add_argument('--cluster',  action='store_true',  default=False ,help='cluster trajectory')
 parser.add_argument('--corr',  action='store_true',  default=False ,help='atempt to find correlations with CHI1, CHI2 . Also calculate Ca-Ca distances')
 parser.add_argument('--rama', action='store_true', help='ramchadan')
 parser.add_argument('--instrama', action='store_true', help='ramchadan')
 parser.add_argument('--exch', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--conv', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--c_event', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--toreftop',help='compare ramachadan to a converged simulation')
 parser.add_argument('--torefxyz',help='compare ramachadan to a converged simulation')

 args = parser.parse_args()
 print(args)

 nreps =  args.nreps
 topolPrefix = args.topol
 trajPrefix = args.trajPrefix
 rewrite = args.rewrite
 RMSF = args.rmsf
 #toref = args.toref
 exch = args.exch
 conv = args.conv
 corr=args.corr
 rama = args.rama
 cluster= args.cluster
 ncluster = 10
 c_event=args.c_event
 if os.path.isfile('./0.demux.nc') == False or args.rewrite==True :
    demuxTraj(trajPrefix,topolPrefix, nreps)
 topolPrefix = 'system.wat.leap.strip.prmtop'

 logfile=open('RESTanalysis.log','a')


 relevantcluster=0
 if RMSF == True:
     runRMSF( topolPrefix, trajPrefix , nreps, 'CA',logfile)
     runRMSF( topolPrefix, trajPrefix , nreps, 'CB',logfile)
 if rama == True:
    runRama( topolPrefix, trajPrefix , nreps )
 if args.instrama == True:
    runInstantRama( topolPrefix, trajPrefix , nreps)
 if  cluster == True or corr == True :
    traj = md.join([md.load( str(rep) + '.striped_traj.nc' , top='system.wat.leap.strip.prmtop' ) for rep in range(nreps) ] , check_topology=True, discard_overlapping_frames=False)
    #c, nb_cluster = runCluster(  traj  )
    c, nb_cluster = runCluster_pyemma(  traj  )
    comment= 'There was ' + str(c) +'clusters'

    print( comment)
    logfile.writelines(comment + '\n')
    relevantcluster=2
    for j in nb_cluster:
            if len(j )> 10 :
                relevantcluster+=1
    if  corr== True: find_correlation(nb_cluster, i)
 if  c_event== True:

        count=count_event( 'system.wat.leap.strip.prmtop','demux.nc', nreps)



 if  exch ==True :
     diffusionAnalysis(nreps )
     """markovbuilder( topolPrefix, trajPrefix , i,relevantcluster)"""
 count=0

 #print(count)
 #logfile.write('We observed %s transitions \n' %(count))
 if args.conv == True:
    Convergence(topolPrefix,args.toreftop,args.torefxyz )
 logfile.close()
if __name__ == '__main__' :

    main()
