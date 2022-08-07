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

def kde2D(x, y, bandwidth, xbins=100j, ybins=100j):
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[-math.pi:math.pi:bandwidth/xbins, -math.pi:math.pi:bandwidth/ybins]

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
    traj = md.load( '0.demux.nc' ,  top=topolPrefix )
    bandwidth=0.5

    if reftopol != 0 :
        ref = md.load(  reftraj ,  top=reftopol)
        phiref = np.concatenate(md.compute_phi(ref)[1][0:])
        psiref = np.concatenate(md.compute_psi(ref)[1][0:])
        xx, yy , zzref =  kde2D(phiref, psiref, bandwidth, xbins=100j, ybins=100j)

    else :
        phi = np.concatenate( md.compute_phi(traj)[1][0:])
        psi = np.concatenate(md.compute_psi(traj)[1][0:])
        print(len(psi))
        xx, yy , zzref =  kde2D(phi, psi, bandwidth, xbins=100j, ybins=100j)
    stride =100
    z = np.array(t.n_frames)
    t = np.array(t.n_frames)
    for i in range(100,t.n_frames,stride):
        phi = np.concatenate( md.compute_phi(traj)[1][0:i])
        psi = np.concatenate(md.compute_psi(traj)[1][0:i])
        xx, yy , zz =  kde2D(phi, psi, bandwidth, xbins=100j, ybins=100j)
        t[i]= i
        z[i]= np.sum( np.power( zz- zzref,2))  # to the power to get a sum of positives would work with #np.sum( np.absolute( zz- zzref))

    plt.plot(t,z , linewidth=2, markersize=12)





def runRama(topolPrefix, trajPrefix , nreps) :
    locator = dates.HourLocator(interval=1)
    locator.MAXTICKS = 1000
    #sns.set()
    mcmap =  cm.get_cmap('gist_ncar_r')
    # https://matplotlib.org/tutorials/colors/colormaps.html  https://python-graph-gallery.com/100-calling-a-color-with-seaborn/
    #gist_heat_r   gist_earth_r   gist_stern    inferno_r')   # gnuplot2')# afmhot_r') #nipy_spectral_r')CMRmap_r  magma_r inferno_r
    mcmap.set_under('w')
    sns.set( style='white' ,font_scale=0.4 , rc={'figure.facecolor':'gray'})  #'darkgrey'
    trajs = []
    for replica in range(nreps):
        trajs.append(md.load( str(replica ) + '.demux.nc' , stride = 1 ,  top=topolPrefix ))
    for laps in range(8000, 50000 , 200):
        #os.system( 'rm phipsiall-%s.dat'  %(  str(laps)))

        fig, axes = plt.subplots(ncols=6, nrows=nreps, figsize=(20,10))
        #axes.xaxis.set_minor_locator(locator)

        for replica in range(nreps):
            traj = trajs[replica]
            stride =  1  #laps// 100
            t = traj[0:laps:stride]
                            #CA = t.topology.select('name CA')
            C = [atom.index for atom in t.topology.atoms if atom.name == 'C']
            N = [atom.index for atom in t.topology.atoms if atom.name == 'N']
            CA = [atom.index for atom in t.topology.atoms if atom.name == 'CA']
            #print(N, C , CA)

            #residues = range(2,7)
            residues= range(1,6)
            #print(residues)
            phiall =[]
            psiall =[]
            if phiall : del phiall
            phiall = np.array([])
            if psiall : del psiall
            psiall = np.array([])
            for residue in residues:
                print( 'residue : %s \n replica : %s' %(residue, replica))
                ''' Using CPPtraj --> noo output is writen :-(
                cpptrajin = open('rama' + str(residue) + str(laps), 'w')
                cpptrajin.writelines( ' \
trajin %s.demux.nc   1 %s\n \
dihedral phi :%s@C :%s@N :%s@CA :%s@C out phipsi%s-%s-%s.dat\n\
dihedral psi :%s@N :%s@CA :%s@C :%s@N out phipsi%s-%s-%s.dat \n\
run  \n \
quit\n'  \
%(str(replica), str(laps) ,\
str(residue) , str(residue + 1 )  ,str(residue + 1 )  , str(residue + 1 )  ,str(replica) , str(residue),  str(laps),\
str(residue + 1 )  ,str(residue + 1 )  , str(residue + 1 ) , str(residue)  ,str(replica) , str(residue),  str(laps) ))
                print('running : \n cpptraj %s -i rama%s%s ' %(topolPrefix, str(residue) , str(laps)) )
                #os.system('echo $AMBERHOME')

                os.system('source $AMBERHOME/amber.sh ; cpptraj  %s -i rama%s%s ' %(topolPrefix, str(residue) , str(laps)))
                filename= 'phipsi%s-%s-%s.dat'   %(str(replica) , str(residue),  str(laps))
                os.system('cat %s >> phipsiall-%s.dat'  %(filename ,  str(laps)) )
                t ,phi ,psi  =np.loadtxt(filename , comments='#' , unpack=True )
                '''

                #C = t.topology.select('name C')
                #N = t.topology.select('name N')



                indicesphi = np.array( [[C[residue-1], N[residue-1],CA[residue-1], C[residue]]])
                indicespsi = np.array( [[ N[residue-1],CA[residue-1], C[residue],  N[residue]] ])

                phi= md.compute_dihedrals(t, indicesphi , periodic=True, opt=True  )
                psi= md.compute_dihedrals(t,indicespsi, periodic=True, opt=True )


                df = pd.DataFrame( { 'phi' : phi.T[0] , 'psi' :psi.T[0] } )
                #print(df.phi, df.psi)
                #sns.kdeplot(df.phi, df.psi, cmap='nipy_spectral_r',  n_levels=1000 , shade=True , ax=axes[replica, residue-2 ]) # norm=LogNorm())

                axes[replica][residue-1 ].set(xlim=(-math.pi, math.pi), ylim=(-math.pi,math.pi  ) )
                #ax.xaxis.set_minor_locator(locator)


                #sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[replica][residue-1 ])

                try : sns.kdeplot(df.phi, df.psi,n_levels=1000 ,cmap=mcmap,  shade=True , ax=axes[replica][residue-1 ])
                except (RuntimeError, TypeError, Exception) :  sns.kdeplot(df.phi, df.psi,n_levels=10 ,cmap=mcmap,  shade=True , ax=axes[replica][residue-1 ])

                #plt.show()
            	#'Spectral_r'
                phiall=np.concatenate((phiall , phi.T[0]))
                psiall=np.concatenate((psiall , psi.T[0]))
                #phiall.append(phi.T[0][:])
                #psiall.append(phi.T[0][:])
                #print(psiall)
            #Overall phi/psi values
            #phi =  md.compute_phi(t , periodic=True, opt=True  )
            #phi =  md.compute_psi(t , periodic=True, opt=True  )

            #phiall = np.array(phiall )
            #psiall= np.array(psiall )
            df = pd.DataFrame( { 'phi' : phiall, 'psi' :psiall  } )
            axes[replica][-1 ].set(xlim=(-math.pi, math.pi), ylim=(-math.pi, math.pi) )
            #print(axes[replica][-1 ])
            sns.kdeplot(df.phi, df.psi, cmap=mcmap,  n_levels=1000 , shade=True , ax=axes[replica][-1 ])
        plt.title(laps)
        #plt.colorbar(fig)
        plt.savefig(str(laps)+'-rama.jpg')





def runCluster(topolPrefix , trajPrefix, nreps):
    #for i in range(nreps) :
        #universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
        #cluster_collection = encore.cluster(universe , method= KMeans(100,n_jobs=5 ))

        dataset = []

        t = md.load( str(0) + '.demux.nc' , top='system.wat.leap.strip.prmtop' )
        phi = md.compute_phi(t[::1])
        psi = md.compute_psi(t[::1])
        #chi1 = md.compute_chi1(t[::100])
        metric = np.concatenate(( phi[1] , psi[1])  , axis = 1 )

        #for i in range(len(metric)): metric[i] = np.reshape(metric[i] , (5, 2), order='F')

        #metric2 = np.ndarray( (np.shape(phi[1])[0],  np.shape(phi[1])[1], 2) )
        cutoff = 5

        nb_cluster=0
        while  metric.shape[0] !=0 :

            nb_cluster+=0
            #medoid = find_medoid(metric)
            print(len(metric))
            medoid = metric[random.randint(0,len(metric)-1)]
            #take rqandom frame
            #rand_frame=
            #medoid =


            #for structure in metric.shape[0] :
            diff = metric - medoid
            #dist = np.sum(metric[])
        #for i in range(len(metric2)): metric2[i] = np.reshape(metric[i] , (5, 2), order='F')
            dist = np.sqrt(np.sum(np.power(diff,2) , axis= 1 ))
            list= np.array([])
            for i in range (len(dist)):
                #print(dist[i])
                if dist[i] < cutoff:
                    #metric = np.delete(metric, i, 0)
                    list = np.append(list , [i] )

            if nb_cluster == 0 :
                #clusters = np.array(list , dtype=object)
                clusters = []
                clusters.append(list)
                nb_cluster += 1
            else :
                #clusters = np.append( clusters , list , axis = 1 )
                clusters.append(list)

                nb_cluster += 1
            for j in list[::-1]  :

                metric = np.delete(metric, j, 0)

        return  nb_cluster,  clusters


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


def demuxTraj(trajPrefix, topolPrefix , nreps):
 cpptrajInputFile=open('cpptrajDemux.in' ,'w')
 cpptrajInput='parm %s \nensemble %s_0.nc trajnames ' %(topolPrefix , trajPrefix)
 for i in range(1 ,nreps)  :
  cpptrajInput+='%s_%s.nc,' %(trajPrefix ,i)
 cpptrajInput = cpptrajInput[:-1] +' nosort remlog rem.log \nautoimage\nstrip :WAT \ntrajout  demux.nc \nrun\nparm %s name top \nparmstrip :WAT parm top \nparmwrite parm top out system.wat.leap.strip.prmtop\nrun'%topolPrefix
 cpptrajInputFile.writelines(cpptrajInput)
 cpptrajInputFile.close()
 os.system('cpptraj -i cpptrajDemux.in' )
 for i in range(nreps)  :
  os.rename('demux.nc.' +str(i) , str(i)+'.demux.nc')

def runRMSF(topolPrefix, trajPrefix, nreps):
 for i in range(nreps) :
  repl= md.load(str(i)+'.demux.nc',top=topolPrefix )
  atom_indices= repl.topology.select('name CB')# [atom.index for atom.name== 'CA' in traj.topology.repl.
  #rmsf = md.rmsf(repl,   frame=0, atom_indices=atom_indices, parallel=True, precentered=False)
  #print(rmsf)

  repl.superpose(repl, atom_indices=atom_indices, ref_atom_indices=atom_indices)
  avg_xyz = np.mean(repl.xyz[:, atom_indices, :], axis=0)
  rmsf = np.sqrt(3*np.mean((repl.xyz[:, atom_indices, :] - avg_xyz)**2, axis=(0,2)))

  print(rmsf)
  '''
  print(topolPrefix , str(i) + '.demux.nc')
  universe = MDAnalysis.Universe(topolPrefix , str(i) + '.demux.nc' )
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
 plt.legend()
 plt.savefig('rmsf.png')





def main() :
 print('starting script')
 parser = argparse.ArgumentParser(description=d)

 parser.add_argument('--nreps', type=int,  help='number of replicas', required=True)
 parser.add_argument('-t', '--topol',  default='system.wat.leap.prmtop', help='topology')
 parser.add_argument('-xyz', '--trajPrefix',  default='traj', help='traj_prfix')
 parser.add_argument('--rewrite',  action='store_true',help='Force the script ro rewrite the trajectories if those have alredy been found')
 parser.add_argument('--rmsf', action='store_true',  help='compute rmsf')
 parser.add_argument('--cluster',  action='store_true',  default=False ,help='cluster trajectory')
 parser.add_argument('--rama', action='store_true', help='ramchadan')
 parser.add_argument('--exch', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--conv', action='store_true', help='provide a small analysis of the exchanges')
 parser.add_argument('--toref',action='store_true', default=False ,help='compare ramachadan to a converged simulation')

 args = parser.parse_args()
 print(args)

 nreps =  args.nreps
 topolPrefix = args.topol
 trajPrefix = args.trajPrefix
 rewrite = args.rewrite
 RMSF = args.rmsf
 toref = args.toref
 exch = args.exch
 conv = args.conv
 print(conv)
 rama = args.rama
 cluster= args.cluster
 ncluster = 10

 if os.path.isfile('./0.demux.nc') == False or args.rewrite==True :
    demuxTraj(trajPrefix,topolPrefix, nreps)
 topolPrefix = 'system.wat.leap.strip.prmtop'



 if  exch ==True : diffusionAnalysis(nreps )

 if RMSF == True:
     runRMSF( topolPrefix, trajPrefix , nreps)
 if rama == True:
    runRama( topolPrefix, trajPrefix , nreps)

 if  cluster == True:
    c, n = runCluster(  topolPrefix, trajPrefix , nreps )
    print( n , 'There was ' + str(c) +'clusters')

 if args.conv == True:
    Convergence(topolPrefix)
if __name__ == '__main__' :

    main()
