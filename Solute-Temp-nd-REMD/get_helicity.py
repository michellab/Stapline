#!/usr/bin/env python3

#####! $env python

d= \
"""
===============================================================================
Return a graph to plot the helicity of a MD run

------------------------
Examples:
./get_helicity --traj mdcrd.nc --topol topology.prmtop

Marie Bluntzer
s1772078@ed.ac.uk

version 0.2                                                   24/10/2019
================================================================================
"""


import sys
import argparse
import mdtraj
import numpy
import matplotlib.pyplot as plt
import os
from matplotlib import cm
import numpy as np
import seaborn as sns
import math
import matplotlib.cm as cm
import operator as o


def main():
    cmap=plt.set_cmap('Paired')
    exp=args.experimental
    experimental_data = [exp]
    overalltime= [0]
    helicity=[]
    #for n in range(len(args.traj.split(' '))):
    fig , axs = plt.subplots(len(args.traj))
    fig.set_size_inches(18.5, 20.5)
    for n in range(len(args.traj)):

        #trajfile=args.traj.split(' ')[n]
        trajfile=args.traj[n]

        try : topfile=args.topol[n]
        except :  topfile=args.topol[0]
        print(trajfile,topfile)


		#trajfile="traj%s.xtc" %(str(n))
		#topfile = 'conf%s.gro'%(str(n))
        stride=10
        traj = mdtraj.load(trajfile, top=topfile,stride=stride)

        idxs=traj.topology.select('protein')
        #traj = mdtraj.load(trajfile , atom_indices=idxs,  top=topfile )
        Hy=[]
        Hx=[]
        Hya=[]
        #dssp2 = mdtraj.compute_dssp(traj, simplified=False)
        dssp = mdtraj.compute_dssp(traj,simplified=True)
        #addavgHya=[]
        Hyaddavg=[]
        Hyrunavg=[]
        Hyrunavginbetween=[]
        Hyaddavginbetween=[]
        Hxrunavg=[]
        experimental_data.append(exp)
        overalltime.append(int(len(dssp))*2)
        #numpy.set_printoptions(threshold=sys.maxsize)
        #print(dssp)
        stapled_res= [residue.index for residue in traj.topology.residues if residue.name in  ( 'AS1' ,'AR1', 'AS2' ,'AR2' ,'AKS' , 'AKR' ,'A1R' ,'A1S' ,'A2R' ,'A2S'  )]
        print(stapled_res)

        for i in range(int(len(dssp))):
            unique, counts = numpy.unique(dssp[i : i+1], return_counts=True)
            uniqueinbetween, countsinbetween = numpy.unique(dssp[i : i+1][0][stapled_res[0]:stapled_res[1]+1], return_counts=True)
            print(dssp[i : i+1][0][stapled_res[0]:stapled_res[1]+1])
            #print(unique, counts)inbetween
            dictss = { key:val for key, val in dict(zip(unique, counts)).items() if key != 'NA' }
            dictssinbetween = { key:val for key, val in dict(zip(uniqueinbetween, countsinbetween)).items() if key != 'NA' }
            #print(dictss)
            lengthseq=0

            print(dictssinbetween)
            for val in dictss.values() :
                lengthseq+=val

            Hyrunavg.append(dictss.get('H',0)/(lengthseq)*100)
            Hyrunavginbetween.append(dictssinbetween.get('H',0)/(lengthseq)*100)
            #Hyrunavg.append(dictss.get('E',0)/(lengthseq)*100)
            Hxrunavg.append(i/stride)
            Hyaddavg.append(numpy.average(Hyrunavg[-30:]))
            Hyaddavginbetween.append(numpy.average(Hyrunavginbetween[-30:]))
            '''
            unique, counts = numpy.unique(dssp[0 : i+1], return_counts=True)
            dictss = { key:val for key, val in dict(zip(unique, counts)).items() if key != 'NA' }
            for val in dictss.values() :
                lengthseq+=val
            Hyaddavg.append(100 - dict(zip(unique, counts)).get(' ',0)/(lengthseq)*100'''
            #Hxaddavg.append(i*2*stride)



		    #	unique, counts = numpy.unique(dssp2[i*10 : (i+1)*10], return_counts=True)
		    #	Hya.append(dict(zip(unique, counts)).get('H',0)/10.0)

		    #    addavgHya.append= np.average(Hya)
            '''
            for i in range(int(len(dssp))):
            unique, counts = numpy.unique(dssp[i : (i+1)], return_counts=True)
            Hy.append(dict(zip(unique, counts)).get('H',0)/(len(dssp[0])-4)*100)
            Hx.append(i*2/10)
            print (i)
            '''
            '''
            plt.plot(Hx, Hy ,label='instant helicity' ,color='grey')
            #plt.plot(Hx, Hya ,label='alpha helicity')
            '''



            #plt.plot(Hxrunavg,addavgHy, 'k' )
        axs[n].plot(Hxrunavg,Hyaddavg , label=args.names[n] )

        axs[n].plot(Hxrunavg , Hyrunavg ,  linestyle='dashed', lw=0.1,label=args.names[n])
        axs[n].plot(Hxrunavg ,Hyaddavginbetween ,label=args.names[n]+'in between residue')
        LAST25=Hyrunavg[math.floor(len(Hyrunavg) //1.33):]
        LAST75=Hyrunavg[math.floor(len(Hyrunavg) //4):]
        LAST75=Hyrunavg[math.floor(len(Hyrunavg) //4):]
        helicity.append(numpy.average(LAST75)) #only use last 75 %
        #print(helicity)
        #helicity.append(Hyaddavg[-1])
        axs[n].set_title(args.names[n])
        axs[n].legend()
            #plt.plot(addavgHya, Hya ,label='alpha helicity')
    plt.suptitle('Variation of helicity over time')
    #axs[0].set_title('helicity averaged over 20 ps')
    #axs[1].set_title('helicity of the snapshots')
    #axs[2].set_title('helicity of theresidues inbetween the staple')


    plt.ylabel('helicity / %')
    plt.xlabel('time / ns')
#    plt.legend()
#    plt.plot(overalltime ,experimental_data)
    plt.savefig('%s' %(args.out))
    plt.close()
    newhisto(helicity, exp)
def newhisto(helicity, exp):
    names=args.names
    '''
    dpoints = np.array([['experimental', '1mfq', 9.97],
               ['experimental', '1gid', 27.31],
               ['experimental', '1y26', 5.77],
               ['standard', '1mfq', 5.55],
               ['standard', '1gid', 37.74],
               ['standard', '1y26', 5.77],
               ['lower-sca', '1mfq', 10.32],
               ['lower-sca', '1gid', 31.46],
               ['lower-sca', '1y26', 18.16]])
    '''
    dpoints=[]
    name_exp=args.name_exp

    for a in range(len(exp)):
        dpoints.append(['CD data', names[a], exp[a]])
    for a in range(len(helicity)):
        print(a,helicity)

        #print(a//len(names))
        dpoints.append([name_exp[a//len(names)], names[a%len(names)], helicity[a]])
    dpoints = np.array(dpoints)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    barplot(ax, dpoints)

def barplot(ax, dpoints):
    '''
    Create a barchart for data across different categories with
    multiple conditions for each category.

    @param ax: The plotting axes from matplotlib.
    @param dpoints: The data set as an (n, 3) numpy array
    '''

    # Aggregate the conditions and the categories according to their
    # mean values
    conditions = [(c, np.mean(dpoints[dpoints[:,0] == c][:,2].astype(float)))
                  for c in np.unique(dpoints[:,0])]
    categories = [(c, np.mean(dpoints[dpoints[:,1] == c][:,2].astype(float)))
                  for c in np.unique(dpoints[:,1])]

    # sort the conditions, categories and data so that the bars in
    # the plot will be ordered by category and condition

    #conditions = [c[0] for c in sorted(conditions, key=o.itemgetter(1))]
    #categories = [c[0] for c in sorted(categories, key=o.itemgetter(1))]
    conditions = [c[0] for c in conditions]
    categories = [c[0] for c in categories]

    dpoints = np.array(sorted(dpoints, key=lambda x: categories.index(x[1])))
    print(dpoints)
    # the space between each set of bars
    space = 0.3
    n = len(conditions)
    width = (1 - space) / (len(conditions))

    # Create a set of bars at each position
    for i,cond in enumerate(conditions):
        indeces = range(1, len(categories)+1)
        vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.float)
        pos = [j - (1 - space) / 2. + i * width for j in indeces]
        ax.bar(pos, vals, width=width, label=cond,
               color=cm.Accent(float(i) / n))

    # Set the x-axis tick labels to be equal to the categories
    ax.set_xticks(indeces)
    ax.set_xticklabels(categories)
    plt.setp(plt.xticks()[1], rotation=90)

    # Add the axis labels
    ax.set_ylabel("percentage helicity")
    ax.set_xlabel("Structure")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left')
    plt.savefig('histo.png', bbox_inches='tight')

    plt.show()


def histo(helicity, experimental_data) :
    # width of the bars
    barWidth = 0.3
    # Choose the height of the error bars (bars1)


    # Choose the height of the error bars (bars2)

    histogram=plt.figure()
    # The x position of bars
    r1 = np.arange(len(experimental_data))
    r2 = np.arange(len(helicity))
    #r2 = [ + barWidth for x in r1]
    #print(experimental_data,helicity)
    plt.bars([experimental_data,helicity],bins = int(180/15))
    #plt.hist(helicity, r2, alpha=0.5)

    plt.show()
    # Create blue bars
    plt.bar(r1, experimental_data, width = barWidth, color = 'green', edgecolor = 'black', capsize=7, label='Solute tempering data')

    # Create cyan bars
    plt.bar(r2,helicity , width = barWidth, color = 'cyan', edgecolor = 'black',  capsize=7, label='litterature')

    # general layout
    plt.xticks([r + barWidth for r in range(len(bars1))], ['cond_A', 'cond_B'])
    plt.ylabel('height')
    plt.legend()
    plt.savefig(comp.png)
    # Show graphic
    plt.show()


if __name__=="__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
    parser.add_argument('--experimental', type=float,   nargs='*' ,  help='experimental data')
    parser.add_argument('--topol',default='system.wat.leap+REST2.000.prmtop',  help='topologies files', nargs='*')
    parser.add_argument('--traj',default='traj_comp.xtc', help='trajectories files', nargs='*')
    parser.add_argument('--folders',default=False, help='folders', nargs='*')
    parser.add_argument('--out',default='helicity.png', help='trajectories files', nargs='?')
    parser.add_argument('--names',default=[], help='trajectories files', nargs='*')

    parser.add_argument('--name_exp',default=[], help='trajectories files', nargs='*')
    args = parser.parse_args()
#    colormap = plt.cm.
#    print(colormap)
#    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 1,8)])
    if args.folders :
        args.topol=[]
        args.traj=[]
        name=False
        if args.names==[]:
            name=True

        for i in args.folders:
            args.topol.append('%s/system.wat.leap+REST2.000.prmtop' %i)
            args.traj.append('%s/traj_0.nc' %i)
            if name==True : args.names.append('%s%s' %(i.split('/')[-2], i.split('/')[-1]))
    main()
