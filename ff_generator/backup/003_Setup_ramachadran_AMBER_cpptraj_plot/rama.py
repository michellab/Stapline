#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns


'''
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran

u = mda.Universe('topol.prmtop', 'mdcrd.nc')
r = u.select_atoms("resname LEU")    # selection of residues


R = Ramachandran(r).run()
'''



t ,phi ,psi,one,two ,three =np.loadtxt('dihedrals.dat' , comments='#' , unpack=True )
df = pd.DataFrame( { 'phi' : phi , 'psi' :psi } )

#sns.jointplot(x="phi", y="psi", data=df, kind="kde",  n_levels=1000 ,cmap="Spectral",norm=log_norm)

#cmap = sns.cubehelix_palette(as_cmap=True, dark=0.2, light=1)
#cmap = sns.cubehelix_palette(as_cmap=True, dark=0.2, light=1, start=1.5, rot=-.25)

#cmap = sns.dark_palette((210, 90, 60), input="husl")


plt.xlim(-180, 180)
plt.ylim(-180, 180)
sns.kdeplot(df.phi, df.psi, cmap='nipy_spectral_r',  n_levels=1000 , shade=True ) # norm=LogNorm())
#'Spectral_r'



plt.show()
