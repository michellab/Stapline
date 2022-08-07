import sys
import numpy as np 

ifile = sys.argv[1]
reader = np.loadtxt(ifile,usecols=[0])

energies = []

for val in reader:
    energies.append(val*627.509)

#find the relative
minimum = min(energies)

new_energies = []

for val in energies:
    new_val = val - minimum
    new_energies.append(new_val)

print(new_energies)
outputenergy=open('energyconvert.dat'  ,'w')
for energy in new_energies :     outputenergy.write(str(energy) + '\n')
