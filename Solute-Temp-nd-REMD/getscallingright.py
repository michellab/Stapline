import pickle
import numpy as np
import shutil
import os
import sys

def readremlog():
	exchanges=[]
	f= open('rem.log', 'r')
	lines= f.readlines()
	i=1
	line=lines[-i]
	while lines[-i][0] != '#':

		i+=1


		exchanges.append(float(line.split()[8]))
		line=lines[-i]
	return exchanges

exchanges=readremlog()
nb_replicas=len(exchanges)
if exchanges.pop(0)!=0.00 : sys.exit('Something weird here!')
exchange=np.array(exchanges)

if np.average(exchange)>0.8: print ('not a bad start!')

else:

	if os.path.isfile('previous_sims.pkl'):
		pickle_in = open("previous_sims.pkl","rb")
		scales = pickle.load(pickle_in)
		scale=scales[-1]
		nb_previous=len(scales)+1
	else :
		scales=[]

		scale=np.linspace(1, 0.6, num=nb_replicas)
		nb_previous=1

##### clean!

path = os.getcwd()
print(path+'/sim_%s'%(nb_previous))
os.mkdir(path+'/sim_%s'%(nb_previous) )
for  files in os.walk(path):
	for filename in files:

		if  filename[:3] in ['tra','rem', 'log' ]  :
			print(path+'/'+ filename, path+'/sim_%s/'%(nb_previous)+filename)
			shutil.move(path+'/'+ filename, path+'/sim_%s/'%(nb_previous)+filename)

####

new_scale=np.linspace(1,scale[-1]+0.05, num=nb_replicas)
scales.append(new_scale)
print(new_scale)
printed_scale=''
for element in new_scale:
	printed_scale+=element
	printed_scale+=' '
os.system('python ~/GIT_REST1/Solute-temp_REST2.py -nreps 8 -scale %s' %(printed_scale)))


pickle_out = open("previous_sims.pkl","wb")
pickle.dump(scales, pickle_out)
pickle_out.close()
