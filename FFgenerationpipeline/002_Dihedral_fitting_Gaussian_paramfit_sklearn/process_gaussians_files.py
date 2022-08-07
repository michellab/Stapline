#! /usr/bin/env python
d=' \
############################ \
Reads gaussians calculations files\
    output - pdb of the stationary structures \
        -name of the pdb in a txt file \
        -QM energy related to the structure \
    \n\n     Author : Marie Bluntzer     \
'



import os
import sys


def extract_stationary_structures(prefix, Opt=True):
    hfenergy=''
    atomdict={'1':'H', '6':'C' ,'7':'N' , '8':'O', '16':'S' }
    outputenergy=open('%s-energySP.dat' %(prefix) ,'w')
    outputlist=open('%s-listfiles' %(prefix),'w')
    for filestem in range(-180,180,1):
        file=prefix+str(filestem)+'.log'

        if os.path.isfile( file ):
            print ('debug 1' )
            i=0
            filepointer=open(file, 'r')
            filecontent=filepointer.readlines()
            print(file)
            for linenumber in range(len(filecontent)):
              if '\\HF='  in  filecontent[linenumber]:

                  lineenergy= (filecontent[linenumber]+filecontent[linenumber+1]).replace("\n", "").replace(" ", "").split('\\')
                  k=0
                  while lineenergy[k][:3]!= 'HF=': k+=1
                  else :
                      if lineenergy[k][-2:]!='\n' : hfenergy=lineenergy[k]
                      else : hfenergy=filecontent[linenumber+1].split('\\')[0]
                      print( hfenergy)
            if  hfenergy=='' : break
            SP=False
            for linenumber in range(len(filecontent)):
                if 'Stationary point found' in  filecontent[linenumber] :
                    SP=True
                    print('Stationary point found in file %s' %(file))
                    linecoord= linenumber
                    while 'Center     Atomic      Atomic             Coordinates (Angstroms)' not in filecontent[linecoord]:
                        linecoord+=-1

                    linecoord+=2
                    fileout=open('SP_%s_%s.pdb' %(file[:-4],i),'w')
                    #i+=1 ###Uncomment to save more than one stationanry point it to
                    fileout.write('# Energy = %s \n' %(hfenergy))
                    j=1

                    while '----------------------------------------------------' not in filecontent[linecoord+j]:
                            linesplited=filecontent[linecoord+j].split()
                            #   print (filecontent[linecoord+j].split())
                            fileout.write('ATOM     %s   %s  AAA     1    %s%s%s  1.00  0.00 \n' %(linesplited[0].rjust(2),\
                                    atomdict[linesplited[1]] , \
                                    str("{0:.3f}".format(float(linesplited[3]))).rjust(8), \
                                    str("{0:.3f}".format(float(linesplited[4]))).rjust(8), \
                                    str("{0:.3f}".format(float(linesplited[5]))).rjust(8)))
                            j+=1
                if Opt==False:
                    linecoord= linenumber
                    if 'Center     Atomic      Atomic             Coordinates (Angstroms)'  in filecontent[linecoord]:
                        linecoord+=2
                        fileout=open('SP_%s_%s.pdb' %(file[:-4],i),'w')
                        #i+=1 ###Uncomment to save more than one stationanry point it to
                        fileout.write('# Energy = %s \n' %(hfenergy))
                        j=1

                        while '----------------------------------------------------' not in filecontent[linecoord+j]:
                                linesplited=filecontent[linecoord+j].split()
                                #print (filecontent[linecoord+j].split())
                                fileout.write('ATOM     %s   %s  AAA     1    %s%s%s  1.00  0.00 \n' %(linesplited[0].rjust(2),\
                                        atomdict[linesplited[1]] , \
                                        str("{0:.3f}".format(float(linesplited[3]))).rjust(8), \
                                        str("{0:.3f}".format(float(linesplited[4]))).rjust(8), \
                                        str("{0:.3f}".format(float(linesplited[5]))).rjust(8)))
                                j+=1
                    if 'Normal termination'  in  filecontent[linenumber]:
                        SP=True
                        #print(filecontent[linenumber] , SP)

        #now write t  he energy of this structure
            if SP==True:
                print(hfenergy)
                print(hfenergy.split("HF=")[1], hfenergy.split("HF=")[0])
                outputenergy.write("%.8f\n" % float(hfenergy.split("HF=")[1]))
                outputlist.write('SP_%s_%s.pdb \n' %(file[:-4],i))


if __name__ == "__main__" :
    if len(sys.argv) != 2 : sys.exit(' Usage : ExtractStationaryStructures.py prefix_  . With your files being prefix_angle.log ')
    else :extract_stationary_structures(  sys.argv[1])
