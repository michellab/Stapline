d='                       \
############################################### \n  \
### The PyRED server can be used online instead of this scprit\n  \
### This script mean to wrap scripts of the REDIII software \n \n \
### F.-Y. Dupradeau, A. Pigache, T. Zaffran, C. Savineau, \n \
###R. Lelong, N. Grivel, D. Lelong, W. Rosanski & P. Cieplak,\n  \
###The R.E.D. tools: Advances in RESP and ESP charge derivation \n \
### and force field library building, Phys. Chem. Chem. Phys. 2010, 12, 7821-7839 \n\n \
### of which tiny modifications has been made especially for amino acids developpement \n \
###  #NoMorePERL  \n \
###############################################  \n \
Usage : python run_resp.py confs.sdf conf.in \n \
with confs.sdf being the sdf obtain on the last step \n \
and conf.in containing a list of atoms to remove from calculations : generally ACE and NME \n \
you can access those with pymol labels : others property : atoms identifiers : index  \n \
 '

import sys
import os
import shutil

def main() :


    sdffile = sys.argv[1][:-4]

    print (sdffile)

    os.mkdir(sdffile)
    shutil.copy ( sys.argv[2] ,sdffile + '/conf.in')
    os.chdir(sdffile)
    os.system('/usr/bin/babel  -isdf ../%s.sdf -opdb %s.pdb '  %(sdffile ,sdffile ) )

    PDB_text = open(sdffile+'.pdb' , 'r')

    model_number = 1
    new_file_text = ""
    for line in  PDB_text.readlines():
         line = line.strip () #for better control of ends of lines
         if line == "ENDMDL":
             # save file with file number in name
             shutil.copy ('conf.in' , "Mol_red" + str(model_number) + ".pdb")
             output_file = open("Mol_red" + str(model_number) + ".pdb", "a")
             output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
             output_file.close()
             # reset everything for next model
             model_number += 1
             new_file_text = ""
         elif not line.startswith("MODEL"):
             new_file_text += line.replace('HETATM' , 'ATOM  ') + '\n'

    for file in range(model_number) :
        os.system('perl %s/Ante_RED-1.5.pl Mol_red%s.pdb' %( sys.path[0], file) )





    os.chdir('..')


if __name__ == '__main__':


    main()
