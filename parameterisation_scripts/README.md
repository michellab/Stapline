**resp_charges**

Contains:

-Jupyter notebook ``resp_charges.ipynb`` for RESP charge calculations. The required input is a SMILES string of the residue to be parameterised, with acetyl- and N-methyl capping groups. There are constraints on capping group atoms to have a charge sum of 0, so when they are cleaved, the resulting residue has an integer charge. The calculated partial charges will be output in a file called ``residue_resp_charges.csv``

-helper scripts inside the parameteriser/ directory

The number of conformers used can be modified within the ``parameteriser/runner.py`` script. If parameterising a residue that needs to have a total charge other than 0, the following code can be added to the ``runner.py`` script in the ``get_resp_job`` function (example for a residue with +1 charge):
``constraints.add_charge_sum_constraint_for_molecule(psirespmol, charge=+1, indices=backbone_list+sidechain_list)``

The ``resp_charges.ipynb`` notebook can also be used to calculate RESP charges for other non-proteinogenic, unstapled residues, provided they have the usual protein backbone C(C=O)N([H]). The Calpha atom can be mono- or di-substituted.

As the RESP charge calculations can take a few hours to complete, they can be run in the background with:
``nohup jupyter nbconvert --to notebook --execute resp_charges.ipynb --output output.ipynb > output.log 2>&1 &``

**fragment_resp_charges**

Contains:

-Jupyter notebook ``resp_charges_fragments.ipynb`` for RESP charge calculations of fragments. It performs the same tasks as the full residue RESP calculation script, with the exception that there are no constraints on capping group atoms to have a charge sum of 0. The number of conformers used can be directly modified within the notebook. The calculated partial charges will be output in a file called ``fragment_resp_charges.csv``

Calculations can be run in the background with:
``nohup jupyter nbconvert --to notebook --execute resp_charges_fragments.ipynb --output output.ipynb > output.log 2>&1 &``

**torsion_fitting**

Contains:

-molecule_fragmenter.ipynb (and helper scripts utils.py and utils_toolkit.py)

(run with the bespokefit conda environment) Takes as input the SMILES string of the capped stapled residue and produces a list of fragments that can be used in torsion scans. We recommend manually selecting those fragments that best describe the parent stapled residue and are large enough to cover the chemical complexity of the staple but not too large so the QM scans are slow and with limited throughput. When choosing multiple fragments, ideally they should have overlapping atoms.

-qm_scans_forward.ipynb and qm_scans_backward.ipynb

(run with the qm_scans conda environment) Takes as input the SMILES string of a selected fragment (obtained from the fragmentation process above). The atom indices making up the torsions to be scanned can be identified from the fragment representation within the ``qm_scans_forward.ipynb`` or ``qm_scans_backward.ipynb`` notebooks, and can be used as inputs for the ``run_papermill_forward.py`` and ``run_papermill_backward.py`` scipts (``idx_sets`` variable). The number of CPUs used can be modified from within the notebooks (``nthread`` variable). RESP charges for the fragments need to be passed in the notebooks so that mol2 files with the desired charges can be created. To do this, pass the path to the ``fragment_resp_charges.csv`` file to the ``charges_filename`` variable.

As the QM scans can take a while to complete (time depends on the size of the fragment, the number of conformers used and the number of torsions to be scanned), the notebooks can be run in the background with the help of the scripts ``run_papermill_forward.py`` and ``run_papermill_backward.py``. Within these scripts, specify the number of conformers (``conf_ids``) and the torsion list with the atom indices participating in each torsion (``idx_sets``). The scripts can simply be run with:

```
python run_papermill_forward.py
python run_papermill_backward.py
```

-plot_qm_energies_forward.ipynb and plot_qm_energies_backward.ipynb

(run with the qm_scans conda environment) Average, normalise and plot the QM energy profiles from the forward and backward scans for all conformers. The resulting csv files with the average forward and backward QM profiles are used later in the ``prepare_fitting_data.ipynb`` notebook to prepare the residual energy profile on which the torsion parameters will be fitted.

-mol2_lib_writer.ipynb

(run with the qm_scans conda environment) Takes as input the SMILES string of the capped stapled residue. Use the desired conformation of the staple (e.g. S/S, S/R etc.) to produce the correct lib files. A mol2 file of the staple is prepared with the calculated RESP charges (``capped_staple_charges.mol2``). To use the pre-calculated RESP charges of the stapled residue, pass the path to the ``residue_resp_charges.csv`` file to the ``charges_filename`` variable. Make sure you agree with the atom types assigned in the mol2 file and make any necessary changes.

The initial frcmod (force field parameter) file with GAFF2 parameters for the closed staple can be created with:

``parmchk2 -i capped_staple_charges.mol2 -f mol2 -o staple_gaff2.frcmod -a Y -p /path/to/your/conda/ambertools/environment/dat/leap/parm/gaff2.dat``

where the correct path to the AmberTools installation must be specified.

To prepare lib files, the stapled residue is split in two residues, usually asymmetric. We recommend avoiding splitting double bonds between the two residues, as attempting to create a double bond later in tleap (see the ``force_field_library`` directory) can produce either a cis- or trans- conformation (this is important for the hydrocarbon staples). The atom indices of the two atoms participating in the split bond are passed in the notebook and create the two lib files (e.g. ``SC1.lib`` and ``SC2.lib``). Note that you have to specify the head (N) and tail (C) atoms for each residue, so that the two split residues can be correctly connected to the rest of the protein backbone. To find which atoms they are, you can visualise the ``capped_staple_charges.mol2`` file in e.g. PyMOL and select the appropriate atoms for the two split residues. Again, check for any atom types wrongly assigned and correct them, and change the atom type of the Calpha atoms to CX (or CJ in the case of di-substituted Calphas). The atom names (not atom types) may also be modified. For example, the main chain carbons can be named CA, CB, CG, CD, CE, CZ etc., and their hydrogens HA, HB1, HB2, HG1, HG2 etc. Have a look at the ``template.pdb`` files to see how atoms are named in the staples we parameterised. You can ignore any errors shown when using tleap to make the lib files.

_______________________________________________________________________________________________
**The following scripts should be applied in the order they are explained for one torsion at a time**

1. mm_parameterisation.ipynb

(run with the qm_scans conda environment) This notebook prepares topology (prmtop) and coordinate (inpcrd) files for the mol2 files produced during the QM scans, in order to calculate their MM energy. Any atom types that have been wrongly assigned in the mol2 files produced during the QM scans can be modified within the ``mm_parameterisation.ipynb`` notebook. For fragments containing Calpha atoms, their atom types can also be updated to CX (or CJ in the case of di-substituted Calphas). 

For the first torsion being fitted, the frcmod file used is the one made in the ``mol2_lib_writer.ipynb`` notebook (``staple_gaff2.frcmod``). Before using it, copy any bonds, angles and torsions involving the Calpha atom in new lines and replace the c3 atom type with the CX (or CJ) atom type, and also include its mass and Lennard-Jones parameters. Additionally, copy the phi, psi, phi' and psi' torsions from the AMBER force field. For example, for the bond between the Calpha and the backbone N, copy the following line:

``c3-ns  263.80   1.462``

to:

``CX-ns  263.80   1.462``

For more examples, see the frcmod files for the pentenyl glycine, cysteine and Click staple (lines involving the CX atom). 

Note that in this step, tleap may show some errors about missing parameters (e.g. angles, torsions) - this is usually because of the addition of extra hydrogen atoms on the fragments that do not exist in the original unfragmented stapled residue. The missing parameters can be located in the ``/path/to/your/conda/ambertools/environment/dat/leap/parm/gaff2.dat`` file. Alternatively, you can use the ``parmchk`` command with a mol2 file of the fragment as input, and get the missing parameters from the resulting frcmod file. An example command is:

``parmchk2 -i fragment.mol2 -f mol2 -o fragment.frcmod -a Y -p /path/to/your/conda/ambertools/environment/dat/leap/parm/gaff2.dat``

For example, for the ``CCCCn1nncc1`` fragment from the Click staple, the following parameters are reported as missing from the frcmod file produced for the full stapled residue (``staple_name_gaff2.frcmod``):
```
/home/eva/anaconda3/envs/qm_resp_3.8/bin/teLeap: Error!
Could not find angle parameter: nc - cc - h4
Building proper torsion parameters.

/home/eva/anaconda3/envs/qm_resp_3.8/bin/teLeap: Error!
 ** No torsion terms for  na-cc-cc-h4

/home/eva/anaconda3/envs/qm_resp_3.8/bin/teLeap: Error!
 ** No torsion terms for  nc-nc-cc-h4

/home/eva/anaconda3/envs/qm_resp_3.8/bin/teLeap: Error!
 ** No torsion terms for  h4-cc-cc-h4
```

This is due to the addition of an extra h4 hydrogen in this fragment, which does not exist in the closed stapled residue. The missing parameters should be copied from the ``gaff2.dat`` file:

```
h4-cc-nc   61.700     121.140
h4-cc-cc-na   4   16.000       180.000           2.000
h4-cc-nc-nc   2    9.500       180.000           2.000
h4-cc-cc-h4   4   16.000       180.000           2.000
```

In this repository, you can see all necessary modifications (parameters involving the CX atom and missing torsions for the fragments) in the ``click_staple_gaff2.frcmod``file (compare it with the ``staple_name_gaff2.frcmod``file made with ``parmchk``). We then use the modified file to obtain the prmtop and inpcrd files for the MM energy decomposition step.

2. mm_energy_decomposition.ipynb

(run with the openbiosim conda environment) Computes the MM energy of the QM conformers (prmtop and inpcrd files prepared in the previous step). In particular, the MM energy of the torsion being fitted is subtracted from the total MM energy of the fragment. The number of CPUs used can be modified from within the notebook (``n_cpus)`` variable). The torsion potential expression from the GAFF2 force field can be accessed here, along with the torsion potential expression of any other torsion. This is also useful for double-checking the potential expressions of previous torsions that have been already fitted. 

3. prepare_fitting_data.ipynb

(run with the qm_scans conda environment) Computes the residual energy (QM - (MM,total - MM,torsion)) profile, which is used for torsion fitting.

4. torsion_fitting.ipynb

(run with the any conda environment that has ``scikit-learn`` installed) Performs torsion fitting on the residual energy (QM - (MM,total - MM,torsion)) profile. The max. number of periods used during fitting can be specified with the ``n_periods`` variable, we recommend using max. 4 periods as is done in the AMBER force fields. For each periodicity and phase combination (e.g. 1 period and phase 0, 1 period and phase 180, 2 periods and phases 0,0, 2 periods and phases 0,180, ..., 4 periods and phases 180,180,180,180), the code identifies the barrier heights k that best fit the residual potential energy. 

The quality of each fit is quantified with the mean squared error (MSE) between the fitted potential and the residual potential. Choosing the fit with the lowest MSE is usually a good option, however we strongly recommend to check how the fitted potential describes the residual potential, and particularly whether it describes the minima from the QM profile (not the residual profile) well. If necessary, a different fit can be selected. 

It is possible that some of the k constants will be negative. To be compatible with the AMBER force fields, they are converted to positive values, accompanied by a phase flip. E.g. if the original k constant is -0.19 and has a phase of 0, it will be changed to 0.19 with a phase of 180.

Torsions generally have the following form in the frcmod files:

```
c3-c3-c3-c3   1    0.130         0.000          -3.000
c3-c3-c3-c3   1    0.290       180.000          -2.000
c3-c3-c3-c3   1    0.110         0.000           1.000
```

where ``c3-c3-c3-c3`` are the four atom types participating in the torsion. This torsion has 3 periods, with the formatting having the 1st period last (``1.000``), and the other periods placed above it with a ``-`` symbol. The phase for the 1st period is ``0.000`` and the force constant is ``0.110``, the phase for the 2nd period is ``180.000`` and the force constant is ``0.290`` etc. 

The potential expression for the above torsion would be ``0.11 cos(phi) + 0.13 cos(3 phi) + 0.29 cos(2 phi - 3.14159) + 0.53``

In this example, ``1`` is the multiplicity of the torsion. 

Other torsions have higher multiplicities, for example the torsion below has a multiplicity of 9:

```
hc-c3-c3-ns   9    1.400         0.000           3.000
```

The potential expression for this torsion would be ``0.156 cos(3 phi) + 0.156`` (the force constant gets divided by the multiplicity). Therefore, for torsions with multiplicities > 1, the force constants obtained from the fitting procedure need to be multiplied with the multiplicity number before being updated in the frcmod file. The multiplicity of the torsion is specified with the ``multiplicity`` variable.

Assuming that the parameters you select from the fitting procedure for the ``c3-c3-c3-c3`` torsion are:

```
fitted k [0.12034869649420327, 0.43640321713044394, 0.5014242310974617, 0.19221578171410886]
fitted periods [1, 2, 3, 4]
fitted phases [0, 0, 0, 0]
```

This torsion has a multiplicity of 1, so the force constants get multiplied with 1. The torsion entry in the frcmod file should be substituted with the following:

```
c3-c3-c3-c3   1    0.192         0.000          -4.000
c3-c3-c3-c3   1    0.501         0.000          -3.000
c3-c3-c3-c3   1    0.436         0.000          -2.000
c3-c3-c3-c3   1    0.120         0.000           1.000
```
______________________________________________________________________________________________________________________

After a torsion has been parameterised and the frcmod file updated, the above steps can be followed again to parameterise a new torsion. Note that in the ``mm_parameterisation.ipynb`` notebook, the frcmod file used should be the **updated** one with the new parameters for the first torsion.

Once all torsions have been fitted and their parameters updated in the frcmod file, there will be missing parameters connecting the stapled residue backbone to the backbones of other residues. These parameters are copied from the AMBER force field and modified to contain the correct mix of uppercase and lowercase atom types. These torsions are common among all staples and can be copied from the force fields we have already prepared. See the pentenyl glycine, cysteine xylene and Click staple frcmod files. These are the parameters needed:

```
BOND
c -N   490.0    1.335       JCC,7,(1986),230; AA
C -ns  490.0    1.335       JCC,7,(1986),230; AA

ANGLE
ns-C -O      80.0      122.90   AA general
N -c -o      80.0      122.90   AA general
C-ns-CX      50.0      121.90   AA general  (was C-N-CT)
c-N-XC       50.0      121.90   AA general  (was C-N-CT)
c-N-CX       50.0      121.90   AA general  (was C-N-CT)
C -ns-hn     50.0      120.00   AA general, gln, asn,changed based on NMA nmodes
c -N -H      50.0      120.00   AA general, gln, asn,changed based on NMA nmodes
CX- C-ns     70.0      116.60   AA general  (was CT-C-N)
XC- C-ns     70.0      116.60   AA general  (was CT-C-N)
CX- c -N     70.0      116.60   AA general  (was CT-C-N)

DIHE
o -c -ns-CX   4   10.00          180.0           2.         AA,NMA, X -C -N -X
CX-c -ns-CX   4   10.00          180.0           2.         AA,NMA, X -C -N -X
CX-c -ns -hn  4   10.00          180.0           2.         AA,NMA, X -C -N -X
O -C -ns-CX   4   10.00          180.0           2.         AA,NMA, X -C -N -X
o -c -N-XC    4   10.00          180.0           2.         AA,NMA, X -C -N -X
O -C-ns-hn    1    2.50          180.0          -2.         JCC,7,(1986),230
O -C-ns-hn    1    2.00            0.0           1.         J.C.cistrans-NMA DE
o -c- N- H    1    2.50          180.0          -2.         JCC,7,(1986),230
o -c- N- H    1    2.00            0.0           1.         J.C.cistrans-NMA DE
C-ns-CX-c     1    0.00          0.0            -4.         
C-ns-CX-c     1    0.42          0.0            -3.         
C-ns-CX-c     1    0.27          0.0            -2.
C-ns-CX-c     1    0.00          0.0             1.         four amplitudes and phases for phi, C -N -CX-C from parm10.dat
C-ns-CX-c3    1    0.00          0.0            -4.
C-ns-CX-c3    1    0.80          0.0            -3.
C-ns-CX-c3    1    1.80          0.0            -2.        
C-ns-CX-c3    1    2.00          0.0             1.         four amplitudes and phases for phi', CT-CX-N -C from frcmod.ff14SB
ns-CX-c-N     1    0.00          0.0            -4.           
ns-CX-c-N     1    0.55        180.0            -3.
ns-CX-c-N     1    1.58        180.0            -2.
ns-CX-c-N     1    0.45        180.0             1.         four amplitudes and phases for psi, N -CX-C -N from parm10.dat
c3-CX-c-N     1    0.000         0.0            -4.
c3-CX-c-N     1    0.400         0.0            -3.
c3-CX-c-N     1    0.200         0.0            -2.
c3-CX-c-N     1    0.200         0.0             1.         four amplitudes and phases for psi', N -C -CX-CT from frcmod.ff14SB
h1-CX-ns-C    6    0.00          0.0             2.         JCC,7,(1986),230 (was X -CT-N -X )
N -c -CX-h1   6    0.00          0.0             2.         JCC,7,(1986),230 (was X-C-CT-X)
XC-C -ns-CX   4   10.00          180.0           2.         AA,NMA, X -C -N -X
XC-C -ns-hn   4   10.00          180.0           2.         AA,NMA, X -C -N -X
CX-c -N -H    4   10.00          180.0           2.         AA,NMA, X -C -N -X
CX-c -N-XC    4   10.00          180.0           2.         AA,NMA, X -C -N -X
```