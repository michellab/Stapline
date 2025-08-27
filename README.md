# Stapline

## A pipeline for the parameterisation of stapled peptide residues and for the analysis of their secondary structure preferences via molecular dynamics simulations

### Requirements

We recommend creating separate conda environments for the residue fragmentation, QM torsion scans and MM energy decomposition. For details on how to create and work with conda environments, refer to: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

**1. Residue Fragmentation**

Refer to https://docs.openforcefield.org/projects/bespokefit/en/stable/getting-started/installation.html for instructions on installing bespokefit.

```
conda create -n bespokefit-env
conda activate bespokefit-env
conda install mamba
mamba install -c conda-forge openff-bespokefit
mamba install -c conda-forge ambertools
mamba install -c conda-forge nglview
```

**2. QM Scans**

Execute the commands below in the specified order:

```
conda create -n qm_scans python=3.8
conda activate qm_scans
```

[ambertools](https://ambermd.org/GetAmber.php#ambertools) 18.0 

``conda install -c omnia ambertools=18.0``

[rdkit](https://github.com/rdkit/rdkit) 2022.09.5 

``conda install -c conda-forge rdkit=2022.09.5``

``conda install mamba``

[psi4](https://github.com/psi4/psi4) 1.6.1 

``mamba install -c psi4 psi4=1.6.1``

[psiresp](https://github.com/lilyminium/psiresp) 0.4.2 

``mamba install -c conda-forge psiresp=0.4.2``

[papermill](https://github.com/nteract/papermill) 2.5.0

``pip install papermill``

After installing the above packages, pydantic will need to be downgraded to 1.10.8

``conda install -c conda-forge pydantic=1.10.8``

**3. MM Energy Decomposition**

Refer to https://sire.openbiosim.org/install.html for instructions on creating a conda environment for Sire. The version of Sire used in this work is 2023.5.1

```
conda create -n openbiosim "python<3.13"
conda activate openbiosim
conda install -n openbiosim -c conda-forge -c openbiosim sire
```

**4. MD Simulations**

Refer to https://ambermd.org/AmberMD.php for installing the latest version of AMBER.

### Contents

**Force Field Library**

Contains the force field parameters for the stapled residues in this work, along with peptide PDB templates and tleap scripts to generate MD simulation inputs.

**Parameterisation Scripts**

Contains code to prepare force field parameters for other stapled residues and/or non-proteinogenic amino acids.

**Molecular Dynamics**

Contains the SAH-p53 peptide PDB input files and scripts to prepare, run and analyse MD simulations.


### Citation

TBC