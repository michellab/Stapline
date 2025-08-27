**simulation_inputs**

Contains the pdb files for the SAH-p53 peptide series used for tleap parameterisation. 

**simulation_scripts**

Contains AMBER minimisation, equilibration and production .in files, along with the commands needed to execute them (in the sh scripts). An optional slurm job submission script that will run all min and eq steps and then run 3 prod replicates is also included. For details on the slurm job scheduling software, see https://slurm.schedmd.com/documentation.html 

The ``helicity_analysis.ipynb`` notebook analyses the area under the curve and helicity of MD trajectories with the DSSP algorithm. For speed, we use trajectory (``prod_dry.nc``) and topology (``dry.prmtop``) files where the water and ion molecules have been removed.

In order to analyse the helicity of the MD trajectories with the STRIDE algorithm, the trajectories have to be loaded in VMD (see https://www.ks.uiuc.edu/Research/vmd/) and their secondary structure can be analysed with the Timeline tool and the calculation of secondary structure option (see https://www.ks.uiuc.edu/Training/Tutorials/science/timeline/tutorial_timeline.pdf). For each simulation replicate, a ``helicity.tml`` file can be saved from the Timeline GUI and analysed with the ``vmd_analysis.ipynb`` notebook. We recommend loading the trajectories with a stride of 10 to speed up the Timeline calculation.

**md_trajectories**

Contains zipped MD simulation data for all peptides and force field/solvent combinations, as well as the Timeline helicity data obtained from VMD analysis.

To unzip the MD trajectories, run:
``tar -xvf file.tar.xz``