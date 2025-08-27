Each directory contains:

-force field file (staple_name.frcmod) for the closed staple

-GAFF2 force field file (staple_name_gaff2.frcmod) for the closed staple as a comparison

-residue library files (lib) for two individual residues that are joined when generating inputs for MD simulations

-template stapled peptide, showing how heavy atoms are named in the stapled residues

-template tleap script for creating MD simulation inputs, where the two individual residues are joined to create a closed staple with the 'bond' command

Different lib files are provided for S- and R- conformers, and they can be mixed and matched, e.g. to create an R/S pentenyl glycine staple, use the RG1 and SG2 lib files, whereas to create an S/R pentenyl alanine staple, use the SA1 and RA2 lib files. The frcmod file is the same for either conformer for the same stapling chemistry.

The 'bond' command in the tleap scripts joins the two open stapled residues and has the following syntax:

``bond m.residue_number.atom_name m.residue_number.atom_name``

The template tleap scripts provide the atom names that have to be bonded for each staple. Note that the resulting prmtop, inpcrd and pdb files may show weird hydrogen atoms around the bonded atoms. This is ok and will be fixed in the energy minimisation step during the MD simulations.
