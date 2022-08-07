# Custom_residus_parametrization-

Workflow:

 001 Use Knime and Rdkit to produce inputs for R.E.D
	-inputs can be given as sdf file or drawn directly from Marvin interface
	-outputs R.E.D readable files

	
 !!! After 001 some manual steps are requiered : 
	-elimination of unwanted fragments 
	-rename atoms : CA , CB ... 
	-partial charges check and reweighting 

 002 Production of the .lib file containing charges and atoms types 
     Python script to produce gaussians inputs files offering a sampling over each dihedral
     

 003 Paramfit Module 

 004 small Simulations production for ramachadran plotting


 !!! Advisabled manual Step : change atomtypes and RESIDUS name in .lib and .frcmod if overlapping with AMBER atomtypes 

