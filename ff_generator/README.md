# Custom_residus_parametrization 

Workflow:

 Follow example in examples/
	-inputs are given as SMILE strings or mol2 formats.
  The user need to had capping groups to the smiles
  options can be given (basis set /functional... ) see config.json in the example folder

Returns:

 Amber compatible .lib file containing multi RESP charges, and additional frcmod file containing dihadral bond and angle parameters.  antechamber is used to provide atoms types.


Requirements:

  python==3.8
  ambertools==18.0
  rdkit==2022.09.5
  parmed==4.0.0
  psi4==1.6.1
  psiresp==0.4.2
  pytorch==1.12.1
