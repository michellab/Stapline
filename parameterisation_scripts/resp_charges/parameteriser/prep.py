from rdkit.Chem import AllChem

import parameteriser.check_structure as check_structure


def build_molecule_from_smiles(smiles):
    mol, res_type, backbone_list, capping_list, sidechain_list = check_structure.check_smiles(smiles)

    AllChem.EmbedMultipleConfs(mol, numConfs=50, randomSeed=-1) # change the number of conformers here
    rmslist = []
    AllChem.AlignMolConformers(mol, RMSlist=rmslist)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=1000) # ensure all rdkit conformers have been energy minimised

    return mol, res_type, backbone_list, capping_list, sidechain_list






