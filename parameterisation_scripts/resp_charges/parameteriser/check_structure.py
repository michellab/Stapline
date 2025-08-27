from rdkit import Chem

def get_match_bond_indices(mol, match_atom_indices):
    bond_indices = []
    frag_types = []
    for query_bond in mol.GetBonds():
        if (
            query_bond.GetBeginAtomIdx() in match_atom_indices # beginning atom is backbone but ending atom is not
            and query_bond.GetEndAtomIdx() not in match_atom_indices
        ):
            bond_indices.append(query_bond.GetIdx())
            begin_at = query_bond.GetBeginAtomIdx()
            end_at = query_bond.GetEndAtomIdx()
            frag_types.append([end_at, get_type(mol, begin_at)])

        elif (
            query_bond.GetBeginAtomIdx() not in match_atom_indices # beginning atom is not backbone but ending atom is
            and query_bond.GetEndAtomIdx() in match_atom_indices
        ):
            bond_indices.append(query_bond.GetIdx())
            begin_at = query_bond.GetBeginAtomIdx()
            end_at = query_bond.GetEndAtomIdx()
            frag_types.append([begin_at, get_type(mol, end_at)])
    return bond_indices, frag_types


def get_type(mol, at):
    if (
        mol.GetAtomWithIdx(at).GetSymbol() == "C"
        and str(mol.GetAtomWithIdx(at).GetHybridization()) == "SP2"
    ):
        return "capping"
    elif mol.GetAtomWithIdx(at).GetSymbol() == "N":
        return "capping"
    elif (
        mol.GetAtomWithIdx(at).GetSymbol() == "C"
        and str(mol.GetAtomWithIdx(at).GetHybridization()) == "SP3"
    ):
        return "sidechain"
    else:
        raise Exception("Problem parsing the structure")


def check_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    patt_backbone = Chem.MolFromSmarts("C(C=O)N([H])")
    mol = Chem.AddHs(mol)

    if mol.HasSubstructMatch(patt_backbone) is False:
        raise Exception("Could not recognise a backbone pattern in your molecule")

    if (len(mol.GetSubstructMatches(patt_backbone))) == 1:
        print("Treating your molecule as an unstapled residue")

        atom_indices = mol.GetSubstructMatch(patt_backbone)
        bond_indices, frag_types = get_match_bond_indices(mol, atom_indices)
        mol_fragmented = Chem.FragmentOnBonds(mol, bond_indices, addDummies=False)
        fragments = Chem.GetMolFrags(mol_fragmented, sanitizeFrags=False)
        #frag_types.insert(fragments_mols.index(tuple(sorted(atom_indices))), "backbone")
        backbone, capping_groups, sidechain = [], [], []

        # now need to separate the fragments into backbone, sidechain and capping
        for i in atom_indices:
            backbone.append(i)
        
        for frag in fragments:
            for frag_type in frag_types:
                if frag_type[0] in frag:
                    if frag_type[1] == 'capping':
                        for i in frag:
                            capping_groups.append(i)
                    elif frag_type[1] == 'sidechain':
                        for i in frag:
                            sidechain.append(i)

        return (mol, 1, backbone, capping_groups, sidechain)
        
    elif (len(mol.GetSubstructMatches(patt_backbone))) == 2:
        print("Treating your residue as a stapled residue")

        atom_indices_1 = mol.GetSubstructMatches(patt_backbone)[0]
        atom_indices_2 = mol.GetSubstructMatches(patt_backbone)[1]

        bond_indices, frag_types = get_match_bond_indices(
            mol, atom_indices_1 + atom_indices_2
        )
        mol_fragmented = Chem.FragmentOnBonds(mol, bond_indices, addDummies=False)
        fragments = Chem.GetMolFrags(mol_fragmented, sanitizeFrags=False)

        backbone, capping_groups, sidechain = [], [], []

        # now need to separate the fragments into backbone, sidechain and capping
        for i in atom_indices_1+atom_indices_2:
            backbone.append(i)

        for frag in fragments:
            for frag_type in frag_types:
                if frag_type[0] in frag:
                    if frag_type[1] == 'capping':
                        for i in frag:
                            capping_groups.append(i)
                    elif frag_type[1] == 'sidechain':
                        for i in frag:
                            sidechain.append(i)

        # checking for duplicates (stapled sidechain will have duplicates as it is connected to the backbone twice)
        backbone = list(set(backbone))
        sidechain = list(set(sidechain))
        capping_groups = list(set(capping_groups))

        return (mol, 2, backbone, capping_groups, sidechain)

    else:
        raise ExceptionError(
            "3 or more backbone fragments found in your molecule."
            "This is not currently supported."
        )
        return 0

