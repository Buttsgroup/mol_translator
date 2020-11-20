from rdkit import Chem

def rdmol_find_num_bonds(rdmol):
    num_bond_dict = {}
    for atom in rdmol.GetAtoms():
        num_bond_dict[atom.GetSymbol()] = [atom.GetTotalValence()]
    return num_bond_dict
