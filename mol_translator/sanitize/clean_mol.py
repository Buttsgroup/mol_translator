from mol_translator.aemol import aemol
from openbabel import OBAtomAtomIter

'''
Openbabel functions to fix molecules based on predetermined rules
'''

def get_incorrect_atom_idx(pybmol):
    bad_atom_index = []
    for atom in pybmol.atoms:
        if atom.valence != atom.implicitvalence:
            # -1 from index value to account for 0 value
            bad_atom_index.append(atom.idx-1)

    return bad_atom_index

def get_neighbour_atoms(pybmol, atom_idx):
    neighbour_atoms_idx = []
    for neighbour in OBAtomAtomIter(pybmol.atoms[atom_idx].OBAtom):
        neighbour_atoms_idx.append(neighbour.GetIdx()-1)

    return neighbour_atoms_idx

def get_matched_common_bonded_atom(pybmol, bad_atom_idx):
    match_list = []
    adj_atom_list = []
    for bad_idx in bad_atom_idx:
        for neighbour in OBAtomAtomIter(pybmol.atoms[bad_idx].OBAtom):
            if neighbour.GetIdx() in match_list:
                adj_atom_list.append(pybmol.atoms[neighbour.GetIdx()-1].OBAtom)
            else:
                match_list.append(neighbour.GetIdx())

    return adj_atom_list

def get_pair_list(pybmol, bad_atom_idx):
    checked = []
    pair_list = []

    for i in bad_atom_idx:
        if i in checked: continue
        pair = []
        i_neighbours = get_neighbour_atoms(pybmol, i)

        for j in bad_atom_idx:
            if i == j: continue
            j_neighbours = get_neighbour_atoms(pybmol, j)

            if set(i_neighbours) & set(j_neighbours):
                pair.append(i)
                pair.append(j)
                pair_list.append(pair)
                checked.append(j)
                continue

    return pair_list

def fix_delocalised_valency(aemol):
    if aemol.pybmol is None:
        aemol.to_pybel()

    bad_atom_idx = get_incorrect_atom_idx(aemol.pybmol)

    if len(bad_atom_idx) == 2:
        adj_atom_list = get_matched_common_bonded_atom(aemol.pybmol, bad_atom_idx)

        for adj_atom in adj_atom_list:
            for bad_idx in bad_atom_idx:
                bond = aemol.pybmol.atoms[bad_idx].OBAtom.GetBond(adj_atom)
                if bond.GetBondOrder() == 1:
                    bond.SetBondOrder(2)
                elif bond.GetBondOrder() == 2:
                    bond.SetBondOrder(1)

        aemol.from_pybel(aemol.pybmol)

    else:
        pair_list = get_pair_list(aemol.pybmol, bad_atom_idx)

        for pair_idx in pair_list:
            fixed = False
            adj_atom_list = get_matched_common_bonded_atom(aemol.pybmol, pair_idx)

            for adj_atom in adj_atom_list:
                for bad_idx in pair_idx:
                    bond = aemol.pybmol.atoms[bad_idx].OBAtom.GetBond(adj_atom)
                    if bond.GetBondOrder() == 1:
                        bond.SetBondOrder(2)
                    elif bond.GetBondOrder() == 2:
                        bond.SetBondOrder(1)

        aemol.from_pybel(aemol.pybmol)



    print(f"Fixed delocalised valency of aemol {aemol.info['molid']}")