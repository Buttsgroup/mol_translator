import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.aemol as aemol
import pybel as pybel

from mol_translator.pybel_converter.pybel_converter import pybmol_to_aemol, aemol_to_pybmol, rdmol_to_aemol

'''
    Write tests for every function and property of the aemol class
'''
test_mol = 'CCCC'
mol = aemol(test_mol)
pybmol_test = pyb.readstring('smi', test_mol)

def test_create_file(tmp_path):
    testp = tmp_path / "sub"
    testp.mkdir()
    tmp_f = testp / "tmp_file.xyz"

def test_init():
    assert mol.info['molid'] == test_mol
    assert mol.info['filepath'] == ""

    assert mol.structure['xyz'] == []
    assert mol.structure['types'] == []
    assert mol.structure['conn'] == []

    assert mol.atom_properties['shift'] == []
    assert mol.pair_properties['coupling'] == []
    assert mol.mol_properties['energy'] == -404.404

def test_from_pybel():

    mol.to_pybel()
    assert pybmol == pybmol_test

def test_from_rdkit():
    rdmol = Chem.MolFromSmiles(test_mol)
    AllChem.EmbedMolecule(rdmol)
    AllChem.MMFFOptimzeMolecule(rdmol)

    return mol.from_rdkit(rdmol)
    assert mol.structure['types'] == types
    assert mol.structure['xyz'] == xyz
    assert mol.structure['conn'] ==conn

#def test_to_rdkit():

#def test_from_file():

def test_from_string():
    return mol.from_string(test_mol)
    assert pybmol == pybmol_test

#def test_to_file_ae():

    #mol.to_file_ae(xyz, "tmp_file.xyz")
    #check tmp_file contains aemol object

#def test_to_file_pyb():
    #mol.to_file_pyb(xyz, "test_file")
    #check test_file contains pybmol object

#def test_prop_tofile():

#def test_prop_fromfile():

#def test_get_all_paths():

#def test_get_bonds():

#def test_get_path_lengths:
