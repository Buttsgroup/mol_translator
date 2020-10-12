import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.aemol as aemol
import pybel as pyb
from rdkit import Chem
from rdkit.Chem import AllChem

from mol_translator.structure.pybel_converter import pybmol_to_aemol, aemol_to_pybmol
from mol_translator.structure.rdkit_converter import rdmol_to_aemol

'''
    Write tests for every function and property of the aemol class
'''
#def test_create_file(tmp_path):

def test_init():
    test_mol = 'CCCC'
    mol = aemol(test_mol)
    assert mol.info['molid'] == test_mol
    assert mol.info['filepath'] == ""

    assert mol.structure['xyz'] == []
    assert mol.structure['types'] == []
    assert mol.structure['conn'] == []

    assert mol.atom_properties == {}
    assert mol.pair_properties == {}
    assert mol.mol_properties['energy'] == -404.404

#def test_from_pybel():
    #test_mol = 'CCCC'
    #mol = aemol(test_mol)
    #mol.to_pybel()
    #mol.from_pybel()
    #assert pybmol == pybmol_test

def test_from_rdkit():
    test_mol = 'CCCC'
    mol = aemol(test_mol)
    rdmol = Chem.MolFromSmiles(test_mol)
    rdmol = Chem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol)

    return mol.from_rdkit(rdmol)
    assert mol.structure['types'] == types
    assert mol.structure['xyz'] == xyz

#def test_to_rdkit():

#def test_from_file():

def test_from_string():
    test_mol = 'CCCC'
    mol = aemol(test_mol)
    pybmol_test = pyb.readstring('smi', test_mol)
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
