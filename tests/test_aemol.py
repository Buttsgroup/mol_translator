import sys
import os
from pathlib import Path
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.aemol as aemol
import pybel as pyb
from rdkit import Chem
from rdkit.Chem import AllChem
from mol_translator.structure.pybel_converter import pybmol_to_aemol, aemol_to_pybmol
'''
    Write tests for every function and property of the aemol class
'''

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

def test_from_pybel():
    test_mol = 'from_pyb_test'
    test_xyz = './test_dataset/test_xyz.xyz'
    mol = aemol(test_mol)
    pybmol = next(pyb.readfile('xyz', test_xyz))

    mol.from_pybel(pybmol)
    types, xyz, conn = pybmol_to_aemol(pybmol)

    assert (mol.structure['types'] == types).all() == True
    assert (mol.structure['xyz'] == xyz).all() == True
    assert (mol.structure['conn'] == conn).all() == True

def test_to_pybel():
    test_mol = 'to_pyb_test'
    mol_test = aemol(test_mol)

    ref_mol = 'ref'
    mol_ref = aemol(ref_mol)

    test_xyz = './test_dataset/test_xyz.xyz'
    pybmol_ref = next(pyb.readfile('xyz', test_xyz))
    mol_ref.from_pybel(pybmol_ref)

    test_pybmol = mol_ref.to_pybel()
    mol_test.from_pybel(test_pybmol)

    assert (mol_test.structure['xyz'] == mol_ref.structure['xyz']).all() == True
    assert (mol_test.structure['types'] == mol_ref.structure['types']).all() == True
    assert (mol_test.structure['conn'] == mol_ref.structure['conn']).all() == True

def test_from_rdkit():
    test_mol = 'from_rd_test'
    mol = aemol(test_mol)
    smiles = ' CCCC'
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol)

    mol.from_rdkit(rdmol)
    assert mol.structure['types'] is not []
    assert mol.structure['xyz'] is not []

def test_to_rdkit():
    test_mol = 'to_rd_test'
    test_xyz = './test_dataset/test_xyz.xyz'

    mol = aemol(test_mol)
    mol.from_file(test_xyz)
    assert mol.structure['types'] is not []
    assert mol.structure['xyz'] is not []
    assert mol.structure['conn'] is not []

    test_rd = mol.to_rdkit()
    assert test_rd is not None

def test_from_file():
    test_mol = 'to_rd_test'
    test_xyz = './test_dataset/test_xyz.xyz'
    mol = aemol(test_mol)

    mol.from_file(test_xyz)
    assert mol.structure['types'] is not []
    assert mol.structure['xyz'] is not []
    assert mol.structure['conn'] is not []

def test_from_string():
    test_mol = 'CCCC'
    mol = aemol(test_mol)
    pybmol_test = pyb.readstring('smi', test_mol)
    return mol.from_string(test_mol)
    assert pybmol == pybmol_test

def test_to_file_ae():
    test_mol = 'to_file_ae_test'
    test_mol = aemol(test_mol)
    ref_xyz = './test_dataset/test_xyz.xyz'

    test_mol.from_file(ref_xyz)
    tmp_file = Path("./tmp_test.xyz")
    test_mol.to_file_ae('xyz', 'tmp_test.xyz')
    assert tmp_file.is_file() == True
    os.remove(tmp_file)

def test_to_file_pyb():
    test_mol = 'to_file_pyb_test'
    test_mol = aemol(test_mol)
    ref_xyz = './test_dataset/test_xyz.xyz'

    test_mol.from_file(ref_xyz)
    tmp_file = Path("./tmp_pyb.smi")
    test_mol.to_file_pyb('smi', 'tmp_pyb.smi')
    assert tmp_file.is_file() == True
    os.remove(tmp_file)

#def test_prop_tofile():

#def test_prop_fromfile():

#def test_get_all_paths():

#def test_get_bonds():

#def test_get_path_lengths:
