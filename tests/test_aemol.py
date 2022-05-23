from mol_translator.properties.nmr.nmr_ops import get_coupling_types
from mol_translator.structure.openbabel_converter import obmol_to_aemol, aemol_to_obmol
import numpy as np
from rdkit.Chem import AllChem
from rdkit import Chem
import openbabel.pybel as pyb
from mol_translator.aemol import Aemol
import sys
import os
from pathlib import Path
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')


'''
    Write tests for every function and property of the aemol class
'''


def test_init():
    test_mol = 'CCCC'
    mol = Aemol(test_mol)
    assert mol.info['molid'] == test_mol
    assert mol.info['filepath'] == None

    assert mol.structure['xyz'] == []
    assert mol.structure['types'] == []
    assert mol.structure['conn'] == []

    assert mol.atom_properties == {}
    assert mol.pair_properties == {}
    assert mol.mol_properties['energy'] == -404.404


def test_from_pybel():
    test_mol = 'from_pyb_test'
    test_xyz = 'tests/test_dataset/test_xyz.xyz'
    mol = Aemol(test_mol)
    pybmol = next(pyb.readfile('xyz', test_xyz))

    mol.from_ob(pybmol)
    types, xyz, conn = obmol_to_aemol(pybmol)

    assert (mol.structure['types'] == types).all() == True
    assert (mol.structure['xyz'] == xyz).all() == True
    assert (mol.structure['conn'] == conn).all() == True


def test_to_pybel():
    test_mol = 'to_pyb_test'
    mol_test = Aemol(test_mol)

    ref_mol = 'ref'
    mol_ref = Aemol(ref_mol)

    test_xyz = 'tests/test_dataset/test_xyz.xyz'
    pybmol_ref = next(pyb.readfile('xyz', test_xyz))
    mol_ref.from_ob(pybmol_ref)

    mol_ref.to_ob()
    test_pybmol = mol_ref.obmol
    mol_test.from_ob(test_pybmol)

    assert (mol_test.structure['xyz'] ==
            mol_ref.structure['xyz']).all() == True
    assert (mol_test.structure['types'] ==
            mol_ref.structure['types']).all() == True
    assert (mol_test.structure['conn'] ==
            mol_ref.structure['conn']).all() == True


def test_from_rdkit():
    test_mol = 'from_rd_test'
    mol = Aemol(test_mol)
    smiles = ' CCCC'
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol)

    mol.from_rdkit(rdmol)
    assert mol.structure['types'] is not []
    assert mol.structure['xyz'] is not []


def test_to_rdkit():
    test_mol = 'to_rd_test'
    test_xyz = 'tests/test_dataset/test_xyz.xyz'

    mol = Aemol(test_mol)
    mol.from_file_ob(test_xyz)
    assert mol.structure['types'] is not []
    assert mol.structure['xyz'] is not []
    assert mol.structure['conn'] is not []

    mol.to_rdkit()
    test_rd = mol.rdmol
    assert test_rd is not None


def test_from_file():
    test_mol = 'to_rd_test'
    test_xyz = 'tests/test_dataset/test_xyz.xyz'
    mol = Aemol(test_mol)

    mol.from_file_ob(test_xyz)
    assert mol.structure['types'] is not []
    assert mol.structure['xyz'] is not []
    assert mol.structure['conn'] is not []


def test_from_string():
    test_mol = 'CCCC'
    mol = Aemol(test_mol)
    pybmol_test = pyb.readstring('smi', test_mol)
    return mol.from_smiles(test_mol)
    assert pybmol == pybmol_test


def test_to_file_ae():
    test_mol = 'to_file_ae_test'
    test_mol = Aemol(test_mol)
    ref_xyz = 'tests/test_dataset/test_xyz.xyz'

    test_mol.from_file_ob(ref_xyz)
    tmp_file = Path("./tmp_test.xyz")
    test_mol.to_file_ae('xyz', 'tmp_test.xyz')
    assert tmp_file.is_file() == True
    os.remove(tmp_file)


def test_to_file_pyb():
    test_mol = 'to_file_pyb_test'
    test_mol = Aemol(test_mol)
    ref_xyz = 'tests/test_dataset/test_xyz.xyz'

    test_mol.from_file_ob(ref_xyz)
    pybmol = test_mol.obmol
    test_mol.from_ob(pybmol)
    tmp_file = Path("./tmp_pyb.smi")
    test_mol.to_file_ob('smi', 'tmp_pyb.smi')
    assert tmp_file.is_file() == True
    os.remove(tmp_file)


def test_check_mol_aemol():
    test_molid = 'test'
    test_file = 'tests/test_dataset/test_xyz.xyz'
    test_mol = Aemol(test_molid)
    test_mol.from_file_ob(test_file)
    pybmol = test_mol.obmol
    test_mol.from_ob(pybmol)

    test_bad_molid = 'bad_test'
    test_bad_file = 'tests/test_dataset/test_bad_mol.xyz'
    test_bad_mol = Aemol(test_bad_molid)
    test_bad_mol.from_file_ob(test_bad_file)
    test_bad_mol.from_ob(test_bad_mol.obmol)

    assert Aemol.check_mol_aemol(test_mol) == True
    assert Aemol.check_mol_aemol(test_bad_mol) == False


def test_rd_neutralise():
    test_file = ['tests/test_mols/charged_mol1.sdf',
                 'tests/test_mols/charged_mol2.sdf']

    for idx, file in enumerate(test_file):
        test_mol = Aemol(idx)
        test_mol.from_file_rdkit(file)
        test_mol.rd_neutralise(test_mol.rdmol)
        test_mol.from_rdkit(test_mol.rdmol)
        assert test_mol.check_mol_aemol() == True


def test_pyb_neutralise():
    test_file = ['tests/test_mols/charged_mol1.sdf',
                 'tests/test_mols/charged_mol2.sdf']

    for idx, file in enumerate(test_file):
        test_mol = Aemol(idx)
        test_mol.from_file_ob(file, ftype='sdf')
        test_mol.ob_neutralise(test_mol.obmol)
        test_mol.from_ob(test_mol.obmol)
        assert test_mol.check_mol_aemol() == True


def test_prop_to_file():
    test_mol = 'to_file_pyb_test'
    test_mol = Aemol(test_mol)
    ref_xyz = 'tests/test_mols/qm9_1.nmredata.sdf'
    assert os.path.isfile(ref_xyz)
    test_mol.from_file_ob(ref_xyz, 'sdf')
    test_mol.from_ob(test_mol.obmol)
    test_mol.get_bonds()
    test_mol.get_path_lengths()
    get_coupling_types(test_mol)

    atoms = len(test_mol.structure['types'])
    test_mol.pair_properties['coupling'] = np.random.randn(atoms, atoms)
    test_mol.atom_properties['shift'] = np.random.randn(atoms)

    check1 = test_mol.pair_properties['coupling']
    check2 = test_mol.atom_properties['shift']

    test_mol.prop_to_file('test.nmredata.sdf', 'nmr', 'nmredata')

    checkmol = Aemol('test')
    checkmol.from_file_ob('test.nmredata.sdf', 'sdf')
    checkmol.prop_from_file('test.nmredata.sdf', 'nmr', 'nmredata')

    os.remove('test.nmredata.sdf')
# def test_prop_fromfile():

# def test_get_all_paths():

# def test_get_bonds():

# def test_get_path_lengths:
