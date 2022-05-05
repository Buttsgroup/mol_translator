import mol_translator.structure.find_paths as fpaths
from mol_translator.aemol import Aemol
import numpy as np
import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')


def test_findallpaths():
    test_mol = Aemol(0)
    test_mol.from_file_ob('tests/test_mols/qm9_8.nmredata.sdf', ftype='sdf')
    test_mol.from_ob(test_mol.obmol)
    test_mol.get_all_paths(maxlen=5)
    # Checked by eye
    assert test_mol.structure['paths'] == [[0, 1], [0, 2],
                                           [0, 3], [0, 4], [0, 1, 5], [1, 0],
                                           [1, 0, 2], [1, 0, 3], [1, 0, 4],
                                           [1, 5], [2, 0], [
                                               2, 0, 1], [2, 0, 3],
                                           [2, 0, 4], [2, 0, 1, 5], [
                                               3, 0], [3, 0, 1],
                                           [3, 0, 2], [3, 0, 4], [
                                               3, 0, 1, 5], [4, 0],
                                           [4, 0, 1], [4, 0, 2], [
                                               4, 0, 3], [4, 0, 1, 5],
                                           [5, 1, 0], [5, 1], [5, 1, 0, 2], [5, 1, 0, 3], [5, 1, 0, 4]]


def test_find_paths():
    test_mol = Aemol(0)
    test_mol.from_file_ob('tests/test_mols/qm9_8.nmredata.sdf', ftype='sdf')

    paths = fpaths.obmol_find_paths(test_mol.obmol, 0, 1, 1)
    assert paths == [[0, 1]]
    paths = fpaths.obmol_find_paths(test_mol.obmol, 0, 5, 1)
    assert paths == []
    paths = fpaths.obmol_find_paths(test_mol.obmol, 0, 5, 2)
    assert paths == [[0, 1, 5]]
    paths = fpaths.obmol_find_paths(test_mol.obmol, 5, 4, 1)
    assert paths == []
    paths = fpaths.obmol_find_paths(test_mol.obmol, 5, 4, 3)
    assert paths == [[5, 1, 0, 4]]
    paths = fpaths.obmol_find_paths(test_mol.obmol, 5, 4, 4)
    assert paths == []


def test_getbondtable():
    test_mol = Aemol(0)
    test_mol.from_file_ob('tests/test_mols/qm9_8.nmredata.sdf', ftype='sdf')

    bond_table = fpaths.obmol_get_bond_table(test_mol.obmol)
    # Checked by eye
    assert np.array_equal(bond_table, np.asarray([[0, 1, 1, 1, 1, 0],
                                                 [1, 0, 0, 0, 0, 1],
                                                 [1, 0, 0, 0, 0, 0],
                                                 [1, 0, 0, 0, 0, 0],
                                                 [1, 0, 0, 0, 0, 0],
                                                 [0, 1, 0, 0, 0, 0]]))


def test_getpathlengths():
    test_mol = Aemol(0)
    test_mol.from_file_ob('tests/test_mols/qm9_8.nmredata.sdf', ftype='sdf')

    coupling_len = fpaths.obmol_get_path_lengths(test_mol.obmol, maxlen=5)
    # checked by eye
    assert np.array_equal(coupling_len, np.asarray([[0, 1, 1, 1, 1, 2],
                                                    [1, 0, 2, 2, 2, 1],
                                                    [1, 2, 0, 2, 2, 3],
                                                    [1, 2, 2, 0, 2, 3],
                                                    [1, 2, 2, 2, 0, 3],
                                                    [2, 1, 3, 3, 3, 0], ]))
