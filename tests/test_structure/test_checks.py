from mol_translator.aemol import aemol
from mol_translator.structure.checks import run_all_checks, check_valence, check_missing_H, check_bonds_are_plausible_and_atom_overlap

def test_run_all_checks():
    test_molid = 'test'
    test_file = 'tests/test_dataset/test_xyz.xyz'
    test_mol = aemol(test_molid)
    test_mol.from_file_pyb(test_file)
    test_mol.from_pybel(test_mol.pybmol)

    test_bad_molid = 'bad_test'
    test_bad_file = 'tests/test_dataset/test_bad_mol.xyz'
    test_bad_mol = aemol(test_bad_molid)
    test_bad_mol.from_file_pyb(test_bad_file)
    test_bad_mol.from_pybel(test_bad_mol.pybmol)

    assert run_all_checks(test_mol) == True
    assert run_all_checks(test_bad_mol) == False
