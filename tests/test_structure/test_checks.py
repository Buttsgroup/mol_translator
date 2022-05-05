from mol_translator.aemol import Aemol
from mol_translator.cleaning.checks import run_all_checks


def test_run_all_checks():
    test_molid = 'test'
    test_file = 'tests/test_dataset/test_xyz.xyz'
    test_mol = Aemol(test_molid)
    test_mol.from_file_ob(test_file)
    test_mol.from_ob(test_mol.obmol)

    test_bad_molid = 'bad_test'
    test_bad_file = 'tests/test_dataset/test_bad_mol.xyz'
    test_bad_mol = Aemol(test_bad_molid)
    test_bad_mol.from_file_ob(test_bad_file)
    test_bad_mol.from_ob(test_bad_mol.obmol)

    assert run_all_checks(test_mol) == True
    assert run_all_checks(test_bad_mol) == False
