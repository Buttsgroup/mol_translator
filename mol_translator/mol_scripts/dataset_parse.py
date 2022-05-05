import os
import glob
import shutil
from typing import Union

from rdkit import Chem
from rdkit.Chem import AllChem

from tqdm import tqdm


def parse_sdf(folder: str, outfolder: str, id_idx: Union[int, bool] = None, id_SDTag: Union[str, bool] = None):
    files = glob.glob(f'{folder}/*.sdf')
    pbar = tqdm(files)
    for file in pbar:
        suppl = Chem.SDMolSupplier(file)
        if id_idx:
            id = file.split('/')[-1].split('.')[0].split('_')[id_idx]
        for mol in tqdm(suppl):
            if mol is None:
                continue
            if id_SDTag:
                id = mol.GetProp(id_SDTag)
            if not mol.HasProp('_Name'):
                mol.SetProp('_Name', id)
            outfile = f"{outfolder}/{id}.sdf"
            writer = Chem.SDWriter(outfile)
            writer.write(mol)
            writer.flush()


def current_dataset_filter(filepath: str, dataset_folder: list, id_idx: int, outpath: str):
    curr_mols = []
    for folder in dataset_folder:
        files = glob.glob(f'{folder}/*.sdf')
        for file in tqdm(files):
            curr_mols.append(file.split(
                '/')[-1].split('.')[0].split('_')[id_idx])

    to_filter_files = glob.glob(f'{filepath}/*.sdf')
    pbar = tqdm(to_filter_files)

    passed = 0
    filtered = 0

    for file in pbar:
        id = file.split('/')[-1].split('.')[0]
        if all(x != id for x in curr_mols):
            shutil.copy(file, f'{outpath}/{id}.sdf')
            passed += 1
        else:
            filtered += 1

        pbar.set_description(f'passed = {passed}, filtered = {filtered}')

    assert len(curr_mols) == filtered


def make_3D(folder: str, outfolder: str):
    files = glob.glob(f'{folder}/*.sdf')
    pbar = tqdm(files)
    for file in pbar:
        id = file.split('/')[-1].split('.')[0]
        suppl = Chem.SDMolSupplier(file)
        writer = Chem.SDWriter(f'{outfolder}/{id}.sdf')
        for mol in suppl:
            if mol is None:
                continue
            if not mol.HasProp('_Name'):
                mol.SetProp('_Name', id)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            writer.write(mol)
            writer.flush()


def batch_library(folder: str, outfolder: str, batch_size: int = 1000):
    files = glob.glob(f'{folder}/*.sdf')
    pbar = tqdm(files)
    batch = 0
    for idx, file in enumerate(pbar):
        if idx % batch_size == 0:
            batch += 1
        batch_folder = f'{outfolder}/batch{batch}'
        if not os.path.exists(batch_folder):
            os.makedirs(batch_folder)
        filename = file.split('/')[-1]
        outfile = f'{batch_folder}/{filename}'
        shutil.copyfile(file, outfile)
