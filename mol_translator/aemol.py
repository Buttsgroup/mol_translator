# Copyright 2020 Will Gerrard, Calvin Yiu
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

import pybel as pyb

from mol_translator.structure.pybel_converter import pybmol_to_aemol, aemol_to_pybmol
from mol_translator.structure.rdkit_converter import rdmol_to_aemol, aemol_to_rdmol
from mol_translator.structure import structure_write as strucwrt

import mol_translator.properties.property_io as prop_io

from mol_translator.structure import find_paths as pathfind

from mol_translator.properties.nmr.nmr_write import write_nmredata

from mol_translator.structure.checks import run_all_checks

from rdkit.Chem import AllChem
from rdkit import Chem

from mol_translator.properties.charge.charge_ops import rdkit_neutralise, pybel_neutralise

class aemol(object):
    """
        molecule object

        molecular information is stored in a set of dictionaries containing strings and numpy arrays

        contains functions for i/o and conversion between different python package objects

    """

    def __init__(self, molid, filepath=""):
        """
        Generates blank internal variables which stores assigned atomic properties as strings and numpy arrays.

        Args:
            molid: The identifier to the aemol obejct
            filepath: The filepath where the assigned aemol objects attributes are stored

        Returns:
            Null : Generates default internal variables
        """

        self.info = {   'molid': molid,
                        'filepath': filepath}

        self.structure = {  'xyz': [],
                            'types': [],
                            'conn': []}

        self.atom_properties = {}

        self.pair_properties = {}

        self.mol_properties = {'energy': -404.404}


    def from_pybel(self, pybmol):
        """
        Converts molecule from pybel/openbabel object into an aemol object

        Args:
            self: aemol object
            pybmol: pybel/openbabel mol object

        Returns:
            Null: Stores interal information to aemol object
        """
        types, xyz, conn = pybmol_to_aemol(pybmol)
        self.structure['types'] = types
        self.structure['xyz'] = xyz
        self.structure['conn'] = conn

    def to_pybel(self):
        """
        Converts aemol object into pybel/openbabel object

        Args:
            self: aemol object

        Returns:
            pybmol (object): pybel/openbabel mol object
        """
        pybmol = aemol_to_pybmol(self.structure)

        return pybmol

    def from_rdkit(self, rdmol):
        """
        Converts rdkit object into aemol object

        Args:
            self: aemol object
            rdmol: rdkit mol object

        Returns:
            Null: Stores interal information to aemol object
        """
        #assumes rdmol is already 3D with Hs included
        types, xyz, conn = rdmol_to_aemol(rdmol)
        self.structure['types'] = types
        self.structure['xyz'] = xyz
        self.structure['conn'] = conn

    def to_rdkit(self, sanitize=True, removeHs=False):
        """
        Converts aemol object into rdkit objects

        Args:
            self: aemol object
            sanitize: True/False, whether the mol is cleaned via rdkit checks when converted
            removeHs: True/False, whether hydrogen atoms are explicitly stated on the rdkit object

        Returns:
            rdmol (object): rdkit object
        """
        rdmol = aemol_to_rdmol(self, sanitize, removeHs)

        return rdmol

    def from_file(self, file, ftype='xyz'):
        """
        Converts a file of filetype 'ftype' into an aemol object via pybel/openbabel

        Args:
            self: aemol object
            file: filepath to the input file containing molecular information of format ftype
            ftype: the filetype extension of the input file

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = next(pyb.readfile(ftype, file))
        self.from_pybel(pybmol)

    def batch_from_file(self, path):
        """
        Mass converts all filetypes within a folder into a rdkit object list

        Args:
            self: aemol object
            path: file path to folder containing list of molecular data

        Returns:
            suppl (list): list of rdkit mol objects
        """
        suppl = SDMolSupplier(path)
        return suppl

    def from_string(self, string, stype='smi'):
        """
        Converts a SMILES string into an aemol object via pybel/openbabel

        Args:
            self: aemol object
            string: SMILES string of molecule
            stype: The format of the string

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = pyb.readstring(stype, string)
        self.from_pybel(pybmol)

    def to_file_ae(self, format, filename):
        """
        Writes the aemol object into an 'xyz' file, currently only supports 'xyz' format

        Args:
            self: aemol object
            format: The output file extension
            filename: The name of the output file

        Returns:
            Null: saves a file of assigned format containing molecular information
        """
        if format == 'xyz':
            strucwrt.write_mol_toxyz(self.structure, filename)
        else:
            raise_formaterror(format)

    def to_file_pyb(self, format, filename):
        """
        Writes aemol object into a filetype support by pybel/openbabel

        Args:
            self: aemol object
            format: The output file extension
            filename: The name of the output file

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = self.to_pybel()
        pybmol.write(format, filename)

    def prop_tofile(self, filename, prop='nmr', format='nmredata'):
        """
        Writes NMR properties stored within aemol object into an nmredata file

        Args:
            self: aemol object
            filename: The name of the output file
            prop: The input property data type
            format: The output file extension

        Returns:
            Null: Stores interal information to aemol object
        """
        prop_io.prop_write(self, filename, prop=prop, format=format)

    def prop_fromfile(self, filename, format, prop):
        """
        Reads NMR properties into a generated aemol object

        Args:
            self: aemol object
            filename: The name of the input file
            prop: The input property data type
            format: The output file extension


        Returns:
            Null: Stores interal information to aemol object
        """
        prop_io.prop_read(self, filename, prop=prop, format=format)

    def get_all_paths(self, maxlen=5):
        """
        Uses pybel/openbabel to generate all connected paths in the aemol object

        Args:
            self: aemol object
            maxlen = defines the max distance to generate paths

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = self.to_pybel()
        self.structure['paths'] = pathfind.pybmol_find_all_paths(pybmol, maxlen)

    def get_bonds(self):
        """
        Uses pybel/openbabel to generate a connectivity matrix of bond connections within the aemol object

        Args:
            self: aemol object

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = self.to_pybel()
        self.structure['conn'] = pathfind.pybmol_get_bond_table(pybmol)

    def get_path_lengths(self, maxlen=5):
        """
        Uses pybel/openbabel to generate path lengths for the paths found

        Args:
            self: aemol object
            maxlen: defines to the max distance to generate paths

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = self.to_pybel()
        self.structure['path_len'] = pathfind.pybmol_get_path_lengths(pybmol, maxlen)

    def get_pyb_fingerprint(self, fingerprint='ecfp6'):
        """
        Uses pybel/openbabel to generate a fingerprint of the aemol object

        Args:
            self: aemol object
            fingerprint: the type of fingerprint to be generated,
            available fingerprints: ['ecfp0', 'ecfp10', 'ecfp2', 'ecfp4', 'ecfp6', 'ecfp8', 'fp2', 'fp3', 'fp4', 'maccs']

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = self.to_pybel()
        self.mol_properties[fingerprint] = pybmol.calcfp(fingerprint)

    def get_rdkit_fingerprint(self, radius=2, nBits=2048):
        """
        Uses rdkit to generate an ecfp fingerprint of the aemol object

        Args:
            self: aemol object
            radius: The distance to identify unique fragments
            nBits: The total length of the fingerprint generated

        Returns:
            Null: Stores interal information to aemol object
        """
        rdmol = self.to_rdkit()
        fp = AllChem.GetMorganFingerprintAsBitVect(rdmol,radius=radius, nBits=nBits)
        self.mol_properties['ecfp4'] = fp

    def get_rdkit_3D_mol(self):
        """
        Uses rdkit to convert the aemol object into 3D with molecular dynamics energy minimisation

        Args:
            self: aemol object

        Returns:
            Null: Stores interal information to aemol object
        """
        rdmol = self.to_rdkit()
        rdmol = AllChem.AddHs(rdmol)
        AllChem.EmbedMolecule(rdmol)
        self.from_rdkit(rdmol)

    def check_mol(self, post_check=False):
        """
        Checks the validity of aemol via basic checks:
            Valence electrons
            Atom count
            Atom distance
            Charge

        Args:
            self: aemol object
            post_check: A post DFT calculated check for missing Hs

        Returns:
            Null: Stores interal information to aemol object
        """
        return run_all_checks(self, post_check=post_check)

    def rd_neutralise(self, opt=True):
        """
        Uses rdkit to neutralise any charged atoms

        Args:
            self: aemol object
            opt: Quick conversion of molecule to 3D and apply forcefield based energy minimisation via rdkit

        Returns:
            Null: Stores interal information to aemol object
        """
        rdmol = self.to_rdkit()
        rdmol = rdkit_neutralise(rdmol)
        if opt:
            rdmol = AllChem.AddHs(rdmol)
            AllChem.EmbedMolecule(rdmol)
        self.from_rdkit(rdmol)

    def pyb_neutralise(self, opt=True):
        """
        Uses pybel/openbabel to neutralise any charged atoms

        Args:
            self: aemol object
            opt: Quick conversion of molecule to 3D and apply forcefield based energy minimisation via rdkit

        Returns:
            Null: Stores interal information to aemol object
        """
        pybmol = self.to_pybel()
        pybmol = pybel_neutralise(pybmol)
        self.from_pybel(pybmol)
        if opt:
            self.get_rdkit_3D_mol()
