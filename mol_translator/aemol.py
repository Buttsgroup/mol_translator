# Copyright 2020 Will Gerrard, Calvin Yiu
# This file is part of autoenrich.

# autoenrich is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# autoenrich is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

import pybel as pyb
import numpy as np

from os import PathLike
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

# Conversion functions for rdkit, pybel
from mol_translator.structure.pybel_converter import pybmol_to_aemol, aemol_to_pybmol
from mol_translator.structure.rdkit_converter import rdmol_to_aemol, aemol_to_rdmol
# Structure writer
from mol_translator.structure import structure_write as strucwrt
# Property input/output
import mol_translator.properties.property_io as prop_io
# Path finder
from mol_translator.structure import find_paths as pathfind
# Writes NMR data to file
from mol_translator.properties.nmr.nmr_write import write_nmredata
# Run checks based off rdkit and pybel functions
from mol_translator.structure.checks import run_all_checks
# Generates descriptors
from mol_translator.descriptors.descriptor_io import get_all_descriptors

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize as mol_std


from mol_translator.properties.charge.charge_ops import rdkit_neutralise, pybel_neutralise


class aemol(object):
    """
        molecule object class
        molecular information is stored in a set of dictionaries containing strings and numpy arrays

        contains functions for i/o and conversion between different python package objects

        input output functions for popular structure formats

        functions to derive common structural info (file paths, fingerprints, etc)
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
        # object information, trying to ween off filepath usage
        self.info = {'molid': molid,
                     'filepath': filepath}
        # Structural information
        self.structure = {'xyz': [],
                          'types': [],
                          'conn': []}
        #Internal storage of rdkit/openbabel molecule objects
        self.rdmol = None

        self.pybmol = None

        # Atom based properties: Chemical shift, Charge, etc
        self.atom_properties = {}
        # Pair properties: Coupling constant, distances, etc
        self.pair_properties = {}
        # Mol properties: Binding affinity, energy, etc
        self.mol_properties = {'energy': -404.404,
                               }

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
    # Create pybel molecule object from aemol object
    def to_pybel(self):
        """
        Converts aemol object into pybel/openbabel object

        Args:
            self: aemol object

        Returns:
            pybmol (object): pybel/openbabel mol object
        """
        self.pybmol = aemol_to_pybmol(self.structure, self.info['molid'])
        return self.pybmol

    def from_rdkit(self, rdmol):
        """
        Converts rdkit object into aemol object

        Args:
            self: aemol object
            rdmol: rdkit mol object

        Returns:
            Null: Stores interal information to aemol object
        """
        # assumes rdmol is already 3D with Hs included
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
        self.rdmol = aemol_to_rdmol(
            self, self.info['molid'], sanitize, removeHs)
        return self.rdmol

    def from_file_pyb(self, file, ftype='xyz', to_aemol=True):
        """
        Converts a file of filetype 'ftype' into an aemol object via pybel/openbabel

        Args:
            self: aemol object
            file: filepath to the input file containing molecular information of format ftype
            ftype: the filetype extension of the input file

        Returns:
            Null: Stores interal information to aemol object
        """
        self.pybmol = next(pyb.readfile(ftype, file))

        if to_aemol:
            self.from_pybel(self.pybmol)

    def from_file_rdkit(self, path, to_aemol=True):
        """
        Mass converts sdf files within a folder into a rdkit object list.
        If passing a single SDF file containing multiple molecules call rdkit independently

        Args:
            self: aemol object
            path: file path to folder containing list of molecular data

        Returns:
            suppl (list): list of rdkit mol objects
        """
        suppl = Chem.SDMolSupplier(path, removeHs=False)
        for rdmol in suppl:
            if rdmol is not None:
                if rdmol.GetProp('_Name') is None:
                    rdmol.SetProp('_Name', self.info['molid'])
                self.rdmol = rdmol
            else:
                continue

        if to_aemol:
            self.from_rdkit(self.rdmol)

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
        self.pybmol = pyb.readstring(stype, string)
        self.from_pybel(self.pybmol)

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
        self.pybmol.write(format, filename)

    def to_file_rdkit(self, filename):
        """
        Writes rdmol object into sdf format by RDKit

        Args:
            self: aemol object
            rdmol: rdkit molecule object
            format: file type, sdf by default
            filename: name of the ouputfile
        """
        w = Chem.SDWriter(filename)
        w.write(self.rdmol)

    def prop_to_file(self, filename, prop='nmr', format='nmredata'):
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

    def prop_from_file(self, filename, prop='nmr', format='nmredata'):
        """
        Reads NMR properties into a generated aemol object

        Args:
            self: aemol object
            filename: The name of the input file
            prop: The input property data type eg. NMR, scf, ic50, mc
            format: The output file format extension eg. g09, g16, nmredata


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
        self.to_pybel()
        self.structure['paths'] = pathfind.pybmol_find_all_paths(
            self.pybmol, maxlen)

    def get_bonds(self):
        """
        Uses pybel/openbabel to generate a connectivity matrix of bond connections within the aemol object

        Args:
            self: aemol object

        Returns:
            Null: Stores interal information to aemol object
        """
        self.to_pybel()
        self.structure['conn'] = pathfind.pybmol_get_bond_table(self.pybmol)

    def get_path_lengths(self, maxlen=5):
        """
        Uses pybel/openbabel to generate path lengths for the paths found

        Args:
            self: aemol object
            maxlen: defines to the max distance to generate paths

        Returns:
            Null: Stores interal information to aemol object
        """
        self.to_pybel()
        self.structure['path_len'] = pathfind.pybmol_get_path_lengths(
            self.pybmol, maxlen)

    def get_pyb_fingerprint(self, fingerprint='ecfp4'):
        """
        Uses pybel/openbabel to generate a fingerprint of the aemol object

        Args:
            self: aemol object
            pybmol: pybel molecule object
            fingerprint: the type of fingerprint to be generated,
            available fingerprints: ['ecfp0', 'ecfp10', 'ecfp2', 'ecfp4', 'ecfp6', 'ecfp8', 'fp2', 'fp3', 'fp4', 'maccs']

        Returns:
            Null: Stores interal information to aemol object
        """
        self.mol_properties[fingerprint] = self.pybmol.calcfp(fingerprint)

    def get_rdkit_fingerprint(self, radius=2, nBits=2048, to_numpy=True):
        """
        Uses rdkit to generate an ecfp fingerprint of the aemol object

        Args:
            self: aemol object
            rdmol: rdkit molecule object
            radius: The distance to identify unique fragments
            nBits: The total length of the fingerprint generated

        Returns:
            Null: Stores interal information to aemol object
        """
        fp = AllChem.GetMorganFingerprintAsBitVect(
            self.rdmol, radius=radius, nBits=nBits)
        if to_numpy:
            arr = np.zeros(0,)
            DataStructs.ConvertToNumpyArray(fp, arr)
            self.mol_properties['ecfp4'] = arr
        else:
            self.mol_properties['ecfp4'] = fp

    def get_rdkit_3D(self, opt=True, to_aemol=True):
        """
        Uses rdkit to convert the aemol object into 3D with molecular dynamics energy minimisation

        Args:
            self: aemol object
            rdmol: rdkit molecule object
            opt: True/False, whether or not to run quick energy minimization

        Returns:
            Null: Stores interal information to aemol object
        """
        self.rdmol = AllChem.AddHs(self.rdmol)
        if opt:
            AllChem.EmbedMolecule(self.rdmol)
        if to_aemol:
            self.from_rdkit(self.rdmol)

    def get_pyb_3D(self, refine=True, to_aemol=True):
        """
        Use Openbabel/Pybel wrapper to convert aemol object into 3D thorough H addition and MMFF energy minimisation

        Args:
            self: aemol object
            refine: whether to run local energy optimisation

        Returns:
            Null: stores new geometry and atoms in internal aemol object
        """
        self.pybmol.addh()

        if refine:
            self.pybmol.localopt()
        else:
            self.pybmol.make3D()
        if to_aemol:
            self.from_pybel(self.pybmol)

    def check_mol_aemol(self, post_check=False):
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

    def check_mol_rdkit(self, clean=False):
        """
        Checks validity of rdmol molecule against internal rules
            Valence electrons
            Bond order
            Charge
            Fragments

        Args:
            self: aemol object
            rdmol: rdkit molecule object

        Returns:
            if clean = True: returns cleaned molecule
            else: returns null
        """
        if self.rdmol == None:
            self.to_rdkit()

        rdkit_vm = mol_std.RDKitValidation()
        molvs_vm = mol_std.MolVSValidation()

        rdkit_vm.validate(self.rdmol)
        molvs_vm.validate(self.rdmol)

        if clean:
            lfc = mol_stf.LargestFragmentChooser()
            idx = self.rdmol.GetProp('_Name')

            self.rdmol = lfc.choose(rself.rdmol)
            mol_std.Cleanup(self.rdmol)
            self.rdmol.SetProp('_Name', idx)

    def rd_neutralise(self, opt=False):
        """
        Uses rdkit to neutralise any charged atoms

        Args:
            self: aemol object
            opt: Quick conversion of molecule to 3D and apply forcefield based energy minimisation via rdkit

        Returns:
            rdmol: rdkit molecule object
        """
        self.rdmol = rdkit_neutralise(self.rdmol)
        if opt:
            return self.get_rdkit_3D(self.rdmol)

    def pyb_neutralise(self, opt=False):
        """
        Uses pybel/openbabel to neutralise any charged atoms

        Args:
            self: aemol object
            opt: Quick conversion of molecule to 3D and apply forcefield based energy minimisation via rdkit

        Returns:
            pybmol: pybel molecule object
        """
        self.pybmol = pybel_neutralise(self.pybmol)
        if opt:
            return self.get_pyb_3D(self.pybmol)

    def generate_descriptors(self):
        """
        Generate atom and pair descriptors based off openbabel and rdkit functions

        Args:
            self: aemol object

        Returns:
            Null: Updates dictionary embedded in aemol object
        """
        get_all_descriptors(self)
