# Copyright 2022 Will Gerrard, Calvin Yiu
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

import openbabel.pybel as pyb
import numpy as np

from os import PathLike

# Conversion functions for rdkit, pybel
from mol_translator.structure.openbabel_converter import obmol_to_aemol, aemol_to_obmol
from mol_translator.structure.rdkit_converter import rdmol_to_aemol, aemol_to_rdmol
# Structure writer
from mol_translator.structure import structure_write as strucwrt
# Property input/output
import mol_translator.properties.property_io as prop_io
# Path finder
from mol_translator.structure import find_paths as pathfind

from mol_translator.properties.nmr.nmr_write import write_nmredata

from mol_translator.structure.checks import run_all_checks

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize as mol_std


from mol_translator.properties.charge.charge_ops import rdkit_neutralise, openbabel_neutralise

class Aemol(object):
    """
        Base class for mol_translator, contains structural and chemical information plus storage of RDKit and Openbabel molecular objects

        >>> mol = Aemol('butane')
        >>> mol.from_smiles('CCCC')
        >>> mol.structure['conn']
        array([[0, 1, 0, 0],
               [1, 0, 1, 0],
               [0, 1, 0, 1],
               [0, 0, 1, 0]], dtype=int32)
    """

    def __init__(self, molid: str, filepath: str = "") -> object:

        """
        Initilises an Aemol object with blank attributes to be written. A molid string needs to be passed to name the molecule represented internally,
        am optional filepath can be passed to store the location of the file read else it stores a blank string.

        Attributes generated currently are:
            info: dictionary of the molid and filepath
            structure: dictionary of  'xyz': xyz coordinate, 'types': atom type array, and 'conn': connectivity matrix
            rdmol: storage of rdkit representation of the molecule
            obmol: storage of openbabel representation of the molecule
            atom_properties: dictionary of atom properties
            pair_properties: dictionary of pairwise properties
            mol_properties: dictionary of molecular properties

        :param molid: The identifier to the aemol obejct
        :param filepath: The filepath where the assigned aemol objects attributes are stored

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

        self.obmol = None

        # Atom based properties: Chemical shift, Charge, etc
        self.atom_properties = {}
        # Pair properties: Coupling constant, distances, etc
        self.pair_properties = {}
        # Mol properties: Binding affinity, energy, etc
        self.mol_properties = {'energy': -404.404,
                               }

    def from_ob(self, obmol: object) -> None:
        """
        Converts Openbabel molecule object into an Aemol object

        :param obmol: Openbabel instance

        """
        types, xyz, conn = obmol_to_aemol(obmol)
        self.structure['types'] = types
        self.structure['xyz'] = xyz
        self.structure['conn'] = conn
    # Create pybel molecule object from aemol object
    def to_ob(self) -> None:
        """
        Creates an Openbabel instance of the Aemol molecule, stored internally in self.obmol

        """
        self.obmol = aemol_to_obmol(self.structure, self.info['molid'])
        return self.obmol

    def from_rdkit(self, rdmol: object) -> None:
        """
        Converts rdkit mol instance into Aemol object

        :param rdmol: RDKit mol instance

        """
        # assumes rdmol is already 3D with Hs included
        types, xyz, conn = rdmol_to_aemol(rdmol)
        self.structure['types'] = types
        self.structure['xyz'] = xyz
        self.structure['conn'] = conn

    def to_rdkit(self, sanitize: bool = True, removeHs: bool = False) -> None:
        """
        Converts aemol object into RDKit objects

        :param sanitize: bool, whether the mol is cleaned via rdkit checks when converted
        :param removeHs: bool, whether hydrogen atoms are explicitly stated on the rdkit object

        """
        self.rdmol = aemol_to_rdmol(
            self, self.info['molid'], sanitize, removeHs)
        return self.rdmol

    def from_file_ob(self, file: str, ftype: str = 'xyz', to_aemol: bool = True) -> None:
        """
        Reads a file of filetype 'ftype' into an aemol object through Openbabel

        :param file: str, filepath to the input file containing molecular information of format ftype
        :param ftype: str, the filetype extension of the input file

        """
        self.obmol = next(pyb.readfile(ftype, file))

        if to_aemol:
            self.from_ob(self.obmol)

    def from_file_rdkit(self, path: str, to_aemol: bool = True) -> None:
        """
        Reads in sdf file through RDKit, only functions properly if sdf file only contains one molecule

        :param path: str, file path to folder containing list of molecular data
        :param to_aemol: bool, whether to convert the rdmol into an aemol object

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

    def from_smiles(self, smiles: str) -> None:
        """
        Converts a SMILES string into an aemol object via Openbabel

        :param smiles: str, SMILES string of molecule

        """
        self.obmol = pyb.readstring('smi', smiles)
        self.from_ob(self.obmol)

    def to_file_ae(self, format: str, filename: str) -> None:
        """
        Writes the aemol object into an 'xyz' file through openbabel, currently only supports 'xyz' format

        :param format: str, output file extension
        :param filename: str, name of the output file

        """
        if format == 'xyz':
            strucwrt.write_mol_toxyz(self.structure, filename)
        else:
            raise_formaterror(format)

    def to_file_ob(self, format: str, filename: str) -> None:
        """
        Writes aemol object into a file with a filetype support by openbabel

        :param format: str, output file extension
        :param filename: str, output file name

        """
        self.obmol.write(format, filename)

    def to_file_rdkit(self, filename: str) -> None:
        """
        Writes rdmol object into sdf file by RDKit

        :param rdmol: object, rdkit molecule object
        :param filename: str, output file name
        """
        w = Chem.SDWriter(filename)
        w.write(self.rdmol)

    def prop_to_file(self, filename: str, prop: str = 'nmr', format: str = 'nmredata') -> None:
        """
        Writes NMR properties stored within aemol object into an nmredata file by default

        :param filename: str, output file name
        :param prop: str, property being read
        :param format: str, output file extention

        """
        prop_io.prop_write(self, filename, prop=prop, format=format)

    def prop_from_file(self, filename: str, prop: str = 'nmr', format: str = 'nmredata') -> None:
        """
        Reads NMR properties and stores internally into aemol object attributes

        :param filename: str, output file name
        :param prop: str, property being read ('nmr', 'scf', 'ic50', 'mc')
        :param format: str, input file type ('g09', 'g16', 'nmredata')

        """
        prop_io.prop_read(self, filename, prop=prop, format=format)

    def get_all_paths(self, maxlen: int = 5) -> None:
        """
        Generates all connected paths within a molecule at a given radius maxlen using openbabel

        :param maxlen: int, max distance to find paths

        """
        self.to_ob()
        self.structure['paths'] = pathfind.obmol_find_all_paths(
            self.obmol, maxlen)

    def get_bonds(self) -> None:
        """
        Generates a connectivity matrix based of the openbabel mol object

        """
        self.to_ob()
        self.structure['conn'] = pathfind.obmol_get_bond_table(self.obmol)

    def get_path_lengths(self, maxlen: int = 5) -> None:
        """
        Get the lengths of the shortest bond connections from one atom to another through openbabel

        :param maxlen: int, max distance to search for paths

        """
        self.to_ob()
        self.structure['path_len'] = pathfind.obmol_get_path_lengths(
            self.obmol, maxlen)

    def get_ob_fingerprint(self, fingerprint: str = 'ecfp4') -> None:
        """
        Generates fingerprints through openbabel mol object, stored internally under mol_properties attribute

        :param fingerprint: str, type of fingerprint to be generated ('ecfp0', 'ecfp10', 'ecfp2', 'ecfp4', 'ecfp6', 'ecfp8', 'fp2', 'fp3', 'fp4', 'maccs')

        """
        self.mol_properties[fingerprint] = self.obmol.calcfp(fingerprint)

    def get_rdkit_fingerprint(self, radius: int = 2, nBits: int = 2048, to_numpy: bool = True) -> None:
        """
        Generates fingerprints through RDKit mol object, stored internally under mol_properties attribute

        :param radius: int, radius of unique fragments
        :param nBits: int, total bit length of the fingerprint
        :param to_numpy: bool, converts fingerprint to numpy array

        """
        fp = AllChem.GetMorganFingerprintAsBitVect(
            self.rdmol, radius=radius, nBits=nBits)
        if to_numpy:
            arr = np.zeros(0,)
            DataStructs.ConvertToNumpyArray(fp, arr)
            self.mol_properties['ecfp4'] = arr
        else:
            self.mol_properties['ecfp4'] = fp

    def get_rdkit_3D(self, opt: bool = True, to_aemol: bool = True) -> None:
        """
        Converts aemol object to 3D using RDKit to add implicit hydrogens then MMFF energy minimise to find a possible conformer

        :param opt: bool, flag to run energy minimisation
        :param to_aemol: bool, flag to convert back to aemol object after function call

        """
        self.rdmol = AllChem.AddHs(self.rdmol)
        if opt:
            AllChem.EmbedMolecule(self.rdmol)
        if to_aemol:
            self.from_rdkit(self.rdmol)

    def get_ob_3D(self, opt: bool = True, to_aemol: bool = True) -> None:
        """
        Converts aemol object to 3D using Openbabel to add implicit hydrogens then MMFF energy minimise to find a possible conformer

        :param opt: bool, flag to run energy minimisation
        :param to_aemol: bool, flag to convert back to aemol object after function call

        """
        self.obmol.addh()

        if opt:
            self.obmol.localopt()
        else:
            self.obmol.make3D()
        if to_aemol:
            self.from_ob(self.obmol)

    def check_mol_aemol(self, post_check: bool = False) -> bool:
        """
        Checks the validity of aemol via basic checks:
            Valence electrons
            Atom count
            Atom distance
            Charge

        :param post_check: bool, flag to assess post DFT optimised geometries
        :return: bool, True/False for test pass
        """
        return run_all_checks(self, post_check=post_check)

    def check_mol_rdkit(self, clean: bool = False) -> bool:
        """
        Checks validity of rdmol molecule against internal rules
            Valence electrons
            Bond order
            Charge
            Fragments

        :param clean: bool, flag to clean molecule if test fails
        :return: bool, True/False for test pass
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

    def rd_neutralise(self, opt: bool = False) -> None:
        """
        Neutralises charged molecules with RDKit

        :param opt: bool, Flag to convert to 3D and energy minimise through RDKit

        """
        self.rdmol = rdkit_neutralise(self.rdmol)

        if opt:
            self.get_rdkit_3D()

    def ob_neutralise(self, opt: bool = False) -> None:
        """
        Neutralises charged molecules with Openbabel

        :param opt: bool, Flag to convert to 3D and energy minimise through Openbabel

        """
        self.obmol = openbabel_neutralise(self.obmol)
        if opt:
            self.get_ob_3D()
