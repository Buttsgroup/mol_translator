# mol_translator

Openbabel 3 variant of mol_translator. Aims to remove rdkit usage to increase molecular compatibility

Requirements:
  Openbabel v3.1.1 - Install in desired conda environment with "conda install -c conda-forge openbabel"

Convert molecules from files to python objects to other files.

Takes any molecular structure file: sdf, pdb, xyz . . . [complete list coming]
Creates a python object with the core structural information and any other information
Takes the standard python object and produces any number of structure/input files
