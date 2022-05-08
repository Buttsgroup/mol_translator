# mol_translator

mol_translator is a python package built to process common molecular data files (.sdf, .mol2, .pdb) for computational framework including conversion to gaussian calculation files and extraction of gaussian log data, conversion to dataframes for machine learning pipeline by extracting the raw structural information. Alongside this it contains rudimentary checks to ensure molecules are feasible and Boltzmann population averaging for conformational analysis.

Applications :
  - Useful in scripts to automate dataset generation
  - Boltzmann population analysis based off NMR gaussian calculation
  - Data preprocessing for machine learning pipeline


Requirements :
  - Python == 3.9.*
  - RDKit
  - openbabel/pybel => 3.1.1
  - scikit-learn
  - numpy
  - pandas
  - tqdm