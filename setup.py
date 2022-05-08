from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

install_requires = ["openbabel >= 3.1.1",
                    "rdkit-pypi >= 2021.03.5",
                    "numpy >= 1.21.2",
                    "pandas >= 1.3.4",
                    "matplotlib",
                    "scikit-learn >= 1.0.2",
                    "seaborn >= 0.11.2",
                    "tqdm >= 4.62.3"]

setup(
    name='mol_translator',
    version='0.1.0',
    description='Molecule processing package for working with NMR data. Used for pipeline with DFT calculation & machine learning preprocessing',
    long_description=readme,
    author='Calvin Yiu, Will Gerrard, Butts Research Group',
    author_email='calvin.yiu@bristol.ac.uk',
    url='https://github.com/Buttsgroup/mol_translator',
    license=license,
    packages=find_packages(where="mol_translator/"),
    package_dir={"": "mol_translator"},
    install_requires=install_requires,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Computational Chemist',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.9',
    ],
)
