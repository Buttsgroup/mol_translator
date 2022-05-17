from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

install_requires = ["openbabel",
                    "rdkit-pypi >= 2022.3.2.1",
                    "numpy",
                    "pandas >= 1.4.2",
                    "matplotlib >=3.5.2",
                    "scikit-learn >= 1.0.2",
                    "tqdm >= 4.64.0"]

setup(
    name='mol_translator',
    version='0.1.8',
    description='Molecule processing package for working with NMR data. Used for pipeline with DFT calculation & machine learning preprocessing',
    long_description=readme,
    author='Calvin Yiu, Will Gerrard, Butts Research Group',
    author_email='calvin.yiu@bristol.ac.uk',
    url='https://github.com/Buttsgroup/mol_translator',
    license=license,
    packages=find_packages(exclude=('tests',)),
    install_requires=install_requires,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires="==3.9.*",
)
