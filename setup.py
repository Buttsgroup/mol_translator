from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='mol_translator',
    version='1.0.0',
    description='Molecule processing package for working with NMR data. Used for pipeline with DFT calculation & machine learning preprocessing',
    long_description=readme,
    author='Calvin Yiu, Will Gerrard, Ben Honore, Butts Research Group',
    author_email='calvin.yiu@bristol.ac.uk',
    url='https://github.com/Buttsgroup/mol_translator',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Computational Chemist',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.7',
    ],
    )