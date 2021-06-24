from setuptools import setup
from setuptools import find_packages
from os import path
this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, "README.md"), encoding='utf-8') as fh:
    long_description = fh.read()

with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    required = f.read().splitlines()

setup(
    name='TwinCons',
    version='0.5.2.dev0',
    description='This projects provides several packages for analysis of MSAs comprised of two sequence groups.',
    long_description = long_description,
    long_description_content_type='text/markdown',
    author='Petar Penev',
    author_email='ppenev@gatech.edu',
    url='https://github.com/petaripenev/AlignmentScore',
    packages=find_packages(exclude=('tests', 'ROC', 'data', '.vscode')),
    python_requires='>=3.5',
    install_requires=required,
    data_files=[
        ('matrices', [
            './matrices/B.dat',
            './matrices/BEHOS.dat',
            './matrices/BH.dat',
            './matrices/BO.dat',
            './matrices/BS.dat',
            './matrices/E.dat',
            './matrices/EH.dat',
            './matrices/EO.dat',
            './matrices/ES.dat',
            './matrices/H.dat',
            './matrices/LG.dat',
            './matrices/O.dat',
            './matrices/S.dat',
            './matrices/WAG.dat',
        ]),
        ('twcPKL', [
            './data/PKL/BBS_cg09_it1_lt3.pkl',
            './data/PKL/BBS_cg09_it1_lt3.pkl.json',
        ]),
    ],
    scripts = [
        './twincons/TwinCons.py',
        './twincons/twcCalculateSegments.py',
        './twincons/twcCrossValidate.py',
        './twincons/twcParseAlnDatabase.py',
        './twincons/twcSVMtest.py',
        './twincons/twcSVMtrain.py',
        './twincons/twcWrapper.py',
    ],
    classifiers = [
        "Programming Language :: Python :: 3", 
        "License :: OSI Approved :: MIT License", 
        "Operating System :: OS Independent", 
    ],
)

#python3 setup.py sdist bdist_wheel
#python3 -m twine check dist/*
#python3 -m twine upload --verbose --repository testpypi dist/*
