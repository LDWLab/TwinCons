import setuptools
from setuptools import find_packages
from os import path
this_directory = path.abspath(path.dirname(__file__))

with open("README.md", "r") as fh:
    long_description = fh.read()

with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    required = f.read().splitlines()

setuptools.setup(
    name='TwinCons',
    version='0.4.6dev',
    description='This projects provides several packages for analysis of MSAs comprised of two sequence groups.',
    long_description = long_description,
    long_description_content_type='text/markdown',
    author='Petar Penev',
    author_email='ppenev@gatech.edu',
    url='https://github.com/petaripenev/AlignmentScore',
    packages=find_packages(exclude=('tests', 'ROC', 'data', '.vscode')),
    python_requires='>=3.5',
    install_requires=required,
    data_files=[('matrices', [
        './twincons/matrices/B.dat',
        './twincons/matrices/BEHOS.dat',
        './twincons/matrices/BH.dat',
        './twincons/matrices/BO.dat',
        './twincons/matrices/BS.dat',
        './twincons/matrices/E.dat',
        './twincons/matrices/EH.dat',
        './twincons/matrices/EO.dat',
        './twincons/matrices/ES.dat',
        './twincons/matrices/H.dat',
        './twincons/matrices/LG.dat',
        './twincons/matrices/O.dat',
        './twincons/matrices/S.dat',
        './twincons/matrices/WAG.dat'
        ])
    ],
    scripts=['./twincons/TwinCons.py'],
        classifiers =( 
            "Programming Language :: Python :: 3", 
            "License :: OSI Approved :: MIT License", 
            "Operating System :: OS Independent", 
        ),
)

#python3 setup.py sdist bdist_wheel
#python3 -m twine check dist/*
#python3 -m twine upload --verbose --repository testpypi dist/*
