from os import path
from glob import glob
from setuptools import setup, find_packages
this_directory = path.abspath(path.dirname(__file__))

#https://docs.python.org/2/extending/building.html
#https://stackoverflow.com/questions/24146840/how-do-i-package-for-distribution-a-python-module-that-uses-a-shared-library
# extensions = [Extension("twincons.ext_library",
#                        ["src/library.c"],
#                        depends=["src/library.h"],
#                        include_dirs=["src"],
#               ),
# ]

with open(path.join(this_directory, "README_twc.md"), encoding='utf-8') as fh:
    long_description = fh.read()

with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    required = f.read().splitlines()

setup(
    name='TwinCons',
    version='0.7.0.dev0',
    description='TwinCons script for analysis of MSAs comprised of two sequence groups.',
    long_description = long_description,
    long_description_content_type='text/markdown',
    author='Petar Penev',
    author_email='ppenev@gatech.edu',
    url='https://github.com/petaripenev/AlignmentScore',
    packages=find_packages(include=['twincons', 'twincons.*']),
    package_data={
        'twincons': ['newton_direct_solve'],
    },
    python_requires='>=3.5',
    install_requires=required,
    #ext_modules=extensions,
    data_files=[
        ('matrices', ['./matrices/LG.dat']),
        ('matrices/BLOSUM', glob('./matrices/BLOSUM/*.out', recursive=True)),
        ('matrices/structureDerived', glob('./matrices/structureDerived/*.dat', recursive=True)),
        ('matrices/jp', glob('./matrices/jp/*.dat', recursive=True)),
    ],
    scripts = [
        './twincons/TwinCons.py'
    ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: MIT",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "twcons=twincons.entry:main",
        ],
    },
)

#python3 setup.py sdist bdist_wheel
#python3 -m twine check dist/*
#python3 -m twine upload --verbose --repository testpypi dist/*
