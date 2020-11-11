from distutils.core import setup
from setuptools import find_packages
from os import path
this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

with open(path.join(this_directory, 'LICENSE.txt'), encoding='utf-8') as f:
    license = f.read()

with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    required = f.read().splitlines()

setup(
	name='TwinCons',
	version='0.2dev',
	description='This projects provides several packages for analysis of MSAs comprised of two sequence groups.',
	long_description = readme,
	long_description_content_type='text/markdown',
	author='Petar Penev',
	author_email='ppenev@gatech.edu',
	url='https://github.com/petaripenev/AlignmentScore',
	license=license,
	packages=find_packages(exclude=('tests', 'docs')),
	python_requires='>=3.5',
	install_requires=required,
)