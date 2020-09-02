from distutils.core import setup
from setuptools import find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

setup(
	name='TwinCons',
	version='0.2dev',
	description='',
	long_description=open('README.txt').read(),
	author='Petar Penev',
	author_email='ppenev@gatech.edu',
	url='https://github.com/petaripenev/AlignmentScore',
	license=license,
	packages=find_packages(exclude=('tests', 'docs')),
)