'''
Contains class for loading substitution matrices
'''
import os
from functools import lru_cache

import numpy as np

MATRIX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'matrices')


def matrix_path(*parts):
	'''Returns the path to a matrix file shipped with the package.'''
	return os.path.join(MATRIX_DIR, *parts)


@lru_cache(maxsize=None)
def load_paml_matrix(path):
	'''Returns a PAMLmatrix for the given path, parsing each file only once per process.'''
	return PAMLmatrix(path)


class PAMLmatrix:
	'''
	Class for constructing a matrix from a PAML dat file
	'''
	_lodd = None
	_piFreq = None

	def __init__(self,matrix_path):
		self.matrix_path = matrix_path

	@property
	def getPiFreqs(self):
		'''Loads and returns the background frequency of the selected matrix.'''
		if self._piFreq is None:
			with open ( self.matrix_path , 'r') as f:
				triangular_mx = [[num for num in line.rstrip('\n').split(' ') ] for line in f if line.strip() != "" ]
				pi_frequencies = triangular_mx.pop()
				pi_frequencies.pop()
			self._piFreq = [float(x) for x in pi_frequencies]
		return np.array(self._piFreq)

	@property
	def lodd(self):
		'''
		Calculates (only once) and returns a log odds representation of the PAML matrix provided.
		'''
		if self._lodd is None:
			with open ( self.matrix_path , 'r') as f:
				triangular_mx = [[num for num in line.rstrip('\n').split(' ') ] for line in f if line.strip() != "" ]
			pi_frequencies = triangular_mx.pop()
			pi_frequencies.pop()

			triangular_mx.insert(0,[''])
			sym_mx=np.zeros((20,20))
			for i in range(len(triangular_mx)):
				for j in range(len(triangular_mx[i])):
					if triangular_mx[i][j] != '':
						sym_mx[j][i]=triangular_mx[i][j]

			i_lower = np.tril_indices(len(sym_mx), -1)
			sym_mx[i_lower] = sym_mx.T[i_lower]

			for i in range(len(sym_mx)):
				row_sum=0
				for j in range(len(sym_mx[i])):
					row_sum= row_sum + 2*float(sym_mx[i][j])*float(pi_frequencies[i])*float(pi_frequencies[j])
				sym_mx[i][i] = float(row_sum)/float(pi_frequencies[i])**2

			self._lodd = np.log2(sym_mx)
		return self._lodd
