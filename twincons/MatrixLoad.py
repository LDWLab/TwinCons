import numpy as np
'''
Contains class for loading substitution matrices
'''

class PAMLmatrix:
	'''
	Class for constructing a matrix from a PAML dat file
	'''
	_lodd = None
	_dict_lodd = None
	_aa_sequence = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
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
			f = open ( self.matrix_path , 'r')
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

			#print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in np.log2(sym_mx)]))
			
			self._lodd = np.log2(sym_mx)
		return self._lodd

	@property
	def dict_lodd(self):
		'''
		Returns a dictionary form of the log odds matrix with keys tuples of aa combinations.
		'''
		if self._lodd is None:
			self._lodd = self.lodd
		if self._dict_lodd is None:
			self._dict_lodd={}
			for i in range(len(self._lodd)):
				for j in range (len(self._lodd[i])):
					self._dict_lodd[(self._aa_sequence[i],self._aa_sequence[j])] = self._lodd[i][j]
		return self._dict_lodd



