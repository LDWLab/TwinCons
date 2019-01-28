import numpy as np
"""
Contains class for loading substitution matrices
"""

class PhymlMatrix:
	"""
	Class for a matrix from dat file
	"""

	def __init__(self,matrix_path):
		self.matrix_path = matrix_path

	def calculate_Sij(self):
		
		f = open ( self.matrix_path , 'r')
		l = [[num for num in line.rstrip('\n').split(' ') ] for line in f if line.strip() != "" ]
		for x in l:
			print(x)



	

