#!/usr/bin/env python3
from Bio.SubsMat import MatrixInfo
from Bio.SubsMat import read_text_matrix
import numpy as np
import itertools, math
import pandas as pd
#pd.set_option('precision',3)
pd.set_option('display.expand_frame_repr', False)

aa_tup = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
aa_tup_h = aa_tup[:-1]
aa_tup_v = aa_tup[1:]

def matrix_reader(file_path):
	with open(file_path) as f:
		content = f.readlines()
		content = [x.strip() for x in content]
		#content = [x.split() for x in content]
	return content[:-2],content[-1].split()

def main(aaRate_path):
	rate_matrix,freq_vector = matrix_reader(aaRate_path)
	aarate_dict = dict()
	for v_counter,line in enumerate(rate_matrix):
		for h_counter, element in enumerate(line.split()):
			aarate_dict[(aa_tup_h[h_counter], aa_tup_v[v_counter])] = float(element)
			aarate_dict[(aa_tup_v[v_counter],aa_tup_h[h_counter])] = float(element)
	matrix_list = []
	for hnum,ah in enumerate(aa_tup):
		rowlist = []
		idsum = 0
		for vnum,av in enumerate(aa_tup):
			if av is not ah:
				rowlist.insert(vnum,aarate_dict[(av,ah)])
				#print(idsum,av,ah,aarate_dict[(av,ah)],vnum+1,hnum+1)
				idsum = idsum + aarate_dict[(av,ah)]*float(freq_vector[vnum])
		rowlist.insert(hnum,(idsum/float(freq_vector[hnum])))
		matrix_list.append(rowlist)
	df = pd.DataFrame(matrix_list, columns=aa_tup,index = aa_tup)
	#print(df)
	print(df.applymap(lambda x: math.log(x,2)))
	with open(sys.argv[2]) as f:
		my_mx = read_text_matrix(f)
	for k,v in my_mx:
		print(k,v,my_mx[(k,v)])

if __name__ == "__main__":
	import sys
	main(sys.argv[1])
