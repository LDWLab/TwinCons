#!/usr/bin/env python3
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align
import numpy as np
import getopt
from Bio import AlignIO
from Bio.Seq import MutableSeq
from Bio.SubsMat import MatrixInfo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.SeqUtils import IUPACData
import seaborn
from textwrap import wrap
import networkx as nx
import scipy as sp
from scipy.spatial import distance
from Bio import BiopythonWarning
from numpy.random import choice
from collections import Counter
from single_cons_comp_old import *


def usage():
	print (\
	"USAGE:\n./entropy.py -a [alignment_file_path] -p [protein_struc_path] -o [output_file_path] -h\n\
	-a: defines path to alignment file. Works only on fasta type of alignments.\tREQUIRED\n\
	-o: defines output image path. Default is ./test_set/output/alignment_file_name.png\n\
	-g: columns with gaps greater than this percentage will be removed; default 0.3(decimal between 0 and 1)\n\
	-m: what substitution matrix to use (blosum45, pam50, etc.) Default is blosum62.\n\
	-h: prints this\
")

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:o:g:m:h', ['alignment=', 'output=', 'percentage=', 'matrix=', 'help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

#output_path='test.png'
my_mx='blosum62'
gap_perc=0.3
for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-a', '--alignment'):
		aln_path = arg
		if re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',aln_path):
			pass
		else:
			print("Fasta file should be in its own folder")
			usage()
			sys.exit(2)
	elif opt in ('-o', '--output'):
		output_path = arg
	elif opt in ('-g'):
		gap_perc = float(arg)
	elif opt in ('-m'):
		my_mx = arg
	else:
		usage()
		sys.exit(2)

#Set up default output for image
try:
	output_path
except:
	output_path = "./test_set/output/"+str(re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',aln_path)[0][1])+".png"

name_file = str(re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',aln_path)[0][1])+str(re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',aln_path)[0][2])

def makeimagenozero(out_dict,group_names):
	'''
	Prints out a scatter of the scores.
	'''
	seaborn.set(style='ticks')
	ylist=[]
	xlist=[]
	for pos in sorted(out_dict.keys()):
		xlist.append(pos)
		ylist.append(out_dict[pos])
	ax = plt.subplot()
	#rint(type(ylist))
	negative_mask = np.array(ylist) < 0
	positive_mask = np.array(ylist) >= 0
	plt.bar(np.array(xlist)[positive_mask],np.array(ylist)[positive_mask], 2,color='red', edgecolor=['black']*len(xlist))
	plt.bar(np.array(xlist)[negative_mask],np.array(ylist)[negative_mask], 2,color='blue', edgecolor=['black']*len(xlist))
	plt.yticks(np.arange(min(getattr(MatrixInfo,my_mx).values()),max(getattr(MatrixInfo,my_mx).values()), step=1))
	ax.set_ylim(min(getattr(MatrixInfo,my_mx).values()),max(getattr(MatrixInfo,my_mx).values()))
	if True in negative_mask and True in positive_mask:
		plt.legend(['more likely than random', 'less likely than random'], loc='upper left')
	plt.xlabel('Alignment position')
	plt.ylabel('Transformation score')
	title = "Alignment file "+name_file+" between "+group_names[0]+" and "+group_names[1]+", substitution matrix "+my_mx
	plt.title("\n".join(wrap(title, 60)))
	ax.grid(True, which='both')
	#seaborn.despine(ax=ax, offset=0)
	plt.savefig(output_path, dpi=600)

def main():
	alignIO_out=read_align(aln_path)
	aa_list=uniq_AA_list(alignIO_out)
	if gap_perc == 1:
		alignIO_out=randomize_gaps(alignIO_out,aa_list)
	gp_mapping,cutted_aln,alen=cut_gaps(alignIO_out,gap_perc)
	gp_map2 = dict((y,x) for x,y in sorted(gp_mapping.items(), reverse=True))
	sliced_dict=slice_by_name(cutted_aln)
	baseline=baseline_calc(aa_list,my_mx)
	out_dict,consnam_list = perresi(alen,sliced_dict,aa_list,gp_map2,my_mx,baseline,gap_perc)
	#build_stats(out_dict)
	build_stats(out_dict)
	
	makeimagenozero(out_dict,consnam_list)

if __name__ == "__main__":
	main()
