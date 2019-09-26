#!/usr/bin/env python3
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align
import numpy as np
from Bio import AlignIO
from Bio.Seq import MutableSeq
from Bio.SubsMat import MatrixInfo
from Bio.SeqUtils import IUPACData
from Bio import BiopythonWarning
from numpy.random import choice
from collections import Counter
from textwrap import wrap
import networkx as nx
import scipy as sp



def process_align(aln_obj):
	vert_dict={}
	for y in aln_obj:
		i=0
		vert_dict[y.id] = {} 
		for x in y:
			i+=1
			vert_dict[y.id][i] = x
	return vert_dict

def read_align(aln_name):
	'''
	Reads the fasta file and gets the sequences.
	'''
	alignments = AlignIO.read(open(aln_name), "fasta")
	return alignments
	
def cut_gaps(align,gap_perc):
	'''
	Cuts columns of the MSA with too many gaps (specififed with -g)
	'''
	n = float(len(align[0]))
	i = 0
	x=0
	length=1
	gap_mapping={}
	while i < n:
		x=align[:, i].count('-')/len(align)					#Get percentage of gaps in column
		if x > float(gap_perc):
			if i == 0:
				align = align[:, 1:]
			elif i+1 == n:
				align = align[:, :i]
			else:
				align = align[:, :i] + align[:, i+1:]
			n -= 1											#  seq. 1 shorter
		else:												#  nothing to delete, proceed
			i += 1
		gap_mapping[int(length)]=i
		length+=1
	return gap_mapping,align,len(align[0])

def randomize_gaps(align, aa_list):
	'''
	Substitutes gaps in the alignment with a random choice from the present AAs.
	Should be an option to either use absolute random or random sampled from the distribution of the sequence.
	'''
	for aln in align:
		i = 0
		newaln=MutableSeq(str(aln.seq))
		aln_no_dash = str(aln.seq).replace('-', '')
		distribution = Counter(aln_no_dash)
		choice_distr = []
		for aa in aa_list:
			choice_distr.append(distribution[aa]/len(aln_no_dash))
		for resi in aln.seq:
			if resi == '-':
				newaln[i] = choice(aa_list, p=choice_distr)
				#newaln[i] = random.choice(aa_list)
			i+=1
		aln.seq=newaln
	return align

def slice_by_name(unsliced_aln_obj):
	'''
	Slices an alignment into different alignments
	by first string (should hold protein name)
	Returns a dictionary with keys being the alignment
	names pointing to the alignIO objects
	'''
	prot_list=[]
	sliced_dict={}
	for entry in unsliced_aln_obj:
		prot_list.append(entry.id.split("_")[0])
	uniq_prot_list=set(prot_list)
	for prot in uniq_prot_list:
		what = Bio.Align.MultipleSeqAlignment([])
		for entry in unsliced_aln_obj:
			if re.match(prot,entry.id.split("_")[0]):
				what.append(entry)
		sliced_dict[prot]=what
	return sliced_dict

def uniq_AA_list(aln_obj):
	'''
	Creates list of unique AA residues in the given MSA to be used for frequency iterator.
	Also checks if the alignment has AA letters from the IUPAC extended_protein_letters.
	'''
	hash_AA=dict()
	for alignment in aln_obj:
		for amac in alignment.seq:
			if re.match(r'-',amac):
				pass
			else:
				hash_AA[amac]='null'
	if all (x in IUPACData.extended_protein_letters for x in hash_AA.keys()):
		pass
	else:
		raise ValueError("Alignment has AA letters not found in the IUPAC extended list!")
	return list(hash_AA.keys())

def freq_iterator(column, aa_list,gap_perc):
	'''
	Calculates frequency of each AA in the column.
	'''
	col_aalist=[]
	for aa in aa_list:
		#print(aa, column.count(aa)/len(column.replace("-", "")))
		if gap_perc == '1':
			col_aalist.append(column.count(aa)/len(column))
		else:
			col_aalist.append(column.count(aa)/len(column.replace("-", "")))
	#print()
	return col_aalist

def edge_helper(G,my_mx):
	'''
	Helper function for edge construction of the two networks.
	'''
	iter_nodes=list()
	iter_nodes=list(itertools.combinations_with_replacement(list(G.nodes()), 2))
	for node_pair in iter_nodes:
		if node_pair in getattr(MatrixInfo,my_mx):
			G.add_edge(node_pair[0], node_pair[1], weight=getattr(MatrixInfo,my_mx)[node_pair])
		else:
			G.add_edge(node_pair[0], node_pair[1], weight=getattr(MatrixInfo,my_mx)[node_pair[1],node_pair[0]])
	return True

def baseline_calc(aa_list, my_mx):
	'''
	Function to calculate the baseline of a given substitution matrix
	Returns integer value which will be added to each element of the
	substitution matrix. This way a random sequence should have score
	of 0.
	'''
	testG=nx.Graph()
	for aa in aa_list:
		testG.add_node(aa)
	edge_helper(testG, my_mx)
	testA = nx.to_numpy_matrix(testG)
	testvr = np.repeat(1/len(aa_list),len(aa_list))
	baseline = float(testvr@testA@testvr.T)
	revtestA=np.add(testA, abs(baseline))
	if int(testvr@revtestA@testvr.T) != 0:
		raise ValueError("Wasn't able to baseline the substitution matrix correctly. Please use BL62")
	else:
		return baseline



def compute_score(perclist1, perclist2, aalist,my_mx, baseline):
	G1 = nx.Graph()
	G2 = nx.Graph()
	i=0
	for aa in aalist:
		#print (aa, i, perclist1[i], perclist2[i])
		if perclist1[i] != 0:
			G1.add_node(aa,value=perclist1[i])
			G2.add_node(aa,value=perclist2[i])
		if perclist2[i] != 0:
			G2.add_node(aa,value=perclist2[i])
			G1.add_node(aa,value=perclist1[i])
		i+=1
	edge_helper(G1,my_mx)
	edge_helper(G2,my_mx)
	A1 = nx.to_numpy_matrix(G1)
	A1 = np.add(A1,abs(baseline))
	#A2 = nx.to_numpy_matrix(G2)				#A1 and A2 are identical

	G1valist=[]
	G2valist=[]
	for x in G1.nodes():
		G1valist.append(nx.get_node_attributes(G1,'value')[x])
		G2valist.append(nx.get_node_attributes(G2,'value')[x])
	A1=np.array(A1)
	#A2=np.array(A2)
	G1vr=np.array(G1valist)[np.newaxis]
	G2vr=np.array(G2valist)[np.newaxis]
	#print('priceG1->G2\n',G2vr,'\n',A1,'\n',G1vr)
	#print('price\t',float(G1vr@A1@G2vr.T))
	return float(G1vr@A1@G2vr.T)

def build_stats(dict):
	'''
	Builds statistics for the entire alignment
	'''
	scorelist=[]
	for x in dict:
		scorelist.append(dict[x])
	print(np.average(scorelist), np.std(scorelist))

def perresi(alen,sliced_dict,aa_list,gp_map2,my_mx,baseline,gap_perc):
	out_dict={}
	i=0
	while i < alen:
		conscorcol_dict={}
		consnam_list=[]
		col_aalist=[]
		for aln_g in sliced_dict:
			#print (gp_map2[i], aln_g)
			col_aalist=freq_iterator(sliced_dict[aln_g][:, i], sorted(aa_list),gap_perc)
			conscorcol_dict[aln_g]=col_aalist
			consnam_list.append(aln_g)
		iter_nams=list(itertools.combinations(consnam_list, 2))
		i+=1
		if len(consnam_list) > 2:
			raise ValueError('Please use only two groups within single alignment.')
		for name_pair in iter_nams:
			#print (gp_map2[i],name_pair[0], name_pair[1])
			pos_score=compute_score(conscorcol_dict[name_pair[0]], conscorcol_dict[name_pair[1]],sorted(aa_list),my_mx,baseline)
			out_dict[gp_map2[i]]=(pos_score)
			#print (gp_map2[i], conscorcol_dict[name_pair[0]], conscorcol_dict[name_pair[1]])
	return out_dict, consnam_list

def main(aln_path,my_mx,gap_perc):
	'''
	Calulate conservation score based on a substitution matrix
	'''
	alignIO_out=read_align(aln_path)
	aa_list=uniq_AA_list(alignIO_out)
	if gap_perc == '1':
		alignIO_out=randomize_gaps(alignIO_out,aa_list)
	gp_mapping,cutted_aln,alen=cut_gaps(alignIO_out,gap_perc)
	gp_map2 = dict((y,x) for x,y in sorted(gp_mapping.items(), reverse=True))
	sliced_dict=slice_by_name(cutted_aln)
	baseline=baseline_calc(aa_list,my_mx)
	out_dict,consnam_list = perresi(alen,sliced_dict,aa_list,gp_map2,my_mx,baseline,gap_perc)
	#build_stats(out_dict)
	return out_dict,consnam_list

