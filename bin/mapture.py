#!/usr/bin/env python3
import os, re, sys, math, copy, random, getopt, warnings, statistics, matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import single_cons_comp
import matplotlib.pyplot as plt
from Bio.SeqUtils import IUPACData
from biopandas.pdb import PandasPdb
from Bio import AlignIO, BiopythonWarning

def usage():
	print (\
	"USAGE:\nmapture.py -a [alignment_file_path] -p [protein_struc_path] -o [output_file_path] -h\n\
	-a: defines path to alignment file. Works only on fasta type of alignments. Must have the same specie as the one in pdb name.\tREQUIRED\n\
	-p: defines path to pdb structure file. Has to be of the format [path]/[5 letter specie name]_[protein_name].pdb\t\tREQUIRED\n\
	-o: defines output pdb or pml path. \n\
	-m: defines substitution matrix to use (blosum45, pam50, etc.) Default is blosum62.\n\
	\tif set to 1 will substitute gaps with random choice of aa letters.\n\
	-t: modifying the bfactor column by the following options (choose one):\t\t\t\t\t\t\t\tREQUIRED\n\
	\t gaps:\tgaps in the column of the alignment file.\n\
	\t entr:\tShannon's entropy for each column.\n\
	\t cons:\tConservation for each column (reflected Shannon).\n\
	\t blos:\tTransformation score for each column (requires two sets of groups in the alignment).\n\
	-h: prints this\
	")

my_mx='blosum62'

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:p:o:m:t:h', ['alignment=', 'structure=', 'output=', 'matrix=', 'option=','help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-a', '--alignment'):
		if re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',arg):
			aln_path = arg
		else:
			usage()
			raise TypeError ("Input alignment file should be fasta format!")
			sys.exit(2)
	elif opt in ('-p', '--protein_structure'):
		protein_struc_path = arg
	elif opt in ('-m', '--matrix'):
		my_mx = arg
	elif opt in ('-o', '--output'):
		output_path = arg
	elif opt in ('-t', '--option'):
		option = arg
	else:
		usage()
		sys.exit(2)
		
#Add check for options


#Set up default output for pdb
try:
	output_path
except:
	output_path = "./output/"+str(re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',aln_path)[0][1])+".pdb"
		
def read_align(aln_path):
	alignments = AlignIO.read(open(aln_path), "fasta") 		# Reading the Fasta file
	return alignments 										# getting the sequences
	

alignDict = {}
def get_align(file, alignments):
	'''
	Make a dictionary complementing the locations in 
	the alignment fasta file and the structure pdb file
	'''
	file = re.findall(r'(.*\/)(.*)(_.*.pdb)',file)[0][1]
	for alignment in alignments: 							# Iterating through each sequence
		alignDict[alignment.id] = alignment.seq 			# Making a dictionary of the sequences
	dictseq={} 												# gives the complement between the location in alns (key) and the location in pdb (value)
	i=0
	a=1
	for alignment in alignments:
		if re.search(file, alignment.id) is not None:
			anchor_seq=alignment 							# anchor_seq will hold the sequence data from the name of the pdb file
			for x in anchor_seq:
				if (x=="-"):
					dictseq[a] = []
					dictseq[a].append(0)
					a+=1
				else:
					i+=1
					dictseq[a] = []
					dictseq[a].append(i)
					a+=1
	dictList={} 							# gives the complement between the location in pdb (key) and the location in alns (value) after removing gaps
	for k,v in dictseq.items():
		if v[0]==0:
			next
		else:
			dictList[v[0]]=k
	return(dictseq, dictList)

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

def replace_bfactor(file):	
	linesList = [line.split() for line in open(file)]	# List of all the lines in the pdb file
	linesList=linesList[:-1]
	residueList = [] 									# List of the residue numbers
	pdbDict={} 											# Has items in col 6 (col 23-26 acc to pdb format) as the key and items in col 4 (columns 18-20) as value
	for items in linesList:
		pdbDict[items[5]]=items[3]
	for i in linesList:
		if i[5] in residueList:
			next
		else:
			residueList.append(i[5])
	return pdbDict, linesList, residueList	
	
def vertical_residues(alignIO_out, pdbDict, dictList):
	'''
	Function that takes in an alignment object, pdbDict, dictList.
	It returns a dictionary with key the position from the anchor sequence and
	values a list of all residues in this column.
	'''
	alnDict={} 						# Gives the struc location and the aa residues at that location from the alignments file vertically. The values (residues) are in a list
	alnList=[] 						# has the struc location and the aa residue for each sequence at that location
	for alignment in alignIO_out:
		for key in pdbDict:
			if int(key) <= len(dictList):
				temp = key,(alignment[int(dictList[int(key)]-1)])
				alnList.append(temp)
			else:
				next
	for x in alnList:
		if x[0] in alnDict.keys():
			alnDict[x[0]].append(x[1])
		else:
			alnDict[x[0]]=[x[1]]
	return alnDict


# Source code for this function-https://github.com/ffrancis/Multiple-sequence-alignment-Shannon-s-entropy/blob/master/msa_shannon_entropy012915.py
def shannon_entropy(alnDict, aa_list, option):
	'''
	Function to calcuate the reflected Shannon entropy per alignment column.
	Fully random sequence will have reflected entropy of 0, while fully conserved
	column will be around 4.322
	H=-sum_{i=1}^{M} P_i,log_2,P_i	(http://imed.med.ucm.es/Tools/svs_help.html)
	Hrefl = abs(log_2,M) +sum_{i=1}^{M} P_i,log_2,P_i
	For gapped regions a random selection is made from the present AAs.
	'''
	entropyDict={}
	for k,v in alnDict.items():
		i=0
		z=copy.deepcopy(v)
		for aa in v:
			if aa == '-':
				z[i]=random.choice(aa_list)
			i+=1
		unique_base = set(z)							# Get only the unique bases in a column
		M = len(z)
		entropy_list = []								# Number of residues in column
		for base in unique_base:
			n_i = z.count(base)							# Number of residues of type i
			P_i = n_i/float(M)							# n_i(Number of residues of type i) / M(Number of residues in column)
			entropy_i = P_i*(math.log(P_i,2))
			entropy_list.append(entropy_i)
		
		if option == 'cons':
			refl_shentr = abs(math.log(1/len(aa_list),2))+sum(entropy_list)
			entropyDict[k]=refl_shentr
		elif option == 'entr':
			sh_entropy = abs(sum(entropy_list))
			entropyDict[k]=sh_entropy
		else:
			refl_shentr = abs(math.log(1/len(aa_list),2))+sum(entropy_list)
			entropyDict[k]=refl_shentr
	return entropyDict

def blscore_builder(blscor_dict, alnDict,dictList):
	'''
	Flip the dictionary used for the blscore using the anchor-alignment position dictionary.
	'''
	blscor_resn={}
	anpos_alpos = {v:k for k,v in dictList.items()}
	for k in blscor_dict:
		if k in anpos_alpos:
			blscor_resn[anpos_alpos[k]]=blscor_dict[k]
	return blscor_resn

def templist_builder(ppdb, option, alnDict, entropyDict, blscor_resn):
	entr_templist=[]
	gaps_templist=[]
	blos_templist=[]
	for resn in ppdb.df['ATOM']['residue_number']:
		i=0
		for v in alnDict[str(resn)]:
			if v == '-':
				i+=1
		gaps_templist.append(i)
		entr_templist.append(entropyDict[resn])
		if resn in blscor_resn.keys():
			blos_templist.append(blscor_resn[resn])
		else:
			blos_templist.append(0)
	return gaps_templist, entr_templist, blos_templist

def plotter (entropyDict, blscor_resn):
	plot_entr=[]
	plot_scor=[]
	for x in sorted(entropyDict.keys(), key=abs):
		#print(x, blscor_resn[int(x)])
		if x in blscor_resn.keys():
			#print (x, blscor_resn[int(x)], entropyDict[x])
			plot_entr.append(entropyDict[x])
			plot_scor.append(blscor_resn[int(x)])
		else:
			raise ValueError ("Non-equal entropy and bl score lists!")
	ax = plt.subplot()
	sns.set(style='ticks')
	plt.plot(plot_scor, label="Transformation score", linewidth=1)
	plt.plot(plot_entr, label="Conservation", linewidth=1)
	plt.legend()
	ax.grid(True, which='both')
	sns.despine(ax=ax, offset=0)
	plt.savefig('./test.png', dpi=600)

def main():
	alignIO_out=read_align(aln_path)
	aa_list=uniq_AA_list(alignIO_out)
	dictseq, dictList=get_align(protein_struc_path, alignIO_out)
	pdbDict, linesList, residueList=replace_bfactor(protein_struc_path)
	alnDict=vertical_residues(alignIO_out, pdbDict, dictList)
	entropyDict = shannon_entropy(alnDict, aa_list,option)
	out_dict,consnam_list=single_cons_comp.main(aln_path,my_mx)
	blscor_resn = blscore_builder(out_dict, alnDict,dictList)
	entropyDict = {int(k):v for k,v in entropyDict.items()}
	
	#Create lineplot figure for conservation and transformation score
	plotter(entropyDict, blscor_resn)
	
	if '.pml' in output_path:
		print("Working here")

	elif '.pdb' in output_path:
		# FOLLOWING NEEDS TO BE A FUNCTION
		# Using BioPandas to modify a column (bfactor)in the pdb structure file.
		# Source-https://rasbt.github.io/biopandas/tutorials/Working_with_PDB_Structures_in_DataFrames/
		ppdb = PandasPdb().read_pdb(protein_struc_path)  		#load pdb files from local directories
		#Option checks and building the bfactor field list
		gaps_templist, entr_templist, blos_templist = templist_builder(ppdb, option, alnDict, entropyDict, blscor_resn)

		#Then substitute this list into the dataframe column called b_factor 
		if option == 'gaps':									# Checking for gaps in the vertical column
			ppdb._df['ATOM']['b_factor'] = gaps_templist
		elif option == "entr" or option == "cons":				# Shannon's entropy (or conservation) for each column
			ppdb._df['ATOM']['b_factor'] = entr_templist
		elif option == "blos":									# Transformation score between two groups
			ppdb._df['ATOM']['b_factor'] = blos_templist
		
		#Finally output the pdb file
		ppdb.to_pdb(path=output_path,
						records=None,
						gz=False,
						append_newline=True)

		#Print some success
		print("Wrote "+output_path+" with option "+option+" for alignment "\
		+aln_path+" and groups "+consnam_list[0]+" "+consnam_list[1])
	else:
		usage()
		sys.exit(2)

if __name__ == "__main__":
	main()
