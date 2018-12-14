"""Script to calculate and visualize conservation between two groups of sequences from one alignment"""
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align
from Bio import AlignIO
from Bio.SeqUtils import IUPACData

def read_align(aln_path):
	'''
	Reads the fasta file and gets the sequences.
	'''
	alignments = AlignIO.read(open(aln_path), "fasta")
	return alignments

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

def get_align(file, alignments):
	'''
	Make a dictionary complementing the locations in 
	the alignment fasta file and the structure pdb file
	'''
	alignDict = {}
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
	return dictList

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

def main():
	"""Main entry for the script"""
	aln_path = sys.argv[1]
	protein_struc_path = sys.argv[2]
	alignIO_out=read_align(aln_path)
	aa_list=uniq_AA_list(alignIO_out)
	struc_to_aln_index_mapping=get_align(protein_struc_path, alignIO_out)
	pdbDict, linesList, residueList=replace_bfactor(protein_struc_path)
	print(pdbDict, residueList)

if __name__ == '__main__':
	sys.exit(main())