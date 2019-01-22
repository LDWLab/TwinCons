"""Script to calculate and visualize conservation between two groups of sequences from one alignment"""
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align
from Bio import AlignIO
from Bio.SeqUtils import IUPACData
from AlignmentGroup import AlginmentGroup

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
	return sliced_dict			#Iterate over the dict and create instances of AlignmentGroup

def parse_pdb(file):
	'''
	Simplifies a pdb file into 3 objects
	'''
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

	sliced_alns = slice_by_name(alignIO_out)
	aa_list=uniq_AA_list(alignIO_out)
	pdbDict, linesList, residueList=parse_pdb(protein_struc_path)

	for alngroup_name in sliced_alns:
		print(alngroup_name)
		alngroup_name = AlginmentGroup(sliced_alns[alngroup_name],protein_struc_path)
		struc_to_aln_index_mapping=AlginmentGroup.create_struc_aln_mapping(alngroup_name)
		AlginmentGroup.randomize_gaps(alngroup_name, aa_list)
		print(struc_to_aln_index_mapping,AlginmentGroup._return_alignment_obj(alngroup_name))


if __name__ == '__main__':
	sys.exit(main())

