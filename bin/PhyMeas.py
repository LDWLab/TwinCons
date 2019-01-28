"""Script to calculate and visualize conservation between two groups of sequences from one alignment"""
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align
from Bio import AlignIO
from collections import defaultdict
from Bio.SeqUtils import IUPACData
from AlignmentGroup import AlginmentGroup
from MatrixLoad import PhymlMatrix

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
	pdb_paths = []
	for pdb_path in sys.argv[2:]:
		pdb_paths.append(pdb_path)
	protein_struc_path = sys.argv[2]
	alignIO_out=read_align(aln_path)

	sliced_alns = slice_by_name(alignIO_out)
	aa_list=uniq_AA_list(alignIO_out)
	parsed_struct={}
	for pdb_file in pdb_paths:
		pdbDict, linesList, residueList=parse_pdb(pdb_file)
		parsed_struct[pdb_file]=[pdbDict, linesList, residueList]

	alngroup_dict={}							#Dictionary combining class outputs
	aln_index_dict=defaultdict(dict)			#Dictionary combining class outputs with first key the aln index
	for alngroup_name in sliced_alns:
		current_path = [s for s in pdb_paths if alngroup_name in s]				#Add a check to ensure it is a single match and that there is a match!
		print(alngroup_name, current_path)
		###
		alngroup_name_object = AlginmentGroup(sliced_alns[alngroup_name],current_path[0])
		struc_to_aln_index_mapping=AlginmentGroup.create_aln_struc_mapping(alngroup_name_object)
		AlginmentGroup.randomize_gaps(alngroup_name_object, aa_list)
		alnindex_col_distr = AlginmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
		ss_aln_index_map,res_depth_aln_index_map = AlginmentGroup.ss_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
		###
		alngroup_dict[alngroup_name] = [struc_to_aln_index_mapping,alnindex_col_distr,AlginmentGroup._return_alignment_obj(alngroup_name_object),ss_aln_index_map,res_depth_aln_index_map]
		
		for aln_index in alnindex_col_distr:
			aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
		

	for aln_index in aln_index_dict:
		#print(aln_index, aln_index_dict[aln_index])
		temp_column_dict=defaultdict(dict)
		for groupname in aln_index_dict[aln_index]:
			if aln_index in alngroup_dict[groupname][3]:
				if aln_index in temp_column_dict:
					temp_column_dict[aln_index].append(groupname)
				else:
					temp_column_dict[aln_index]=[]
					temp_column_dict[aln_index].append(groupname)
				#print(groupname, aln_index, alngroup_dict[groupname][3][aln_index], alngroup_dict[groupname][4][aln_index],aln_index_dict[aln_index][groupname])
		if len(temp_column_dict[aln_index]) > 1:
			if alngroup_dict[temp_column_dict[aln_index][0]][3][aln_index] == alngroup_dict[temp_column_dict[aln_index][1]][3][aln_index]:
				print (aln_index, alngroup_dict[temp_column_dict[aln_index][0]][3][aln_index], '\n\t',temp_column_dict[aln_index][0], '\t',alngroup_dict[temp_column_dict[aln_index][0]][1][aln_index], '\n\t',temp_column_dict[aln_index][1], '\t',alngroup_dict[temp_column_dict[aln_index][1]][1][aln_index])

	lg_matrix=PhymlMatrix('../test_data/LG.dat')
	lg_matrix.calculate_Sij()

if __name__ == '__main__':
	sys.exit(main())

