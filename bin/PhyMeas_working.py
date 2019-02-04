"""Calculate and visualize conservation between two groups of sequences from one alignment"""
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align, argparse, random, math
import numpy as np
from Bio import AlignIO
from collections import defaultdict, Counter
from Bio.SeqUtils import IUPACData
from AlignmentGroup import AlignmentGroup
from MatrixLoad import PAMLmatrix


###Argument block; might want to make it into a function
parser = argparse.ArgumentParser(description='Calculate and visualize conservation between two groups of sequences from one alignment')
parser.add_argument('alignment_path', help='Path to alignment file')
parser.add_argument('-s','--structure_paths', nargs='+', help='Paths to structure files, can be one or many.')
entropy_group = parser.add_mutually_exclusive_group()
entropy_group.add_argument('-e','--shannon_entropy', help='Use shannon entropy for conservation calculation.', action="store_true")
entropy_group.add_argument('-c','--reflected_shannon', help='Use shannon entropy for conservation calculation and reflect the result so that a fully random sequence will be scored as 0.', action="store_true")
structure_option = parser.add_mutually_exclusive_group()
structure_option.add_argument('-ss','--secondary_structure', help = 'Use substitution matrices derived from data dependent on the secondary structure assignment.', action="store_true")
structure_option.add_argument('-be','--burried_exposed', help = 'Use substitution matrices derived from data dependent on the solvent accessability of a residue.', action="store_true")
structure_option.add_argument('-ssbe','--both', help = 'Use substitution matrices derived from data dependent on both the secondary structure and the solvent accessability of a residue.', action="store_true")
commandline_args = parser.parse_args()
###

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
	default_aa_sequence = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
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
	if Counter(list(hash_AA.keys())) == Counter(default_aa_sequence):
		return default_aa_sequence
	else:
		raise ValueError("Alignment has non-standard AAs!")

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

def shannon_entropy(alnObject, aa_list, commandline_args):
	'''
	Function to calcuate the reflected Shannon entropy per alignment column.
	Fully random sequence will have reflected entropy of 0, while fully conserved
	column will be around 4.322
	H=-sum_{i=1}^{M} P_i,log_2,P_i	(http://imed.med.ucm.es/Tools/svs_help.html)
	Hrefl = abs(log_2,M) +sum_{i=1}^{M} P_i,log_2,P_i
	For gapped regions a random selection is made from the present AAs.
	'''
	entropy_alignment = AlignmentGroup(alnObject)
	AlignmentGroup.randomize_gaps(entropy_alignment, aa_list)
	alnObject = AlignmentGroup._return_alignment_obj(entropy_alignment)
	entropyDict={}
	for i in range(0, alnObject.get_alignment_length()):
		column_string=alnObject[:, i]
		unique_base = set(column_string)				# Get only the unique bases in a column
		M = len(column_string)							# Number of residues in column
		entropy_list = []
		for base in unique_base:
			n_i = column_string.count(base)				# Number of residues of type i
			P_i = n_i/float(M)							# n_i(Number of residues of type i) / M(Number of residues in column)
			entropy_i = P_i*(math.log(P_i,2))
			entropy_list.append(entropy_i)
		if commandline_args.reflected_shannon:
			refl_shentr = abs(math.log(1/len(aa_list),2))+sum(entropy_list)
			entropyDict[i]=refl_shentr
		elif commandline_args.shannon_entropy:
			sh_entropy = abs(sum(entropy_list))
			entropyDict[i]=sh_entropy
	return entropyDict

def compute_score(aln_index_dict, alngroup_dict):
	'''
	Computes score for each alignment position, 
	accounting for secondary structure.'''
	lg_matrix=PAMLmatrix('../test_data/LG.dat')
	ext_matrix=PAMLmatrix('../test_data/S.dat')
	hel_matrix=PAMLmatrix('../test_data/H.dat')
	oth_matrix=PAMLmatrix('../test_data/O.dat')
	
	lgmx=np.array(lg_matrix.lodd)
	extmx=np.array(ext_matrix.lodd)
	helmx=np.array(hel_matrix.lodd)
	othmx=np.array(oth_matrix.lodd)
	
	#Problematic code follows
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
			ss_indicator1=alngroup_dict[temp_column_dict[aln_index][0]][3][aln_index]
			ss_indicator2=alngroup_dict[temp_column_dict[aln_index][1]][3][aln_index]
			vr1=np.array(alngroup_dict[temp_column_dict[aln_index][0]][1][aln_index])
			vr2=np.array(alngroup_dict[temp_column_dict[aln_index][1]][1][aln_index])
			if ss_indicator1 == ss_indicator2:
				if ss_indicator1 == 'E' or ss_indicator1 == 'B':
					print(aln_index, vr1@extmx@vr2.T)
				elif ss_indicator1 == 'H' or ss_indicator1 == 'G':
					print(aln_index, vr1@helmx@vr2.T)
				else:
					print(aln_index, vr1@othmx@vr2.T)
			else:
				print(aln_index, vr1@lgmx@vr2.T)
			#print (aln_index, alngroup_dict[temp_column_dict[aln_index][0]][3][aln_index], '\n\t',temp_column_dict[aln_index][0], '\t',alngroup_dict[temp_column_dict[aln_index][0]][1][aln_index], '\n\t',temp_column_dict[aln_index][1], '\t',alngroup_dict[temp_column_dict[aln_index][1]][1][aln_index])
		else:								#In case of only one structure or no structure
			groupnames = list(alngroup_dict.keys())
			vr1=np.array(aln_index_dict[aln_index][groupnames[0]])
			vr2=np.array(aln_index_dict[aln_index][groupnames[1]])
			print(aln_index, vr1@lgmx@vr2.T)

def better_compute_score(aln_index_dict, list_groupnames, commandline_args, *args):
	'''Should handle all structure options'''
	lg_matrix=PAMLmatrix('../test_data/LG.dat')
	lgmx=np.array(lg_matrix.lodd)
	if len(args) == 0:
		for aln_index in aln_index_dict:
			vr1=np.array(aln_index_dict[aln_index][list_groupnames[0]])
			vr2=np.array(aln_index_dict[aln_index][list_groupnames[1]])
			print(aln_index, vr1@lgmx@vr2.T)
	if commandline_args.both:
		#Gotta add all matrices
		pass
	elif commandline_args.burried_exposed:
		#Gotta add burried/exposed matrices
		pass
	elif commandline_args.secondary_structure:
		ss_aln_index_map = args[0][0]
		ext_matrix=PAMLmatrix('../test_data/EXT.dat')
		hel_matrix=PAMLmatrix('../test_data/HEL.dat')
		oth_matrix=PAMLmatrix('../test_data/OTH.dat')
		extmx=np.array(ext_matrix.lodd)
		helmx=np.array(hel_matrix.lodd)
		othmx=np.array(oth_matrix.lodd)
		for aln_index in aln_index_dict:
			if aln_index in ss_aln_index_map[list_groupnames[0]] and aln_index in ss_aln_index_map[list_groupnames[1]]:
				ss_indicator1=ss_aln_index_map[list_groupnames[0]][aln_index]
				ss_indicator2= ss_aln_index_map[list_groupnames[1]][aln_index]
				vr1=np.array(aln_index_dict[aln_index][list_groupnames[0]])
				vr2=np.array(aln_index_dict[aln_index][list_groupnames[1]])
				
				


def main():
	'''Main entry point'''
	aln_path = commandline_args.alignment_path
	alignIO_out=read_align(aln_path)
	aa_list=uniq_AA_list(alignIO_out)
	sliced_alns = slice_by_name(alignIO_out)
	
	if commandline_args.shannon_entropy or commandline_args.reflected_shannon:
		entropyDict=shannon_entropy(alignIO_out, aa_list, commandline_args)
		for x in entropyDict:
			print(x, entropyDict[x])
		sys.exit(0)
	
	alngroup_dict={}							#Dictionary combining class outputs
	aln_index_dict=defaultdict(dict)			#Dictionary combining class outputs with first key the aln index and second key the group name

	if commandline_args.structure_paths is None:
		groupname_list=[]
		for alngroup_name in sliced_alns:
			groupname_list.append(alngroup_name)
			alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name])
			AlignmentGroup.randomize_gaps(alngroup_name_object, aa_list)
			alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
			for aln_index in alnindex_col_distr:
					aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
		better_compute_score(aln_index_dict, groupname_list)
		sys.exit(0)

	elif commandline_args.structure_paths:
		if len(commandline_args.structure_paths) < len(sliced_alns):
			print("Add code to handle case of fewer structures than groups in alignment")
			sys.exit(0)

		groupname_list=[]
		ss_aln_index_map = defaultdict(dict)
		res_depth_aln_index_map = defaultdict(dict)
		for alngroup_name in sliced_alns:
			groupname_list.append(alngroup_name)
			current_path = [s for s in commandline_args.structure_paths if alngroup_name in s]				#Add a check to ensure it is a single match and that there is a match!
			print(alngroup_name, current_path)
			alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name],current_path[0])
			struc_to_aln_index_mapping=AlignmentGroup.create_aln_struc_mapping(alngroup_name_object)
			AlignmentGroup.randomize_gaps(alngroup_name_object, aa_list)
			alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
			for aln_index in alnindex_col_distr:
					aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
			if commandline_args.secondary_structure:
				ss_aln_index_map[alngroup_name] = AlignmentGroup.ss_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				alngroup_dict[alngroup_name] = [struc_to_aln_index_mapping,ss_aln_index_map]
			elif commandline_args.burried_exposed:
				res_depth_aln_index_map[alngroup_name] = AlignmentGroup.depth_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				alngroup_dict[alngroup_name] = [struc_to_aln_index_mapping,res_depth_aln_index_map]
			elif commandline_args.both:
				ss_aln_index_map[alngroup_name] = AlignmentGroup.ss_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				res_depth_aln_index_map[alngroup_name] = AlignmentGroup.depth_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				#Fix this ridiculousness
				alngroup_dict[alngroup_name] = [struc_to_aln_index_mapping,alnindex_col_distr,AlignmentGroup._return_alignment_obj(alngroup_name_object),ss_aln_index_map,res_depth_aln_index_map]
			else:
				parser.print_help()
				raise ValueError("When a structure is defined, one of the matrix options are rquired!")
		structure_derived_data = []
		if commandline_args.both:
			structure_derived_data.append(ss_aln_index_map)
			structure_derived_data.append(res_depth_aln_index_map)
		elif commandline_args.burried_exposed:
			structure_derived_data.append(res_depth_aln_index_map)
		elif commandline_args.secondary_structure:
			structure_derived_data.append(ss_aln_index_map)
		
		better_compute_score(aln_index_dict, groupname_list, commandline_args, structure_derived_data)
		compute_score(aln_index_dict, alngroup_dict)


if __name__ == '__main__':
	sys.exit(main())

