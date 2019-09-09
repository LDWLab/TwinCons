"""Calculate and visualize conservation between two groups of sequences from one alignment"""
import re, sys, warnings, statistics, copy, itertools, random, Bio.Align, argparse, random, math, matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from Bio import AlignIO
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
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
		print(list(hash_AA.keys()))
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
			entropyDict[i+1]=refl_shentr
		elif commandline_args.shannon_entropy:
			sh_entropy = abs(sum(entropy_list))
			entropyDict[i+1]=sh_entropy
	return entropyDict

def struc_anno_matrices (struc_anno):
	'''Returns a log odds matrix from a given name of a PAML type matrix'''
	return np.array(PAMLmatrix('../test_data/'+struc_anno+'.dat').lodd)

def compute_score(aln_index_dict, *args):
	'''
	Computes transformation score between two groups, using substitution 
	matrices on the common structural elements between two groups.
	Returns a dictionary with key the alignment index and value the computed score.
	'''
	alnindex_score=defaultdict(dict)
	if commandline_args.structure_paths:
		struc_annotation = args[0]
		groupnames = list(struc_annotation.keys())
		for aln_index in aln_index_dict:
			vr1 = np.array(aln_index_dict[aln_index][groupnames[0]])
			vr2 = np.array(aln_index_dict[aln_index][groupnames[1]])
			if aln_index in struc_annotation[groupnames[0]] and aln_index in struc_annotation[groupnames[1]]:
				common_chars = sorted(set(struc_annotation[groupnames[0]][aln_index]) & set (struc_annotation[groupnames[1]][aln_index]))
				if len(common_chars) > 0:
					alnindex_score[aln_index] = vr1@struc_anno_matrices(''.join(common_chars))@vr2.T
					#print(aln_index, vr1@struc_anno_matrices(''.join(common_chars))@vr2.T)
				else:
					lgmx = np.array(PAMLmatrix('../test_data/LG.dat').lodd)
					alnindex_score[aln_index] = vr1@lgmx@vr2.T
					#print(aln_index,vr1@lgmx@vr2.T)
			else:
				lgmx = np.array(PAMLmatrix('../test_data/LG.dat').lodd)
				alnindex_score[aln_index] = vr1@lgmx@vr2.T
				#print(aln_index,vr1@lgmx@vr2.T)
	else:											#Case of no structure defined inputs
		groupnames = args[0]
		for aln_index in aln_index_dict:
			vr1 = np.array(aln_index_dict[aln_index][groupnames[0]])
			vr2 = np.array(aln_index_dict[aln_index][groupnames[1]])
			lgmx = np.array(PAMLmatrix('../test_data/LG.dat').lodd)
			alnindex_score[aln_index] = vr1@lgmx@vr2.T
			#print(aln_index,vr1@lgmx@vr2.T)
	return alnindex_score

def lookahead(iterable):
	"""Pass through all values from the given iterable, augmented by the
	information if there are more values to come after the current one
	(True), or if it is the last value (False).
	"""
	# Get an iterator and pull the first value.
	it = iter(iterable)
	last = next(it)
	# Run the iterator to exhaustion (starting from the second value).
	for val in it:
		# Report the *previous* value (more to come).
		yield last, True
		last = val
	# Report the last value.
	yield last, False

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
	sns.set(style='ticks',rc={'figure.figsize':(20,15)})
	plt.plot(plot_scor, label="Transformation score", linewidth=0.5)
	plt.plot(plot_entr, label="Conservation", linewidth=0.5)
	
	plt.legend()
	ax.grid(True, which='both')
	sns.despine(ax=ax, offset=0)
	dpi_scaling = 3*len(blscor_resn)
	plt.savefig('./test.svg', dpi=dpi_scaling)

def plotter2(out_dict,group_names):
	'''
	Prints out a scatter of the scores.
	'''
	sns.set(style='ticks')
	ylist=[]
	xlist=[]
	for pos in sorted(out_dict.keys()):
		xlist.append(pos)
		ylist.append(out_dict[pos])
	f = plt.figure(figsize=(20,15))
	ax = f.add_subplot(111)
	negative_mask = np.array(ylist) < -1
	null_mask=[]
	for x in ylist:
		if x > -1 and x < 1:
			null_mask.append(True)
		else:
			null_mask.append(False)
	positive_mask = np.array(ylist) >= 1
	plt.bar(np.array(xlist)[positive_mask],np.array(ylist)[positive_mask], 1,edgecolor='None',color='red')
	plt.bar(np.array(xlist)[null_mask],np.array(ylist)[null_mask], 1,edgecolor='None',color='gray')
	plt.bar(np.array(xlist)[negative_mask],np.array(ylist)[negative_mask], 1,edgecolor='None',color='blue')
	#plt.yticks(np.arange(min(getattr(MatrixInfo,my_mx).values()),max(getattr(MatrixInfo,my_mx).values()), step=1))
	#ax.set_ylim(min(getattr(MatrixInfo,my_mx).values()),max(getattr(MatrixInfo,my_mx).values()))
	if True in negative_mask and True in positive_mask:
		plt.legend(['more likely than random', 'random','less likely than random'], loc='upper left')
	plt.xlabel('Alignment position')
	plt.ylabel('Transformation score')
	title = "Alignment file "+commandline_args.alignment_path+" between "+group_names[0]+" and "+group_names[1]
	plt.title("\n".join(wrap(title, 60)))
	ax.grid(True, which='both')
	sns.despine(ax=ax, offset=0)
	dpi_scaling = 3*len(out_dict)
	plt.savefig('./test.svg', dpi=dpi_scaling)

#Use this one - works with negative values
#Important: scaling should be done so that colors are comparable between negative and  /All this might be wrong
#positive values... max(data)+min(abs(data)) (maybe negatives will be too white then)? /see how it looks
def upsidedown_horizontal_gradient_bar(out_dict,group_names):
	fig, ax = plt.subplots()
	data=[]
	for x in sorted(out_dict.keys()):
		data.append(out_dict[x])
	bar = ax.bar(range(len(data)),data)
	def gradientbars(bars):
		ax = bars[0].axes
		lim = ax.get_xlim()+ax.get_ylim()
		for bar in bars:
			bar.set_zorder(1)
			bar.set_facecolor("none")
			x,y = bar.get_xy()
			w, h = bar.get_width(), bar.get_height()
			if h > 0:
				grad = np.atleast_2d(np.linspace(0,h/max(data),256)).T
				ax.imshow(grad, extent=[x,x+w,y,y+h], cmap=plt.get_cmap('Blues'), aspect="auto", zorder=0, norm=matplotlib.colors.NoNorm(vmin=0,vmax=1))
			else:			# can add different gradient for negative values
				grad = np.atleast_2d(np.linspace(0,h/min(data),256)).T
				ax.imshow(grad, extent=[x,x+w,y,y+h], cmap=plt.get_cmap('Reds'), aspect="auto", zorder=0, norm=matplotlib.colors.NoNorm(vmin=0,vmax=1))
		#ax.set_facecolor('Gray')
		ax.axis(lim)
	gradientbars(bar)
	dpi_scaling = 3*len(out_dict)
	plt.savefig('./test.svg',dpi=dpi_scaling)

def uninterrupted_stretches(alnindex, alnindex_score):
	"""Calculates lengths of uninterrupted lengths of positive and negative scores;
	Also associates these scores with the last position of the alignment index.
	For now uses > 0 as positive and < 0 as negative. Can be a variable or changed as appropriate.
	"""
	posdata={}
	negdata={}
	pos=0
	neg=0
	for x,has_more in lookahead(range(0,len(alnindex))):
		if alnindex_score[alnindex[x]] > 1:
			#print('pos',alnindex[x], alnindex_score[alnindex[x]])
			pos+=1
			if neg != 0:
				negdata[alnindex[x-1]]=neg
				neg = 0
		elif alnindex_score[alnindex[x]] < -1:
			#print('neg',alnindex[x], alnindex_score[alnindex[x]])
			neg+=1
			if pos != 0:
				posdata[alnindex[x-1]]=pos
				pos = 0
		else:				#in case of using some range between positive and negative scores for random
			#print('rand',alnindex[x], alnindex_score[alnindex[x]])
			if pos != 0:
				posdata[alnindex[x-1]]=pos
				pos = 0
			if neg != 0:
				negdata[alnindex[x-1]]=neg
				neg = 0
		if has_more is False:
			if pos != 0:
				posdata[alnindex[x]]=pos
			if neg != 0:
				negdata[alnindex[x]]=neg
	return posdata, negdata

def main():
	'''Main entry point'''
	aln_path = commandline_args.alignment_path
	alignIO_out=read_align(aln_path)
	aa_list=uniq_AA_list(alignIO_out)
	sliced_alns = slice_by_name(alignIO_out)
	if commandline_args.shannon_entropy or commandline_args.reflected_shannon:
		#print("Conservation/entropy result:")
		entropyDict = shannon_entropy(alignIO_out, aa_list, commandline_args)
		#for x in entropyDict:
		#	print(x, entropyDict[x])
	aln_index_dict=defaultdict(dict)
	struc_annotation = defaultdict(dict)
	if commandline_args.structure_paths:
		for alngroup_name in sliced_alns:
			print(alngroup_name)
			current_path = [s for s in commandline_args.structure_paths if alngroup_name in s]
			if len(current_path) < 1:
				alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name])
				AlignmentGroup.randomize_gaps(alngroup_name_object, aa_list)
				alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
				for aln_index in alnindex_col_distr:
					struc_annotation[alngroup_name][aln_index] = 'BEHOS'
				for aln_index in alnindex_col_distr:
						aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
			else:
				alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name],current_path[0])
				struc_to_aln_index_mapping=AlignmentGroup.create_aln_struc_mapping(alngroup_name_object)
				AlignmentGroup.randomize_gaps(alngroup_name_object, aa_list)
				alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
				for aln_index in alnindex_col_distr:
						aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
				if commandline_args.secondary_structure:
					struc_annotation[alngroup_name] = AlignmentGroup.ss_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				elif commandline_args.burried_exposed:
					struc_annotation[alngroup_name] = AlignmentGroup.depth_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				elif commandline_args.both:
					struc_annotation[alngroup_name] = AlignmentGroup.both_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
				else:
					parser.print_help()
					raise ValueError("When a structure is defined, one of the matrix options are rquired!")
		alnindex_score = compute_score(aln_index_dict, struc_annotation)
	else:										#Case of no structure defined outputs
		for alngroup_name in sliced_alns:
			alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name])
			AlignmentGroup.randomize_gaps(alngroup_name_object, aa_list)
			alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
			for aln_index in alnindex_col_distr:
				aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
		alnindex_score = compute_score(aln_index_dict,list(sliced_alns.keys()))
	#if commandline_args.shannon_entropy or commandline_args.reflected_shannon:		#temporary for plotting
	#	upsidedown_horizontal_gradient_bar(alnindex_score, list(sliced_alns.keys()))
	
	alnindex = sorted(alnindex_score.keys())
	posdata,negdata = uninterrupted_stretches(alnindex, alnindex_score)

	for x in sorted(posdata.keys()):
		print(x, posdata[x])

if __name__ == '__main__':
	sys.exit(main())

