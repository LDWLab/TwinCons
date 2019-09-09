#!/usr/bin/env python3
"""Calculate and visualize conservation between two groups of sequences from one alignment"""
import re, sys, random, Bio.Align, argparse, random, math, matplotlib
matplotlib.use('Agg')
import numpy as np
from datetime import date
from Bio import AlignIO
from io import StringIO
from textwrap import wrap
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from Bio.SeqUtils import IUPACData
from AlignmentGroup import AlignmentGroup
import Sequence_Weight_from_Tree
from MatrixLoad import PAMLmatrix
from Bio.SubsMat import MatrixInfo

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Calculate and visualize conservation between two groups of sequences from one alignment')
	parser.add_argument('-o','--output_path', help='Temp output path for svgs.')
	input_file = parser.add_mutually_exclusive_group(required=True)
	input_file.add_argument('-a','--alignment_path', help='Path to alignment file')
	input_file.add_argument('-as','--alignment_string', help='Alignment string', type=str)
	parser.add_argument('-cg','--cut_gaps', help='Algorithm will cut alignment gaps with more than specified value gaps. Decimal format, eg 0.9', type=float)
	parser.add_argument('-s','--structure_paths', nargs='+', help='Paths to structure files, for score calculation. Does not work with --nucleotide!')
	parser.add_argument('-sy','--structure_pymol', nargs='+', help='Paths to structure files, for plotting a pml.')
	parser.add_argument('-phy','--phylo_split', help='Split the alignment in two groups by constructing a tree instead of looking for _ separated strings.', action="store_true")
	output_type_group = parser.add_mutually_exclusive_group(required=True)
	output_type_group.add_argument('-p', '--plotit', help='Plots the calculated score as a bar graph for each alignment position.', action="store_true")
	output_type_group.add_argument('-pml', '--write_pml_script', help='Writes out a PyMOL coloring script for any structure files that have been defined', action="store_true")
	output_type_group.add_argument('-r', '--return_within', help='To be used from within other python programs. Returns dictionary of alnpos->score.', action="store_true")
	output_type_group.add_argument('-csv', '--return_csv', help='Saves a csv with alignment position -> score.', action="store_true")
	output_type_group.add_argument('-jv', '--jalview_output', help='Saves an annotation file for Jalview.', action="store_true")
	entropy_group = parser.add_mutually_exclusive_group(required=True)
	entropy_group.add_argument('-lg','--leegascuel', help='Use LG matrix for score calculation', action="store_true")
	entropy_group.add_argument('-bl','--blosum', help='Use Blosum62 matrix for score calculation', action="store_true")
	entropy_group.add_argument('-nc','--nucleotide', help='Use nucleotide matrix for score calculation', action="store_true")
	entropy_group.add_argument('-e','--shannon_entropy', help='Use shannon entropy for conservation calculation.', action="store_true")
	entropy_group.add_argument('-c','--reflected_shannon', help='Use shannon entropy for conservation calculation and reflect the result so that a fully random sequence will be scored as 0.', action="store_true")
	structure_option = parser.add_mutually_exclusive_group()
	structure_option.add_argument('-ss','--secondary_structure', help = 'Use substitution matrices derived from data dependent on the secondary structure assignment.', action="store_true")
	structure_option.add_argument('-be','--burried_exposed', help = 'Use substitution matrices derived from data dependent on the solvent accessability of a residue.', action="store_true")
	structure_option.add_argument('-ssbe','--both', help = 'Use substitution matrices derived from data dependent on both the secondary structure and the solvent accessability of a residue.', action="store_true")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def read_align(aln_path):
	'''
	Reads the fasta file and gets the sequences.
	'''
	alignments = AlignIO.read(open(aln_path), "fasta")
	return alignments

def count_aligned_positions(aln_obj):
	'''
	Counts how many positions are aligned (less than 10% gaps)
	'''
	number_seqs = len(aln_obj)
	aligned_positions = 0
	for i in range(0,aln_obj.get_alignment_length()):
		if aln_obj[:,i].count('-')/number_seqs < 0.1:
			aligned_positions+=1
	if aligned_positions == 0:
		raise ValueError("Alignment:\n"+str(aln_obj)+"\nhas no positions with less than 90% gaps!")
	return aligned_positions

def remove_extremely_gapped_regions(align,gap_perc):
	'''
	Removes columns of alignment with more than gap_perc gaps.
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

def uniq_resi_list(aln_obj):
	'''
	Creates list of unique AA or nucleotide residues in the given MSA to be used for frequency iterator.
	Also checks if the alignment has AA letters from the IUPAC extended_protein_letters.
	'''
	default_aa_sequence = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	hash_resi=dict()
	for alignment in aln_obj:
		for resi in alignment.seq:
			if not re.match(r'-|X',resi):
				hash_resi[resi]='null'
	if not all (x in IUPACData.extended_protein_letters for x in hash_resi.keys()):
		raise ValueError("Alignment:\n"+str(aln_obj)+"\nhas AA letters not found in the IUPAC extended list!")
	if len(Counter(hash_resi.keys())) > len(Counter(default_aa_sequence)):
		raise ValueError("Alignment has non-standard AAs:\n"+' ,'.join(hash_resi.keys()))
	return default_aa_sequence


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
	if len(uniq_prot_list) != 2:
		raise ValueError("For now does not support more than two groups! Offending groups are "+str(uniq_prot_list))
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
	return np.array(PAMLmatrix('../matrices/'+struc_anno+'.dat').lodd)

def nucl_matrix():
	'''Return a substitution matrix for nucleotides.
	Should be merged with blos_matrix() intoone general matrix creation method
	'''
	nucl_sequence = ['A','U','G','C']
	nuc_mx = np.array([[4,-5,-5,-5],[-5,4,-5,-5],[-5,-5,4,-5],[-5,-5,-5,4],])
	testvr = np.repeat(1/len(nucl_sequence),len(nucl_sequence))
	baseline = float(testvr@np.array(nuc_mx)@testvr.T)
	revtestA=np.add(np.array(nuc_mx), abs(baseline))
	if int(testvr@revtestA@testvr.T) != 0:
		raise ValueError("Wasn't able to baseline the substitution matrix correctly!")
	return np.add(np.array(nuc_mx),abs(baseline))

def blos_matrix():
	'''Baseline and return a numpy form of the BLOSUM62 matrix.
	Need to make it more general (eg use BLOSUM50 and so on)
	'''
	aa_sequence = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	loddmx = []
	for aa in aa_sequence:
		linemx=[]
		for aa2 in aa_sequence:

			if (aa,aa2) in MatrixInfo.blosum62:
				linemx.append(MatrixInfo.blosum62[(aa,aa2)])
				#print(aa,aa2,MatrixInfo.blosum62[(aa,aa2)])
			else:
				linemx.append(MatrixInfo.blosum62[(aa2,aa)])
				#print(aa2,aa,MatrixInfo.blosum62[(aa2,aa)])
		loddmx.append(linemx)
	testvr = np.repeat(1/len(aa_sequence),len(aa_sequence))
	baseline = float(testvr@np.array(loddmx)@testvr.T)
	revtestA=np.add(np.array(loddmx), abs(baseline))
	if int(testvr@revtestA@testvr.T) != 0:
		raise ValueError("Wasn't able to baseline the substitution matrix correctly!")
	return np.add(np.array(loddmx),abs(baseline))

def compute_score(commandline_args,aln_index_dict, *args):
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
					lgmx = np.array(PAMLmatrix('../matrices/LG.dat').lodd)
					alnindex_score[aln_index] = vr1@lgmx@vr2.T
					#print(aln_index,vr1@lgmx@vr2.T)
			else:
				lgmx = np.array(PAMLmatrix('../matrices/LG.dat').lodd)
				alnindex_score[aln_index] = vr1@lgmx@vr2.T
				#print(aln_index,vr1@lgmx@vr2.T)
	elif commandline_args.leegascuel:							#Case of no structure defined inputs; improve this block...
		groupnames = args[0]
		for aln_index in aln_index_dict:
			vr1 = np.array(aln_index_dict[aln_index][groupnames[0]])
			vr2 = np.array(aln_index_dict[aln_index][groupnames[1]])
			lgmx = np.array(PAMLmatrix('../matrices/LG.dat').lodd)
			alnindex_score[aln_index] = vr1@lgmx@vr2.T
			#print(aln_index,vr1@lgmx@vr2.T)
	elif commandline_args.blosum:
		groupnames = args[0]
		for aln_index in aln_index_dict:
			vr1 = np.array(aln_index_dict[aln_index][groupnames[0]])
			vr2 = np.array(aln_index_dict[aln_index][groupnames[1]])
			alnindex_score[aln_index] = vr1@blos_matrix()@vr2.T
	elif commandline_args.nucleotide:
		groupnames = args[0]
		for aln_index in aln_index_dict:
			vr1 = np.array(aln_index_dict[aln_index][groupnames[0]])
			vr2 = np.array(aln_index_dict[aln_index][groupnames[1]])
			alnindex_score[aln_index] = vr1@nucl_matrix()@vr2.T
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
	"""Used for comparison of blscore and entropy/conservation
	"""
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

def upsidedown_horizontal_gradient_bar(out_dict,group_names,comm_args):
	plt.rcParams['image.composite_image'] = False					#So that bars are separate images
	plt.rcParams["figure.figsize"] = (8,8)
	fig, ax = plt.subplots()
	data=[]
	stdevdata=[]
	for x in sorted(out_dict.keys()):
		data.append(out_dict[x][0])
		stdevdata.append(out_dict[x][1])
	bar = ax.bar(range(1,len(data)+1),data, yerr=stdevdata,error_kw=dict(ecolor='gray', lw=0.25))

	def gradientbars(bars,positivegradient,negativegradient):
		ax = bars[0].axes
		lim = ax.get_xlim()+ax.get_ylim()
		for bar in bars:
			bar.set_zorder(1)
			bar.set_facecolor("none")
			x,y = bar.get_xy()
			w, h = bar.get_width(), bar.get_height()
			if h >= 0:
				grad = np.atleast_2d(np.linspace(0,h/max(data),256)).T
				ax.imshow(grad, extent=[x,x+w,y,y+h], cmap=plt.get_cmap(positivegradient), aspect="auto", norm=matplotlib.colors.NoNorm(vmin=0,vmax=1))
			else:			#different gradient for negative values
				grad = np.atleast_2d(np.linspace(0,h/min(data),256)).T
				ax.imshow(grad, extent=[x,x+w,y,y+h], cmap=plt.get_cmap(negativegradient), aspect="auto", norm=matplotlib.colors.NoNorm(vmin=0,vmax=1))
		#ax.set_facecolor('Gray')
		ax.axis(lim)
	if comm_args.leegascuel or comm_args.blosum:
		pamlarray = np.array(blos_matrix())
		#plt.yticks(np.arange(int(np.min(pamlarray)),int(np.max(pamlarray)+2), step=1))
		#plt.plot((0, len(data)+1), (1, 1), 'k-', linewidth=0.5)				#Horizontal line
		#gradientbars(bar,'Blues','Reds')
		gradientbars(bar,'Greens','Purples')
	#In case of no negative values BUG!
	elif comm_args.reflected_shannon or comm_args.shannon_entropy:
		plt.yticks(np.arange(0,4.2, step=0.5))
		gradientbars(bar,'viridis','binary')
	dpi_scaling = 3*len(out_dict)
	plt.savefig('./outputs/'+'-'.join(sorted(group_names))+'.svg',format = 'svg',dpi=dpi_scaling)
	#plt.savefig(comm_args.output_path,format = 'png',dpi=dpi_scaling)
	return True

def data_to_diverging_gradients(datapoint, datalist, positivegradient, negativegradient):
	'''Maps a data point to a diverging colormap depending on whether its above or bellow 0.
	Returns a hex code.
	'''
	if datapoint >= 0:
		grad = np.atleast_2d(np.linspace(0,datapoint/max(datalist),256)).T
		rgb = plt.get_cmap(positivegradient)(grad[len(grad)-1])[0][:3]
	else:
		grad = np.atleast_2d(np.linspace(0,datapoint/min(datalist),256)).T
		rgb = plt.get_cmap(negativegradient)(grad[len(grad)-1])[0][:3]
	return matplotlib.colors.rgb2hex(rgb)

def pymol_script_writer(out_dict,gapped_sliced_alns,comm_args):
	"""Creates the same gradients used for svg output and writes out a .pml file for PyMOL visualization.
	"""
	#Similarly to upsidedown_horizontal_gradient_bar but it doesn't write out figures or plots.
	fig, ax = plt.subplots()
	data = []
	for x in sorted(out_dict.keys()):
		data.append(out_dict[x][0])
	bar = ax.bar(range(1,len(data)+1),data)
	def gradientbars(bars,positivegradient,negativegradient):
		"""Creates a dictionary with keys the alignment index and values a hex code for the color.
		"""
		aln_index_hexcmap = {}
		aln_index=1
		for bar in bars:
			w, h = bar.get_width(), bar.get_height()
			aln_index_hexcmap[aln_index] = data_to_diverging_gradients(h, data, positivegradient, negativegradient)
			aln_index+=1
		return aln_index_hexcmap
	if comm_args.leegascuel or comm_args.blosum or comm_args.nucleotide:
		alnindex_to_hexcolors = gradientbars(bar,'Blues','Reds')
		#alnindex_to_hexcolors = gradientbars(bar,'Greens','Purples')
	elif comm_args.reflected_shannon or comm_args.shannon_entropy:
		alnindex_to_hexcolors = gradientbars(bar,'viridis','binary')

	group_names = list(gapped_sliced_alns.keys())
	#Open .pml file for structure coloring
	pml_output = open(comm_args.output_path,"w")
	#pml_output = open("./"+'-'.join(sorted(group_names))+'.pml',"w")
	pml_output.write("set hash_max, 500\n\
		set cartoon_loop_radius,0.4\n\
		set cartoon_tube_radius,1\n\
		set cartoon_ladder_radius,0.6\n\
		set cartoon_oval_length,1.4\n\
		set cartoon_oval_width,0.6\n\
		set ray_opaque_background, off\n\
		bg_color black\n\
		set ray_trace_mode,1\n\
		set ray_shadows,0\n")

	#Bellow here needs fixing to properly do structures for plotting
	for alngroup_name in group_names:
		#Match groupnames with structure files
		current_path = [s for s in comm_args.structure_pymol if alngroup_name in s]
		
		if len(current_path) < 1:
			pass		#Gotta fix this
			#raise ValueError("Cannot write PyMOL coloring script without at least single matching structure \
			#	and sequence!\nSequence:\t"+alngroup_name+"\nStructure:\t"+str(current_path))
		else:
			#We have to recalculate the structure to alignment mapping
			alngroup_name_object = AlignmentGroup(gapped_sliced_alns[alngroup_name],current_path[0])
			struc_to_aln_index_mapping=AlignmentGroup.create_aln_struc_mapping(alngroup_name_object)
			#Open the structure file
			pml_output.write("load "+current_path[0]+", "+alngroup_name+"\n")
			#For each alignment position, color the appropriate residue with the hex transformed color from the gradient
			for aln_index in alnindex_to_hexcolors.keys():
				if aln_index in struc_to_aln_index_mapping:
					hexcolors_appropriate_for_pml = alnindex_to_hexcolors[aln_index].replace('#','0x')
					pml_output.write("color "+hexcolors_appropriate_for_pml+", "+alngroup_name+" and resi "+str(struc_to_aln_index_mapping[aln_index])+"\n")
	pml_output.write("super "+group_names[0]+","+group_names[1]+"\n")
	return True

def jalview_output(output_dict, comm_args):
	'''Writes out a Jalview alignment annotation file (http://www.jalview.org/help/html/features/annotationsFormat.html)
	'''
	out_data=list()
	for x in sorted(output_dict.keys()):
		out_data.append(output_dict[x][0])
	
	jv_output = open(comm_args.output_path,"w")
	jv_output.write('JALVIEW_ANNOTATION\n')
	jv_output.write('# Created: '+str(date.today())+"\n")
	jv_output.write('# Contact: ppenev@gatech.edu\n')
	jv_output.write('BAR_GRAPH	TWINCONS	')
	for position in sorted(output_dict.keys(), key=abs):
		color_hex = data_to_diverging_gradients(output_dict[position][0], out_data, 'Greens', 'Purples')
		jv_output.write(str(output_dict[position][0])+'['+str(color_hex).replace('#','')+']|')
	return True

def decision_maker(commandline_args,alignIO_out_gapped,deepestanc_to_child,aa_list):
	"""Checks through the commandline options and does the appropriate calculations; gap randomizations.
	Returns a dictionary of alignment position -> computed score.
	"""
	alignIO_out=alignIO_out_gapped[:,:]
	if commandline_args.phylo_split:
		sliced_alns = Sequence_Weight_from_Tree.slice_by_anc(alignIO_out, deepestanc_to_child)
	else:
		sliced_alns = slice_by_name(alignIO_out)

	if commandline_args.shannon_entropy or commandline_args.reflected_shannon:
		return shannon_entropy(alignIO_out, aa_list, commandline_args)
	
	aln_index_dict=defaultdict(dict)
	struc_annotation = defaultdict(dict)
	if commandline_args.structure_paths:
		for alngroup_name in sliced_alns:
			current_path = [s for s in commandline_args.structure_paths if alngroup_name in s]
			if len(current_path) < 1:					#FIX In case number of structures is fewer than the number of alignment groups
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
					raise ValueError("When a structure is defined, one of the matrix options are required!")
		return compute_score(commandline_args,aln_index_dict, struc_annotation)
	elif commandline_args.leegascuel or commandline_args.blosum or commandline_args.nucleotide:							#Case of no structure defined outputs
		for alngroup_name in sliced_alns:
			alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name])
			AlignmentGroup.randomize_gaps(alngroup_name_object, aa_list)
			alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,aa_list,len(alignIO_out[0]))
			for aln_index in alnindex_col_distr:
				aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
		return compute_score(commandline_args,aln_index_dict,list(sliced_alns.keys()))

def main(commandline_arguments):
	'''Main entry point'''
	comm_args = create_and_parse_argument_options(commandline_arguments)
	if comm_args.alignment_path:
		alignIO_out_gapped=read_align(comm_args.alignment_path)
	elif comm_args.alignment_string:
		alignIO_out_gapped = list(AlignIO.parse(StringIO(comm_args.alignment_string), 'fasta'))[0]
	randindex_norm=defaultdict(dict)
	number_of_aligned_positions = count_aligned_positions(alignIO_out_gapped)
	if comm_args.cut_gaps:
		tempaln = alignIO_out_gapped[:,:]
		alignIO_out_gapped = Bio.Align.MultipleSeqAlignment([])
		gp_mapping,alignIO_out_gapped,alen=remove_extremely_gapped_regions(tempaln,float(comm_args.cut_gaps))
	
	if comm_args.phylo_split:
		tree = Sequence_Weight_from_Tree.tree_contruct(alignIO_out_gapped)
		deepestanc_to_child = Sequence_Weight_from_Tree.find_deepest_ancestors(tree)
		gapped_sliced_alns = Sequence_Weight_from_Tree.slice_by_anc(alignIO_out_gapped, deepestanc_to_child)
	else:
		deepestanc_to_child = {}
		gapped_sliced_alns = slice_by_name(alignIO_out_gapped)

	for rand_index in range(0,10):
		"""Every calculation and gap filling of alignment is performed 10 times.
		This is done to dampen the errors in heavily gapped regions.
		Allows us to add errorbars on the output graph."""
		if comm_args.nucleotide:
			#Add function to Change alignIO_out_gapped so that T -> U
			randindex_norm[rand_index] = decision_maker(comm_args,alignIO_out_gapped,deepestanc_to_child,['A','U','G','C'])
		else:
			randindex_norm[rand_index] = decision_maker(comm_args,alignIO_out_gapped,deepestanc_to_child,uniq_resi_list(alignIO_out_gapped))
	#Calculating mean and stdev per alignment position
	position_defined_scores=defaultdict(dict)
	for x in randindex_norm.keys():
		for pos in randindex_norm[x]:
			if pos not in position_defined_scores:
				position_defined_scores[pos]=[]
			position_defined_scores[pos].append(randindex_norm[x][pos])

	output_dict = {}
	for x in position_defined_scores.keys():				#If standard deviation is too big, set the result as 0
		#print(x, position_defined_scores[x], abs(np.average(position_defined_scores[x])), np.std(position_defined_scores[x]))
		if 2*np.std(position_defined_scores[x]) > abs(np.average(position_defined_scores[x]))/1.5:
			output_dict[x] = (0, np.std(position_defined_scores[x]))
		else:
			output_dict[x] = (np.average(position_defined_scores[x]), np.std(position_defined_scores[x]))
	
	if comm_args.plotit:									#for plotting
		upsidedown_horizontal_gradient_bar(output_dict, list(gapped_sliced_alns.keys()),comm_args)
	elif comm_args.write_pml_script:
		pymol_script_writer(output_dict, gapped_sliced_alns,comm_args)
	elif comm_args.return_within:
		return output_dict, gapped_sliced_alns, number_of_aligned_positions
	elif comm_args.return_csv:
		for x in sorted(output_dict.keys(), key=abs):
			print(str(x)+','+str(output_dict[int(x)][0])+','+str(output_dict[int(x)][1]))
	elif comm_args.jalview_output:
		jalview_output(output_dict, comm_args)

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))