#!/usr/bin/env python3
"""Calculates conservation score for multiple alignments"""
import PhyMeas, SlidingWindow
import re, os, sys, csv, getopt, argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
#from labellines import labelLine, labelLines
from operator import itemgetter

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('alignment_path', help='Path to folder with alignment files.')
	parser.add_argument('output_path', help='Path to image for output.')
	parser.add_argument('-co', '--calculation_options', help='Options for score calculation. See README for details.', required=True, type=str)
	parser.add_argument('-t','--threshold', help='Threshold for number of allowed bad scores when calculating length of positive sections.', type=int, default=1, required=False)
	parser.add_argument('-avew','--average_weight', help='Instead of plotting total weight per segment, plot average weight', action="store_true", required=False)
	parser.add_argument('-it','--intensity_threshold', help='Threshold for intensity over which a score is considered truly positive.', type=float, default=1, required=False)
	parser.add_argument('-s','--structure_path', help='Path to folder with structure files; names should match alignment groups within files.')
	parser.add_argument('-w','--window', help='Window for sliding the two groups', type=int, required=False)
	parser.add_argument('-c','--csv', help='Output length and weight distributions in a csv file. Uses the output file name specified by appending .csv', required=False, action="store_true")
	parser.add_argument('-l','--leg', help='Do not write out a legend', default=False, action="store_true")
	
	#calculation_type = parser.add_mutually_exclusive_group(required=True)
	#calculation_type.add_argument('-bl','--blosum', help='Use BLOSUM62 for calculation of scores', action="store_true")
	#calculation_type.add_argument('-lg','--le_gascuel', help='Use Le & Gascuel matrix for calculation of scores', action="store_true")
	#calculation_type.add_argument('-cons','--conservation', help='Use reverse entropy for calculation of scores.', action="store_true")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def make_length_distr(df,comm_args,group_dict,aln_total_lengths):
	'''Takes in dataframe with values per file and returns a length distribution 
	dictionary with keys files and values a tuple with:
		normalized lengths of uninterrupted positive scoring positions (segments)
		segment weights
		segment lengths
	(Can be interupted by 1 low scoring position (or more set with -t)
	This means that a sequence of ++-+-+-+-- will have a length of 5)
	Also returns a multilevel dictionary with the structure:
	file -> segment length -> array with tuple data of weight, aln position, normalized length
	'''
	weight_distr={}
	length_to_weight={}
	for file in df:
		###   l - alignment position; i - segment length; k - threshold memory; w - segment weight   ###
		l,i,k,w=0,0,0,0
		alignment_length = len(group_dict[file])
		for pos,has_more in PhyMeas.lookahead(df[file]):
			l+=1
			#print(l,file)
			if has_more is False:
				scaled_length = i/aln_total_lengths[file][1]
				#print(i, 'last')
				if i > 0 and w > 1:
					if file not in weight_distr.keys():
						weight_distr[file]=[]
						###
						length_to_weight[file]={}
						length_to_weight[file][i]=[]
						###
					weight_distr[file].append((scaled_length,w,i))
					###
					if i not in length_to_weight[file].keys():
						length_to_weight[file][i]=[]
					length_to_weight[file][i].append((w,l,scaled_length))
					###
			if pos > comm_args.intensity_threshold:
				#print(i, 'greater')
				k=0
				i+=1
				w+=pos
			elif pos <= comm_args.intensity_threshold:
				scaled_length = i/aln_total_lengths[file][1]
				if k == comm_args.threshold:
					#print(i, 'smaller k is threshold')
					if i > 0 and w > 1:
						if file not in weight_distr.keys():
							weight_distr[file]=[]
							###
							length_to_weight[file]={}
							length_to_weight[file][i]=[]
							###
						weight_distr[file].append((scaled_length,w,i))
						###
						if i not in length_to_weight[file].keys():
							length_to_weight[file][i]=[]
						length_to_weight[file][i].append((w,l,scaled_length))
						###
						i,k,w=0,0,0
				elif k < comm_args.threshold:
					#print(i, 'smaller k is not threshold')
					#w+=pos	#Think about these two
					#i+=1	
					k+=1
	return weight_distr, length_to_weight

def scatter_plot(comm_args,weight_distr):
	'''Outputs scatter plot image with colors depending on the number of the input alignments.
	'''
	###   Defining Color scheme   ###
	ax = plt.subplot()
	if len(weight_distr) == 20:
		colors = matplotlib.cm.tab20(np.linspace(0, 1, len(weight_distr)))
	elif len(weight_distr) == 10:
		#In cases of just 10 alignments it assumes two sets of 5 each and uses a seismic/divergent gradient.
		#colors = matplotlib.cm.seismic(np.linspace(0, 1, len(weight_distr)))
		colors = matplotlib.cm.PRGn(np.linspace(0, 1, len(weight_distr)))
	elif len(weight_distr) > 21 and re.match('A_', sorted(weight_distr.keys())[0]):
		#Creates a color mapping for groups of alignments defined with prepended A_, B_, C_ and so on.
		#Used in cases where there are too many alignments to properly discriminate with colors.
		colors=[]
		letters=[]
		for letter in sorted(weight_distr.keys()):
			letters.append(letter.split('_')[0])

		color_set = matplotlib.cm.viridis(np.linspace(0, 1, len(set(letters))))
		letter_color_map = dict(zip(sorted(set(letters)),color_set))

		for letter in letters:
			colors.append(letter_color_map[letter])
	else:
		colors = matplotlib.cm.tab20(np.linspace(0, 1, len(weight_distr)))

	###   Plotting   ###
	sorted_names = sorted(weight_distr.keys())
	file_number = 0
	label_order_tups,label_order = [], []
	for file, color in zip(sorted_names,colors):
		scaled_lengths = [n[0] for n in weight_distr[file]]
		if comm_args.average_weight:
			segment_weights = [n[1]/n[2] for n in weight_distr[file]]
		else:
			segment_weights = [n[1] for n in weight_distr[file]]
		segment_lengths = [n[2] for n in weight_distr[file]]
		abs_length = [n**2 for n in segment_lengths]
		plt.scatter(scaled_lengths, segment_weights, label=re.sub(r'\.fas.*','',file),marker='.',s=abs_length,color=color)
		label_order_tups.append((file_number,max(segment_weights)+max(segment_lengths)))
		file_number += 1
	
	###   Legend labels ordering   ###
	handles, labels = ax.get_legend_handles_labels()
	for tup in sorted(label_order_tups, key = itemgetter(1), reverse=True):
		label_order.append(tup[0])
	ordered_labels = [labels[idx] for idx in label_order]
	ordered_handles = [handles[idx] for idx in label_order]

	###   Legend block   ###
	if not comm_args.leg:
		lgnd = plt.legend(ordered_handles,ordered_labels,bbox_to_anchor=(1.04,1), borderaxespad=0)
		for n in range(len(weight_distr)):
			lgnd.legendHandles[n]._sizes = [30]
	plt.savefig(comm_args.output_path, dpi=600, bbox_inches='tight')
	return True

def csv_output(comm_args, file_to_data):
	'''Writes out data used for generating the plot in a csv file.
	'''
	if comm_args.output_path.endswith('.png'):
		file_for_writing = re.sub(r'.png','.csv',comm_args.output_path)
	else:
		file_for_writing = comm_args.output_path+'.csv'

	with open(file_for_writing, mode ="w") as output_csv:
		csv_writer = csv.writer(output_csv)
		csv_writer.writerow(['File name', 'Segment length', 'Normalized segment length','Total segment weight', 'Alignment position'])
		for file, length_buckets in sorted(file_to_data.items()):
			for length, weights in sorted(length_buckets.items()):
				for single_weight in sorted(weights,key=itemgetter(0)): 
					csv_writer.writerow([file, length, single_weight[2],single_weight[0], single_weight[1]])
	return True


def main(commandline_args):
	comm_args = create_and_parse_argument_options(commandline_args)
	if comm_args.window:
		if not os.path.isfile(comm_args.alignment_path):
			raise ValueError("In case of specified window for sliding (-w argument), the alignment path  must be a single file!")
		SlidingWindow.main([comm_args.alignment_path,'-w '+str(comm_args.window)])
		sys.exit()

	group_dict={}
	aln_lengths={}
	for file in os.listdir(comm_args.alignment_path):
		if re.findall(r'(.*\/)(.*)(\.fasta|\.fas|\.fa)',comm_args.alignment_path+file):
			###   Constructing arguments for PhyMeas   ###
			out_dict={}
			list_for_phymeas = ['-a',comm_args.alignment_path+file, '-r']
			for opt in comm_args.calculation_options.split(" "):
				if re.search(r'_', opt):
					list_for_phymeas.append('-'+opt.split('_')[0])
					list_for_phymeas.append(opt.split('_')[1])
				else:
					list_for_phymeas.append('-'+opt)
			print(list_for_phymeas)
			###   Executing PhyMeas   ###
			alnindex_score,sliced_alns,number_of_aligned_positions=PhyMeas.main(list_for_phymeas)
			###   Alignment lengths and organizing PhyMeas outputs   ###
			for x in sliced_alns:
				aln_lengths[file]=(sliced_alns[x].get_alignment_length(),number_of_aligned_positions)
				break
			for x in alnindex_score.keys():
				out_dict[x] = alnindex_score[x][0]
			group_dict[file] = out_dict
		else:
			raise ValueError("Directory must have only .fa, .fas or .fasta alignment files!")
	
	df = pd.DataFrame.from_dict(group_dict)
	weight_distr, length_to_weight = make_length_distr(df,comm_args,group_dict,aln_lengths)
	scatter_plot(comm_args,weight_distr)
	
	if comm_args.csv:
		csv_output(comm_args, length_to_weight)

if __name__ == "__main__":
	main(sys.argv[1:])