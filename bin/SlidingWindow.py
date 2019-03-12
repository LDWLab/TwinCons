#!/usr/bin/env python3
"""Calculates conservation score for sliding window of one alignment"""
import re, sys, argparse
import numpy as np
from Bio import AlignIO
import pandas as pd

import PhyMeas

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Slide two groups of an alignment and calculate a score for each sliding position')
	parser.add_argument('alignment_path', help='Path to folder with alignment files.')
	parser.add_argument('-w','--window', help='Window for sliding the two groups', type=int, default=1)
	parser.add_argument('-t','--threshold', help='Threshold for number of allowed bad scores when calculating length of positive sections.', type=int, default=1, required=False)
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def uninterrupted_stretches(alnindex, alnindex_score):
	"""Calculates lengths of uninterrupted lengths of positive and negative scores;
	Also associates these scores with the last position of the alignment index.
	For now uses > 1 as positive and < -1 as negative. Can be a variable or changed as appropriate.
	"""
	posdata={}
	negdata={}
	pos,neg,wei=0,0,0
	for x,has_more in PhyMeas.lookahead(range(0,len(alnindex))):
		if alnindex_score[alnindex[x]] > 0.5:
			#print('pos',alnindex[x], alnindex_score[alnindex[x]][0])
			wei+=alnindex_score[alnindex[x]]
			pos+=1
			if neg != 0:
				negdata[alnindex[x-1]]=neg
				neg = 0
		elif alnindex_score[alnindex[x]] < -0.5:
			#print('neg',alnindex[x], alnindex_score[alnindex[x]][0])
			neg+=1
			if pos != 0:
				posdata[alnindex[x-1]]=(pos,wei)
				wei=0
				pos = 0
		else:				#in case of using some range between positive and negative scores for random
			#print('rand',alnindex[x], alnindex_score[alnindex[x]][0])
			if pos != 0:
				posdata[alnindex[x-1]]=(pos,wei)
				wei=0
				pos = 0
			if neg != 0:
				negdata[alnindex[x-1]]=neg
				neg = 0
		if has_more is False:
			if pos != 0:
				posdata[alnindex[x]]=(pos,wei)
			if neg != 0:
				negdata[alnindex[x]]=neg
	return posdata, negdata


def main(commandline_args):
	comm_args = create_and_parse_argument_options(commandline_args)
	alignment_file = PhyMeas.read_align(comm_args.alignment_path)
	uniq_aa = PhyMeas.uniq_AA_list(alignment_file)
	sliced_alignments = PhyMeas.slice_by_name(alignment_file)
	first_aln = sorted(list(sliced_alignments.keys()))[0]
	slided_scores={}				#Sliding increment -> (scores,alignment objects)
	for i in range(0,sliced_alignments[first_aln].get_alignment_length(),comm_args.window):
		#print(i)
		second_aln = AlignIO.MultipleSeqAlignment([])
		for record in sliced_alignments[sorted(list(sliced_alignments.keys()))[1]]:
			second_aln.append(record)
		#Reorders an alignment group using the specified window size
		reordered_aln=sliced_alignments[first_aln][:,-(sliced_alignments[first_aln].get_alignment_length()-i):]+sliced_alignments[first_aln][:,:i]
		for record in reordered_aln:
			second_aln.append(record)
		alnindex_score,sliced_alns=PhyMeas.main(['-as',second_aln.format("fasta"), '-r', '-bl'])
		

		out_dict={}
		for x in alnindex_score.keys():
			out_dict[x] = alnindex_score[x][0]
		slided_scores[i] = out_dict
	

	for file in slided_scores:
		print ("Increment is "+str(file))
		alnindex = sorted(slided_scores[file].keys())
		posdata,negdata = uninterrupted_stretches(alnindex, slided_scores[file])
		for x in sorted(posdata.keys()):
			print(x, posdata[x])

if __name__ == "__main__":
	main(sys.argv[1:])
