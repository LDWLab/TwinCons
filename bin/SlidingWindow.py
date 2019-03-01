#!/usr/bin/env python3
"""Calculates conservation score for sliding window of one alignment"""
import re, sys, argparse, PhyMeas
import numpy as np
from Bio import AlignIO

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Slide two groups of an alignment and calculate a score for each sliding position')
	parser.add_argument('alignment_path', help='Path to folder with alignment files.')
	parser.add_argument('-w','--window', help='Window for sliding the two groups', type=int, default=1)
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def main(commandline_args):
	comm_args = create_and_parse_argument_options(commandline_args)
	alignment_file = PhyMeas.read_align(comm_args.alignment_path)
	uniq_aa = PhyMeas.uniq_AA_list(alignment_file)
	sliced_alignments = PhyMeas.slice_by_name(alignment_file)
	first_aln = sorted(list(sliced_alignments.keys()))[0]
	slided_scores={}				#Sliding increment -> (scores,alignment objects)
	for i in range(comm_args.window,sliced_alignments[first_aln].get_alignment_length(),comm_args.window):
		second_aln = AlignIO.MultipleSeqAlignment([])
		for record in sliced_alignments[sorted(list(sliced_alignments.keys()))[1]]:
			second_aln.append(record)
		newaln=sliced_alignments[first_aln][:,-(sliced_alignments[first_aln].get_alignment_length()-i):]+sliced_alignments[first_aln][:,:i]
		for record in newaln:
			second_aln.append(record)
		alnindex_score,sliced_alns=PhyMeas.main(['-as',second_aln.format("fasta"), '-r', '-bl'])
		slided_scores[i] = (alnindex_score,sliced_alns)

	print(slided_scores)


if __name__ == "__main__":
	main(sys.argv[1:])
