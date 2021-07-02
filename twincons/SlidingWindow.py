#!/usr/bin/env python3
"""Calculates conservation score for sliding window of one alignment"""
import sys, argparse, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Bio import AlignIO

import twincons.TwinCons as TwinCons

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Slide two groups of an alignment and calculate a score for each sliding position')
	parser.add_argument('alignment_path', help='Path to folder with alignment files.')
	parser.add_argument('-w','--window', help='Window for sliding the two groups', type=int, default=1)
	parser.add_argument('-bt','--bad_score_threshold', help='Threshold for number of allowed bad scores when calculating length of positive sections.', type=int, default=0, required=False)
	parser.add_argument('-st','--score_threshold', help='Absolute value of threshold for a position to be considered positive or negative', type=int, default=1, required=False)
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

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

def uninterrupted_stretches(alnindex, alnindex_score,comm_args):
	"""Calculates lengths of uninterrupted lengths of positive and negative scores given a threshold;
	Also associates these scores with the last position of the alignment index.
	"""
	posdata,negdata={},{}
	unint_pos_len,unint_neg_len,wei,bad_length=0,0,0,0
	for x,has_more in lookahead(range(0,len(alnindex))):
		if alnindex_score[alnindex[x]] > comm_args.score_threshold:
			wei+=alnindex_score[alnindex[x]]
			unint_pos_len+=1
			if unint_neg_len != 0:
				negdata[alnindex[x-1]]=unint_neg_len
				unint_neg_len = 0
		elif alnindex_score[alnindex[x]] < -comm_args.score_threshold:
			unint_neg_len+=1
			if bad_length >= comm_args.bad_score_threshold:
				if unint_pos_len != 0:
					posdata[alnindex[x-1]]=(unint_pos_len,wei)
					wei,unint_pos_len,bad_length=0,0,0
			else:
				bad_length+=1
		else:				#in case of score threshold different from 0
			if bad_length >= comm_args.bad_score_threshold:
				if unint_pos_len != 0:
					posdata[alnindex[x-1]]=(unint_pos_len,wei)
					unint_pos_len,wei,bad_length = 0,0,0
			else:
				bad_length+=1
			if unint_neg_len != 0:
				negdata[alnindex[x-1]]=unint_neg_len
				unint_neg_len = 0
		if has_more is False:
			if unint_pos_len != 0:
				posdata[alnindex[x]]=(unint_pos_len,wei)
			if unint_neg_len != 0:
				negdata[alnindex[x]]=unint_neg_len
	return posdata, negdata


def main(commandline_args):
	comm_args = create_and_parse_argument_options(commandline_args)
	alignment_file = TwinCons.read_align(comm_args.alignment_path)
	sliced_alignments = TwinCons.slice_by_name(alignment_file)
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
		alnindex_score, gapped_sliced_alns, number_of_aligned_positions, gp_mapping=TwinCons.main(['-as',format(second_aln, "fasta"), '-r', '-mx', 'blosum62'])
		

		out_dict={}
		for x in alnindex_score.keys():
			out_dict[x] = alnindex_score[x][0]
		slided_scores[i] = out_dict
	

	for file in slided_scores:
		print ("Increment is "+str(file))
		alnindex = sorted(slided_scores[file].keys())
		posdata,negdata = uninterrupted_stretches(alnindex, slided_scores[file],comm_args)
		for x in sorted(posdata.keys()):
			print(x, posdata[x])

if __name__ == "__main__":
	main(sys.argv[1:])
