#!/usr/bin/env python3
"""Tree testing"""

import re, sys, Bio.Align, argparse, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO
from Bio import Phylo
from AlignmentGroup import AlignmentGroup
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import PhyMeas

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Calculate importance of sequences based on phylogenetic tree.')
	parser.add_argument('alignment_file', help='Path to alignment file')
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def main(commandline_arguments):
	comm_args = create_and_parse_argument_options(commandline_arguments)
	alignIO_out = PhyMeas.read_align(comm_args.alignment_file)
	sliced_alns = PhyMeas.slice_by_name(alignIO_out)


	calculator = DistanceCalculator('blosum62')
	for alngroup_name in sliced_alns:
		dist_mx = calculator.get_distance(sliced_alns[alngroup_name])
		constructor = DistanceTreeConstructor()
		tree = constructor.upgma(dist_mx)
		tree.ladderize()
		treedepths_int = tree.depths(unit_branch_lengths=True)
		treedepths = tree.depths()
		print(alngroup_name)
		print(Phylo.draw_ascii(tree))
		for terminal_leaf in treedepths:
			print(terminal_leaf,treedepths[terminal_leaf],treedepths_int[terminal_leaf])

		
		
		

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))