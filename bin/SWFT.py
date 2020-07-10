#!/usr/bin/env python3
"""Calculate importance of sequences based on phylogenetic tree."""

import sys, argparse
from collections import Counter
from numpy.random import choice
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from AlignmentGroup import AlignmentGroup
from TwinCons import slice_by_name, read_align, uniq_resi_list
from Sequence_Weight_from_Tree import tree_construct, pairwise_dist, slice_by_anc, find_deepest_ancestors, calculate_weight_vector
from Bio.Phylo.TreeConstruction import DistanceCalculator


def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', help='Path to alignment file')
    parser.add_argument('-nc','--nucleotide_alignments', help='Alignments are nucleotidic', action="store_true")
    parser.add_argument('-w','--weight_by_tree', help='Weight sequences by tree topology.', type=str, 
                                                        choices=['pairwise', 'voronoi'], default='dist')
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def generate_sequence_sampled_from_alignment(aln_obj):
    outseq = ''
    i = 0
    while i < len(aln_obj[0]):
        aa_list = list(set(aln_obj[:, i]))
        distribution = Counter(aln_obj[:, i])
        choice_distr = []
        for aa in aa_list:
            choice_distr.append(distribution[aa]/len(aln_obj[:, i]))
        outseq += choice(aa_list)
        i+=1
    return outseq

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    aln_name = str(comm_args.alignment_file).split("/")[-1].replace(".fasta", "").replace(".fas", "")

    alignIO_out = read_align(comm_args.alignment_file)
    
    sliced_alns = slice_by_name(alignIO_out)

    if comm_args.weight_by_tree == "voronoi":
        calculate_weight_vector(sliced_alns["bL27"], algorithm='voronoi')

    if comm_args.weight_by_tree == "pairwise":
        calculate_weight_vector(sliced_alns["bL27"], algorithm='pairwise')
        # if len(sliced_alns.keys()) != 2:
        #     print("Using UPGMA tree to assign groups!")
        #     if comm_args.nucleotide_alignments:
        #         for sequence in alignIO_out:
        #             sequence.seq = sequence.seq.back_transcribe()
        #         tree = tree_construct(alignIO_out, nucl=True)
        #     else:
        #         tree = tree_construct(alignIO_out)
        #     deepestanc_to_child = find_deepest_ancestors(tree)
        #     sliced_alns = slice_by_anc(alignIO_out, deepestanc_to_child)

        # tree_dict = dict()
        # for alngroup in sliced_alns:
        #     uniq_resi_list(sliced_alns[alngroup]) #Only to check for non-standard
        #     test_seq = generate_sequence_sampled_from_alignment(sliced_alns[alngroup])
        #     test_entry = SeqRecord(Seq(test_seq), id=alngroup+"_testSEQ")
        #     sliced_alns[alngroup].append(test_entry)
        #     tree_dict[alngroup] = tree_construct(sliced_alns[alngroup])

        # if comm_args.weight_by_tree:
        #     seq_names = list()
        #     group_names = list()
        #     for group in sliced_alns:
        #         group_names.append(group)
        #         seq_names.append([x.id for x in sliced_alns[group]])
        #     named_intra1, intra1 = pairwise_dist(tree_dict[group_names[0]], seq_names[0])
        #     named_intra2, intra2 = pairwise_dist(tree_dict[group_names[1]], seq_names[1])

        # calculator = DistanceCalculator('blosum62')
        # aln_obj = sliced_alns[group_names[0]]
        # for seq_obj in aln_obj:
        #     fastd = calculator._pairwise(seq_obj.seq, aln_obj[len(sliced_alns[group_names[0]])-1].seq)
        #     print(seq_obj.id, group_names[0]+"_testSEQ", fastd)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))