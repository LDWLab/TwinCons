#!/usr/bin/env python3
'''Calculate and visualize conservation between two groups of sequences from one alignment'''
import re, os, sys, argparse
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from itertools import groupby

from twincons import TwinCons, CalculateSegments, SVM_test

decision_boundaries = {
    'BaliBase-BL62':"data/PKL/BBS_cg09_it1_lt3.pkl"
    }

def create_and_parse_argument_options(argument_list):

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o','--output_path', help='Output path')
    input_file = parser.add_mutually_exclusive_group(required=True)
    input_file.add_argument('-a','--alignment_paths', nargs='+', help='Path to alignment files. If given two files it will use mafft --merge to merge them in single alignment. Fasta format.', action=TwinCons.required_length(1,2))
    input_file.add_argument('-as','--alignment_string', help='Alignment in fasta format as a string.', type=str)
    output_type_group = parser.add_argument_group()
    output_type_group.add_argument('-p', '--plotit', help='Plots the calculated score as a bar graph for each alignment position.', action="store_true")
    output_type_group.add_argument('-pml', '--write_pml_script', help='Writes out a PyMOL coloring script for any structure files that have been defined. Choose between unix or windows style paths for the pymol script.', choices=['unix', 'windows'])
    output_type_group.add_argument('-jv', '--jalview_output', help='Saves an annotation file for Jalview.', action="store_true")
    train_group = parser.add_mutually_exclusive_group(required=True)
    train_group.add_argument('-db','--decision_boundary', help='Choose a decision boundary classifier.', choices=decision_boundaries.keys())
    train_group.add_argument('-cdb','--custom_decision_boundary', help='Path to a custom pickled decision boundary classifier.', type=str)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)

    if comm_args.decision_boundary:
        decision_boundary_path = f"./{decision_boundaries[comm_args.decision_boundary]}"
    else:
        decision_boundary_path = f"./{decision_boundaries[comm_args.custom_decision_boundary]}"
        if not os.path.isfile(decision_boundary_path):
            raise IOError(f"Could not find the custom decision boundary pickle file at {comm_args.custom_decision_boundary}")
        if not os.path.isfile(decision_boundary_path+".json"):
            raise IOError(f"Could not find the custom decision boundary json file at {comm_args.custom_decision_boundary}.json")

    calc_args, max_features = SVM_test.read_features(decision_boundary_path+".json")

    g_list=[list(g) for k,g in groupby(calc_args , lambda i : '-twca' in i or '-csa' in i)]
    for i, args in enumerate(g_list[1:]):
        if args[0] == '-twca':
            twincons_args = g_list[i+2]
        if args[0] == '-csa':
            calcSegments_args = g_list[i+2]
    int_thr, length_thr, positive_as_negative = 1, 3, False
    if calcSegments_args is not None:
        for segm_opt in calcSegments_args:
            if re.match('it_|intensity_threshold_', segm_opt):
                int_thr = int(segm_opt.split('_')[1])
            if re.match('lt_|length_threshold_', segm_opt):
                length_thr = int(segm_opt.split('_')[1])
            if re.match('np|treat_highly_negative_as_conserved', segm_opt):
                positive_as_negative = True

    twincons_alns = ["-r"]
    if comm_args.alignment_paths:
        twincons_alns.append("-a")
        for path in comm_args.alignment_paths:
            if not os.path.isfile(path):
                raise IOError(f"Could not find alignment file at {comm_args.alignment_paths}")
            twincons_alns.append(path)
    else:
        twincons_alns.append("-as", comm_args.alignment_string)
    
    twc_args_list = CalculateSegments.parse_arguments_for_twc(twincons_args, twincons_alns)
    alnindex_score, sliced_alns, num_alned_pos, gap_mapping = TwinCons.main(twc_args_list)
    segment_stats = CalculateSegments.calc_segments_for_aln(alnindex_score, num_alned_pos, 
                                                            int_thr=int_thr, 
                                                            length_thr=length_thr, 
                                                            highly_neg_as_pos=positive_as_negative)
    pass


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))