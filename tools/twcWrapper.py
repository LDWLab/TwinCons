#!/usr/bin/env python3
'''Calculate and visualize conservation between two groups of sequences from one alignment'''
import re, os, sys, argparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from itertools import groupby
import pickle as cPickle
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier,
                              AdaBoostClassifier)



from twincons import TwinCons
import twcCalculateSegments, twcSVMtest, twcTreesTest

MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')
classifiers = {
    'SVM-BBS-CumulativeW9-LG': os.path.join(MODEL_DIR, "BBS_lg_bgfreq_cg0p9__cmsW7_nn__ts0p5_normalized.pkl"),
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
    train_group.add_argument('-c','--classifier', help='Choose a decision boundary classifier.', choices=classifiers.keys())
    train_group.add_argument('-cc','--custom_classifier', help='Path to a custom pickled classifier.', type=str)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)

    if comm_args.classifier:
        classifier_path = classifiers[comm_args.classifier]
    else:
        classifier_path = comm_args.custom_classifier
        if not os.path.isfile(classifier_path):
            raise IOError(f"Could not find the custom decision boundary pickle file at {comm_args.custom_classifier}")
        features_path = classifier_path.replace('.pkl', '') + ".json"
        if not os.path.isfile(features_path):
            raise IOError(f"Could not find the custom decision boundary json file at {features_path}")

    calc_args, minmax_features = twcSVMtest.read_features(classifier_path.replace('.pkl','')+".json")

    twincons_args, calcSegments_args = None, None
    g_list=[list(g) for k,g in groupby(calc_args , lambda i : '-twca' in i or '-csa' in i)]
    for i, args in enumerate(g_list[1:]):
        if args[0] == '-twca':
            twincons_args = g_list[i+2]
        if args[0] == '-csa':
            calcSegments_args = g_list[i+2]
    if twincons_args is None:
        raise IOError(f"The classifier metadata {classifier_path.replace('.pkl', '')}.json does not record the TwinCons "
                      "arguments (-twca) it was trained with, so alignments cannot be scored consistently with it.")
    int_thr, length_thr, positive_as_negative, cmsWindow = 1, 3, False, None
    if calcSegments_args is not None:
        for segm_opt in calcSegments_args:
            if re.match('it_|intensity_threshold_', segm_opt):
                int_thr = int(segm_opt.split('_')[1])
            if re.match('lt_|length_threshold_', segm_opt):
                length_thr = int(segm_opt.split('_')[1])
            if re.match('np|treat_highly_negative_as_conserved', segm_opt):
                positive_as_negative = True
            if re.match('cms|cumulative_segments', segm_opt):
                cmsWindow = int(segm_opt.split('cmsW')[1])

    twincons_alns = ["-r"]
    if comm_args.alignment_paths:
        twincons_alns.append("-a")
        for path in comm_args.alignment_paths:
            if not os.path.isfile(path):
                raise IOError(f"Could not find alignment file at {comm_args.alignment_paths}")
            twincons_alns.append(path)
    else:
        twincons_alns.extend(["-as", comm_args.alignment_string])
    
    twc_args_list = twcCalculateSegments.parse_arguments_for_twc(twincons_args, twincons_alns)
    alnindex_score, sliced_alns, num_alned_pos, gap_mapping = TwinCons.main(twc_args_list)
    segment_stats = twcCalculateSegments.calc_segments_for_aln(alnindex_score, num_alned_pos, 
                                                                int_thr=int_thr, 
                                                                length_thr=length_thr, 
                                                                highly_neg_as_pos=positive_as_negative,
                                                                cmsWindow=cmsWindow)
    
    classifier = cPickle.load(open(classifier_path, 'rb'))
    if isinstance(classifier, SVC):
        pass
    elif isinstance(classifier, DecisionTreeClassifier):
        pass
        #Decision tree
    elif isinstance(classifier, RandomForestClassifier):
        pass
    elif isinstance(classifier, ExtraTreesClassifier):
        pass
    elif isinstance(classifier, AdaBoostClassifier):
        pass
    else:
        raise IOError("Unsuported classifier type!")


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))