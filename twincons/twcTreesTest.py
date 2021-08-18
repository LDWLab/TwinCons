#!/usr/bin/env python3
"""
Evaluates alignment entries in csv generated from twcCalculateSegments.
Requires a decision function and max features generated from SVM_train.
Train and test only with the same parameters!
Such parameters can be % cutting gaps, center mass segments, top segments.
"""
import os, re, sys, csv, math, argparse, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import pickle as cPickle

from twincons.twcSVMtest import csv_iterator, \
                                load_csv_data, \
                                read_features, \
                                normalize_features, \
                                bypass_zero_division, \
                                trim_data_by_top_segments, \
                                use_absolute_length_of_segments, \
                                recalculate_data_by_averaging_segments, \
                                flatten_alignment_segment_stats_to_groups
                                


def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('output_path', help='Path to output significant segment results', type=str)
    parser.add_argument('pickle',help='Provide path to classifier pickle binary file. The script will also search for an identically\
                                    \nnamed file with extension ".json" containing comma separated max and min values, for example:\
                                    \npickle file:\trandom_test.pkl\
                                    \nmaximal values:\trandom_test.pkl.json', type=str)
    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to the top segments that cover\
                                    \n this percentage of the total normalized length and weight. Default = 1 (all segments)', 
                                                                                            type=float, default=1)
    parser.add_argument('-l','--length_type_calculation', help='Choose what type of segment calculation should be used.\
        \n\t absolute:   absolute length of the segments.\
        \n\t normalized: length of segments is normalized with the total alignment length.\
        \n\t cms:        average position (center of mass) from all segments per alignment.', choices=['absolute', 'normalized', 'cms'], default='normalized')
    calculate_positive = parser.add_mutually_exclusive_group(required=True)
    calculate_positive.add_argument('-tcp','--test_classifier_precision', help='Provided csv is annotated for testing the classifier.', action="store_true")
    calculate_positive.add_argument('-tqa','--test_query_alignments', help='Provided csv is a query and is not annotated for testing the classifier.', action="store_true")
    parser.add_argument('-pt', '--probability_threshold', nargs=2, metavar=('Significance', 'Step'), 
        help='Probability for identifying significant segments.\
        \n\tSignificance: segment with probability above or equal is considered significant segment. Segment with probability below this number is non-significant.\
        \n\tStep: only used when testing the classifier precision. It makes multiple evaluations of the data moving from 0 to 1 probabilities with the given step.\
        \n\tDefault (0.95, 0.001)', default=[0.95, 0.001], type = float)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def test_function(csv_list, normData, model):
    '''
    Executes prediction and distance calculation on each
    segment from the input for a given decision function.
    '''
    segment_pred_dist = {}
    for test_segment, entry in zip(normData, csv_list):
        segment_pred = model.predict(np.array(test_segment).reshape(1,-1))[0]
        #Predicts probability for both the 0 and 1 classes. Take only prob for positive results.
        segment_prob = model.predict_proba(np.array(test_segment).reshape(1,-1))[0][1]
        if str(entry[0]) not in segment_pred_dist.keys():
            segment_pred_dist[str(entry[0])] = []
        segment_pred_dist[str(entry[0])].append([entry[4],(segment_pred,segment_prob,entry[1])])
    return segment_pred_dist

def bypass_zero_division(x,y):
    if y != 0:
        return x/y
    else:
        return 0

def compare_thr(thr, segments):
    negative, positive = 0, 0
    for i in segments:
        positive += (i[1][1] >= thr)
        negative += (i[1][1] < thr)
    return (negative, positive)

def calc_identity_stats(segment_pred_dist, thr, tp, tn, fp, fn):
    for aln, segments in segment_pred_dist.items():
        if re.match('^A_|^B_', aln) is not None:
            if compare_thr(thr, segments)[1] == 0:
                fn += 1
            if compare_thr(thr, segments)[1] > 0:
                tp += 1
        if re.match('^C_|^D_', aln) is not None:
            if compare_thr(thr, segments)[1] == 0:
                tn += 1
            if compare_thr(thr, segments)[1] > 0:
                fp += 1
    #print(tp, tn, fp, fn)
    return tp, tn, fp, fn

def mass_test(segment_pred_dist, min_threshold=0, max_threshold=2, step=0.1):
    '''For evaluating the decision function with 
    great number of alignments in separate groups.
    '''
    ###   Calculate Sensitivity(tpr), Specificity(tnr) and Precision   ###
    ###   over a range of distances from the decision boundary         ###
    dist_to_stats = {}
    for thr in np.arange(min_threshold,max_threshold,step):
        tp, tn, fp, fn = calc_identity_stats(segment_pred_dist, thr, 0, 0, 0, 0)
        tpr = bypass_zero_division(tp,tp+fn)
        tnr = bypass_zero_division(tn,tn+fp)
        precision = bypass_zero_division(tp,tp+fp)
        dist_to_stats[thr] = (tpr, tnr, precision)
        #print ("Threshold "+str(thr),'\n',"tpr", tpr,"tnr", tnr,'\n',"precision",precision)
    return dist_to_stats

def write_aln_rows(segments, csv_writer, aln, probThr):
    for segment in segments:
        if segment[1][1] >= probThr:
            csv_writer.writerow([aln, segment[1][1], segment[0]])
    return True

def main(commandline_arguments):
    '''Main entry point'''
    comm_args = create_and_parse_argument_options(commandline_arguments)

    ###   Load alignment segment data   ###
    csv_list = csv_iterator(comm_args.csv_path)
    cms = False
    if comm_args.top_segments:
        csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    if comm_args.length_type_calculation == 'cms':
        cms = True
        csv_list = recalculate_data_by_averaging_segments(csv_list)
    if comm_args.length_type_calculation == 'absolute':
        csv_list = use_absolute_length_of_segments(csv_list)

    ###   Test data   ###
    train_args, min_max_features = read_features(comm_args.pickle+".json")
    classifier = cPickle.load(open(comm_args.pickle, 'rb'))
    normalizedData = normalize_features([[float(row[2]),float(row[3])] for row in csv_list], 
                                        min_max_features[0], 
                                        min_max_features[1], 
                                        min_max_features[2], 
                                        min_max_features[3])
    segment_pred_dist = test_function(csv_list, normalizedData, classifier)


    #Fix from here
    if comm_args.test_classifier_precision:
        dist_to_se_sp_pr = mass_test(segment_pred_dist,
            min_threshold=float(comm_args.probability_threshold[0]), 
            max_threshold=float(comm_args.probability_threshold[1])+float(comm_args.probability_threshold[2]), 
            step=float(comm_args.probability_threshold[2]))
        levels_labels_csv = [dist_to_se_sp_pr[thr] for thr in sorted(dist_to_se_sp_pr.keys())]
        roc_stats = [(thr,lev) for lev,thr in zip(levels_labels_csv, sorted(dist_to_se_sp_pr.keys()))]
        with open(comm_args.output_path, mode ="w") as output_csv:
            csv_writer = csv.writer(output_csv)
            csv_writer.writerow(['Alignment name', 'Segment Probability', 'Segment Position in Alignment'])
            for aln, segments in segment_pred_dist.items():
                write_aln_rows(segments, csv_writer, aln)
            csv_writer.writerow(['ROC stats'])
            csv_writer.writerow(['Distance from boundary', 'TPR', 'TNR', 'PREC'])
            for row in roc_stats:
                csv_writer.writerow([row[0],row[1][0],row[1][1],row[1][2]])
    if comm_args.test_query_alignments:
        alnIDwithProb, alnIDtoProb, i = list(), dict(), 0
        if comm_args.length_type_calculation == 'cms':
            grouped_data = flatten_alignment_segment_stats_to_groups(segment_pred_dist)
            for aln in sorted(grouped_data.keys()):
                alnIDtoProb[aln] = grouped_data[aln][0][0]
                alnIDwithProb.append((aln, grouped_data[aln][0][0], i))
                i+=1
        with open(comm_args.output_path, mode ="w") as output_csv:
            csv_writer = csv.writer(output_csv)
            csv_writer.writerow(['Alignment name', 'Segment Probability', 'Segment Position in Alignment'])
            for aln, segments in segment_pred_dist.items():
                segment_id = compare_thr(float(comm_args.probability_threshold[0]), segments)
                alnIDwithProb.append((aln, segment_id[1], i, max(segments, key = itemgetter(1))[1][1] ))
                alnIDtoProb[aln] = segment_id[1]
                i+=1
                write_aln_rows(segments, csv_writer, aln, float(comm_args.probability_threshold[0]))

    ###   Plot the classifier   ###
    if comm_args.plot_df:
        plot_title = re.sub('.csv','',str(re.sub(r'.*/','',comm_args.csv_path))).replace("_", " ")
        X, y, sample_weight, maxX, maxY, minX, minY, aln_names = load_csv_data(csv_list, min_max_features=min_max_features)
        fig, axes = plt.subplots(1, 1, )
        if comm_args.test_classifier_precision:
            plot_decision_function(classifier,X,y, sample_weight, axes, fig,
                                plot_title, aln_names, 
                                decision_levels=dist_to_se_sp_pr, 
                                cms_or_identity=cms)
        elif comm_args.test_query_alignments:
            aln_names_eval = list()
            for alnid in aln_names:
                if cms:
                    aln_names_eval.append(("{:2.1e}".format(alnIDtoProb[alnid])+" "+alnid[:15]))
                else:
                    aln_names_eval.append((str(alnIDtoProb[alnid])+" "+alnid[:15]))
            plot_decision_function(classifier,X,y, sample_weight, axes, fig,
                                plot_title,aln_names_eval,label_order_tups=alnIDwithProb,
                                cms_or_identity=cms,
                                thresholds=comm_args.probability_threshold)
        plt.tight_layout()
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))