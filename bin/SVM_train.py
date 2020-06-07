#!/usr/bin/env python3
"""
Generate SVM from alignment segments.
Computes a decision function from csv generated with MultiPhyMeas
"""
#print(__doc__)
import re, sys, csv, math, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import svm
import _pickle as cPickle
from SVM_test import load_csv_data, trim_data_by_top_segments

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('output_path', help='Output path', type=str)
    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-tp','--penalty', help='Penalty for training algorithm', type=float, default=10)
    parser.add_argument('-g','--gamma', help='Gamma function for training algorithm', \
                                            choices=['auto', 'scale'], default='auto')
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to longest N segments.', type=int)
    parser.add_argument('-cms','--center_mass_segments', help='Use the average position (Center of mass) \
        \nfrom all segments per alignment. Uses total segment length instead of normalized.', action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args, parser

def csv_iterator(csv_location):
    '''Put csv in list'''
    output_list = list()
    first_line = True
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if first_line:
                first_line = False
                continue
            output_list.append(row)
    return output_list

def csv_to_segments_annotated_by_aln(csv_list):
    aln_id_to_data = dict()
    for row in csv_list:
        if row[0] not in aln_id_to_data.keys():
            aln_id_to_data[row[0]] = list()
        aln_id_to_data[row[0]].append((float(row[1]), float(row[2]), float(row[3]), float(row[4])))
    return aln_id_to_data

def average(xs):
    N = float(len(xs))
    return tuple(sum(col)/N for col in zip(*xs))

def recalculate_data_by_averaging_segments(csv_list):
    '''Calculates an average position for segment in the xy.
    Uses total segment lengths, instead of normalized'''
    aln_id_to_data = csv_to_segments_annotated_by_aln(csv_list)
    output_list = list()
    for aln_id, data_items in aln_id_to_data.items():
        averaged_results = average(data_items)
        output_list.append([aln_id, averaged_results[1]*100, averaged_results[0], averaged_results[2], 'NA'])
    return output_list

def plot_decision_function(classifier, X, y, sample_weight, axis, title):
    '''
    Plots the data and the decision function. Besides the classifier function, takes in the sample weights
    and plots it as a size of the datapoints. If they are different than 1.
    '''
    ###   Plots the decision function   ###
    xx, yy = np.meshgrid(np.linspace(0, math.ceil(max(X[:, 0])), 100), np.linspace(0, math.ceil(max(X[:, 1])), 100))
    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    vir_cmap = plt.cm.get_cmap('viridis')

    ###   Draws the decision function as a red line   ###
    axis.contour(xx, yy, Z, levels=[0],colors='r', linestyles=['-'], linewidths=0.5)

    abs_length = [float(n)**2 for n in sample_weight]

    ###   Plot scatter of segments   ###
    scatter = axis.scatter(X[:, 0], X[:, 1], c=y, alpha=0.5, cmap=plt.cm.bone, edgecolors='black', s=abs_length)
    #Size legend needs work
    # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    # size_labels = list()
    # for label in labels:
    #     size_labels.append(round(math.sqrt(int(re.findall(r'\d+', label)[0]))))
    # legend2 = axis.legend(handles, size_labels, loc="upper right", title="Segment length")
    plt.xlim(0, math.ceil(max(X[:, 0])))
    plt.ylim(0, math.ceil(max(X[:, 1])))
    axis.set_title(title)

def train_classifier(X, y, penalty, gamma, sample_weight=''):
    '''Fits the classifier'''
    decision_function = svm.SVC(C=penalty, gamma=gamma)
    if sample_weight != '':
        decision_function.fit(X, y, sample_weight=sample_weight)
    else:
        decision_function.fit(X, y)
    return decision_function

def main(commandline_arguments):
    '''Main entry point'''
    comm_args, parser = create_and_parse_argument_options(commandline_arguments)

    ###   Load alignment segment data   ###
    csv_list = csv_iterator(comm_args.csv_path)
    if comm_args.top_segments:
        csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    if comm_args.center_mass_segments:
        csv_list = recalculate_data_by_averaging_segments(csv_list)

    ###   Train the classifier  ###
    X, y, sample_weight, maxX, maxY, aln_names = load_csv_data(csv_list)
    decision_function = train_classifier(X, y, comm_args.penalty, comm_args.gamma, sample_weight=sample_weight)

    ###   Save the classifier   ###
    with open(comm_args.output_path, 'wb') as classifier_output:
        cPickle.dump(decision_function, classifier_output)

    ###   Save associated max feature values   ###
    with open(str(comm_args.output_path)+".maxvals", 'w') as max_features:
        csv_writer = csv.writer(max_features, delimiter=',')
        csv_writer.writerow([maxX, maxY])
    print("Max on X axis:", maxX, "\nMax on Y axis:", maxY)

    ###   Plot the classifier   ###
    if comm_args.plot_df:
        fig, axis = plt.subplots(1, 1, )
        plot_decision_function(decision_function, X, y, sample_weight, axis, "Decision function")
        plt.tight_layout()
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))