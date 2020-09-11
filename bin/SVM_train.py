#!/usr/bin/env python3
"""
Generate SVM from alignment segments.
Computes a decision function from csv generated with MultiPhyMeas
"""
#print(__doc__)
import re, os, sys, csv, math, argparse, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import svm
import _pickle as cPickle
from bin.SVM_test import load_csv_data, trim_data_by_top_segments, csv_to_segments_annotated_by_aln, recalculate_data_by_averaging_segments, use_absolute_length_of_segments

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('output_path', help='Output path', type=str)
    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-tp','--penalty', help='Penalty for training algorithm. (Default = 1)', type=float, default=1)
    parser.add_argument('-k','--kernel', help='Kernel for the training algorithm', type=str, 
                                            choices=['linear', 'poly', 'rbf', 'sigmoid', 'precomputed'], default='rbf')
    parser.add_argument('-g','--gamma', help='Gamma function for training algorithm', \
                                            choices=['auto', 'scale'], default='auto')
    parser.add_argument('-l','--length_type_calculation', help='Choose what type of segment calculation should be used.\
        \n\t absolute:   absolute length of the segments.\
        \n\t normalized: length of segments is normalized with the total alignment length.\
        \n\t cms:        average position (center of mass) from all segments per alignment.', choices=['absolute', 'normalized', 'cms'], default='normalized')
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to the top segments that cover\
        \nthis percentage of the total normalized length and weight. (Default = 0.5)', type=float, default=0.5)
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

def plot_decision_function(classifier, X, y, sample_weight, axis, title, aln_names):
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
    scatter = axis.scatter( X[:, 0],
                            X[:, 1],
                            c=y,
                            cmap=plt.cm.viridis)
    
    #Data point labels needs work
    # for i, txt in enumerate(aln_names):
    #     if re.search('A_',txt):
    #         print(txt, X[:, 0][i], X[:, 1][i])
    #         axis.annotate(txt.replace("A_",""), (X[:, 0][i], X[:, 1][i]))

    #Size legend needs work
    # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    # size_labels = list()
    # for label in labels:
    #     size_labels.append(round(math.sqrt(int(re.findall(r'\d+', label)[0]))))
    # legend2 = axis.legend(handles, size_labels, loc="upper right", title="Segment length")
    
    plt.xlim(0, math.ceil(max(X[:, 0])))
    plt.ylim(0, math.ceil(max(X[:, 1])))
    axis.set_title(title)

def train_classifier(X, y, penalty, gamma, kernel, sample_weight=''):
    '''Fits the classifier'''
    decision_function = svm.SVC(C=penalty, gamma=gamma, kernel=kernel)
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
    csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    if comm_args.length_type_calculation == 'cms':
        csv_list = recalculate_data_by_averaging_segments(csv_list)
    if comm_args.length_type_calculation == 'absolute':
        csv_list = use_absolute_length_of_segments(csv_list)
    
    ###   Train the classifier  ###
    X, y, sample_weight, maxX, maxY, aln_names = load_csv_data(csv_list)
    if comm_args.length_type_calculation != 'absolute':
        sample_weight = [math.log(x) for x in sample_weight]
    decision_function = train_classifier(X, y, comm_args.penalty, comm_args.gamma, comm_args.kernel, sample_weight=sample_weight)

    ###   Save the classifier   ###
    with open(comm_args.output_path, 'wb') as classifier_output:
        cPickle.dump(decision_function, classifier_output)

    ###   Save associated max feature values   ###
    max_features = {"maxX":maxX, "maxY": maxY}
    data = [max_features, commandline_arguments]
    with open(str(comm_args.output_path)+".json", 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
    # with open(str(comm_args.output_path)+".maxvals", 'w') as max_features:
    #     csv_writer = csv.writer(max_features, delimiter=',')
    #     csv_writer.writerow([maxX, maxY])
    print("Max on X axis:", maxX, "\nMax on Y axis:", maxY)

    ###   Plot the classifier   ###
    if comm_args.plot_df:
        fig, axis = plt.subplots(1, 1, )
        plot_decision_function(decision_function, X, y, sample_weight, axis, "Decision function", aln_names)
        plt.tight_layout()
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))