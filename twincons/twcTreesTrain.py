#!/usr/bin/env python3
"""
Generate ensemble classifier from alignment segments.
Computes a decision function from csv generated with twcCalculateSegments
"""

import os, sys, csv, math, argparse, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier,
                              AdaBoostClassifier)
from sklearn.tree import DecisionTreeClassifier
import pickle as cPickle

from twincons.twcSVMtest import csv_iterator, \
                                load_csv_data, \
                                trim_data_by_top_segments, \
                                recalculate_data_by_averaging_segments, \
                                use_absolute_length_of_segments

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('output_path', help='Output path', type=str)
    parser.add_argument('-mt','--model_type', help='Select the model of the classifier. List multiple for comparison. (Deafult: RandomForest)', 
                        nargs='+', choices=['DecisionTree', 'RandomForest','ExtraTrees', 'AdaBoost'], default=['RandomForest'])
    parser.add_argument('-ne','--n_estimators', help='Parameter for RandomForest, ExtraTrees and AdaBoost. Default = 30', type=int, default=30)
    parser.add_argument('-twca','--twincons_args', help='Arguments used with TwinCons.', nargs='+', type=str)
    parser.add_argument('-csa','--calcsegm_args', help='Arguments used with twcCalculateSegments.', nargs='+', type=str)

    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-l','--length_type_calculation', help='Choose what type of segment calculation should be used.\
        \n\t absolute:   absolute length of the segments.\
        \n\t normalized: length of segments is normalized with the total alignment length.\
        \n\t cms:        average position (center of mass) from all segments per alignment.', choices=['absolute', 'normalized', 'cms'], default='normalized')
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to the top segments that cover\
        \nthis percentage of the total normalized length and weight. (Default = 1 - all data)', type=float, default=1)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def plot_decision_function(model, X, y, model_title, plot_idx, cmap, num_models, aln_names=False, threshold=None, labelOrder=None, weight=None, fullData=False):
    '''
    Plots the data and the decision function. Besides the classifier function, takes in the sample weights
    and plots it as a size of the datapoints. If they are different than 1.
    '''
    axes = plt.subplot(1, num_models, plot_idx, aspect='equal', adjustable='box')
    # Add a title for each subplot
    plt.title(model_title, fontsize=9)
    # Now plot the decision boundary using a fine mesh as input to a
    # filled contour plot
    xx, yy = np.meshgrid(np.linspace(0, math.ceil(max(X[:, 0])), 100),
                         np.linspace(0, math.ceil(max(X[:, 1])), 100))
    # Plot either a single DecisionTreeClassifier or alpha blend the
    # decision surfaces of the ensemble of classifiers
    if isinstance(model, DecisionTreeClassifier):
        Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        cs = plt.contourf(xx, yy, Z, cmap=cmap)
    else:
        # Choose alpha blend level with respect to the number
        # of estimators that are in use
        estimator_alpha = 1.0 / len(model.estimators_)
        for tree in model.estimators_:
            Z = tree.predict(np.c_[xx.ravel(), yy.ravel()])
            Z = Z.reshape(xx.shape)
            cs = plt.contourf(xx, yy, Z, alpha=estimator_alpha, cmap=cmap)

    # Plot the training points, these are clustered together and have a
    # black outline
    if aln_names:
        scatter = plotQuerySegments(X, aln_names, "black", axes, labelOrder, threshold, weight, fullData=fullData)
    else:
        plt.scatter(X[:, 0], X[:, 1], c=y,
                cmap=ListedColormap(['purple', 'green']),
                edgecolor='k', s=20)
    plt.xlim(0)
    return True

def plotQuerySegments(X, aln_names, edgecolor, axis, labelOrder, threshold, sample_weight, fullData=False):
    import seaborn as sns
    from operator import itemgetter
    label_order = list()
    abs_length = [float(n)**1.6 for n in sample_weight]
    scatter = sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=aln_names, 
            palette="tab10", edgecolor=edgecolor, s=abs_length)
    
    ##   Legend labels ordering only if full data  ###
    if fullData:
        handles, labels = axis.get_legend_handles_labels()
        count_bellow_1=0
        for tup in sorted(labelOrder, key = itemgetter(1), reverse=True):
            if tup[1] > threshold:
                count_bellow_1+=1
            label_order.append(tup[2])
        ordered_labels = [labels[idx] for idx in label_order]
        ordered_handles = [handles[idx] for idx in label_order]
        lgnd = plt.legend(ordered_handles[:count_bellow_1],
                            ordered_labels[:count_bellow_1])
    plt.legend(title=f"Alignments with segments\nabove probability {threshold}",
              bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    return scatter

def train_classifier(X, y, model, sample_weight=''):
    '''Fits the classifier'''
    if sample_weight != '':
        model.fit(X, y, sample_weight=sample_weight)
    else:
        model.fit(X, y)
    return model

def main(commandline_arguments):
    '''Main entry point'''
    comm_args = create_and_parse_argument_options(commandline_arguments)
    modelDict = dict(DecisionTree = DecisionTreeClassifier(max_depth=None),
                     RandomForest = RandomForestClassifier(n_estimators=comm_args.n_estimators),
                     ExtraTrees = ExtraTreesClassifier(n_estimators=comm_args.n_estimators),
                     AdaBoost = AdaBoostClassifier(DecisionTreeClassifier(max_depth=7), n_estimators=comm_args.n_estimators))
    ###   Load alignment segment data   ###
    csv_list = csv_iterator(comm_args.csv_path)
    csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    if comm_args.length_type_calculation == 'cms':
        csv_list = recalculate_data_by_averaging_segments(csv_list)
    if comm_args.length_type_calculation == 'absolute':
        csv_list = use_absolute_length_of_segments(csv_list)
    
    ### Normalize data, create features array (X), and identity array (y) ###
    X, y, sample_weight, maxX, maxY, minX, minY, aln_names = load_csv_data(csv_list)
    if comm_args.length_type_calculation != 'absolute':
        sample_weight = [math.log(x) for x in sample_weight]
    
    modelIdx = 1
    for modelName in comm_args.model_type:
        model = modelDict[modelName]
        ###   Train the classifier  ###
        decision_function = train_classifier(X, y, model, sample_weight=sample_weight)
        scores = decision_function.score(X, y)

        ###   Save the classifier   ###
        with open(f'{comm_args.output_path}_{modelName}.pkl', 'wb') as classifier_output:
            cPickle.dump(decision_function, classifier_output)

        ###   Save associated max feature values   ###
        min_max_features = {"maxX":maxX, "maxY": maxY, "minX":minX, "minY":minY}
        data = [min_max_features, commandline_arguments]
        with open(str(comm_args.output_path)+".json", 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)
        print(f"MinMax on X axis: {minX} {maxX} \nMinMax on Y axis: {minY} {maxY}")
        print(f"Model {modelName} has a score of {scores}")

        ###   Plot the classifier   ###
        if comm_args.plot_df:
            plot_decision_function(decision_function, X, y, modelName, modelIdx, plt.cm.gray, len(comm_args.model_type))
        modelIdx += 1
    if comm_args.plot_df:
        plt.tight_layout()
        #plt.tight_layout(h_pad=0.2, w_pad=0.2, pad=2.5)
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))