#!/usr/bin/env python3
'''Cross validate penalty selection for a given file.
'''
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import os, sys, csv, argparse
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier,
                              AdaBoostClassifier)
from sklearn.tree import DecisionTreeClassifier
from twincons.twcSVMtest import normalize_features
from twincons.twcCrossValidate import cv_by_alns, load_data, plot_roc_curve
from twincons.twcTreesTrain import train_classifier
from twincons.twcTreesTest import mass_test

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('-o', '--output_path', help='Path and name for output files. (Default = csv_path_crossval)', type=str, default=None)
    parser.add_argument('-abs', '--absolute_length', help='Use the segment absolute length as X coordinate. Does not weight segments.', 
                                    action="store_true", default=False)
    parser.add_argument('-nf', '--number_folds', help='Number of folds to split the training dataset into. (Default = 3)', type=int, default=3)
    parser.add_argument('-ne', '--n_estimators', nargs='+', help='List of numbers of estimators to evaluate. (Default = 2, 5, 10, 30, 60, 100',
                                    default=[2, 5, 10, 30, 60, 100], type=int)
    parser.add_argument('-ps', '--probability_step', help='Probability step used for ROC calculation when identifying significant segments.\
        \nMakes multiple evaluations of the data moving from 0 to 1 probabilities with the given step.Default (0.001)', default=0.001, type = float)
    parser.add_argument('-mt','--model_type', help='Select the model of the classifier. TODO: List multiple for comparison. (Deafult: RandomForest)', 
                        nargs='+', choices=['RandomForest','ExtraTrees', 'AdaBoost'], default=['RandomForest'])
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def predict_test_set(model, test_segment):
    segment_pred = model.predict(np.array(test_segment).reshape(1,-1))[0]
    segment_prob = model.predict_proba(np.array(test_segment).reshape(1,-1))[0][1]
    return segment_pred, segment_prob

def getModel(modelName, nEst, maxDepth=7):
    modelDict = dict(RandomForest = RandomForestClassifier(n_estimators=nEst),
                     ExtraTrees = ExtraTreesClassifier(n_estimators=nEst),
                     AdaBoost = AdaBoostClassifier(DecisionTreeClassifier(max_depth=maxDepth), n_estimators=nEst))
    return modelDict[modelName]

def calc_stats_by_folds(aln_names, number_folds, X, y, sample_weight, model, probStep):
    tprs, fprs, aucs = list(), list(), list()
    for i, (train_ind, test_ind) in enumerate(cv_by_alns(aln_names, number_folds)):
        tprs.append(list())
        fprs.append(list())
        segment_pred_prob = dict()
        X_test, y_test, X_train, y_train, = X[test_ind], y[test_ind], X[train_ind], y[train_ind]
        sample_weight_train = list(np.asarray(sample_weight)[train_ind])
        
        maxX, maxY, minX, minY = max(X_train[:, 0]), max(X_train[:, 1]), min(X_train[:, 0]), min(X_train[:, 1])
        X_train_norm = np.asarray(normalize_features(list(zip(X_train[:,0], X_train[:,1])), maxX, maxY, minX, minY))
        maxX, maxY, minX, minY = max(X_test[:, 0]), max(X_test[:, 1]), min(X_test[:, 0]), min(X_test[:, 1])
        X_test_norm = np.asarray(normalize_features(list(zip(X_test[:,0], X_test[:,1])), maxX, maxY, minX, minY))

        classifier = train_classifier(X_train_norm, y_train, model, sample_weight=sample_weight_train)

        for test_segment, aln_ind in zip(X_test_norm, test_ind):
            segment_pred, segment_prob = predict_test_set(classifier, test_segment)
            if aln_names[aln_ind] not in segment_pred_prob.keys():
                segment_pred_prob[aln_names[aln_ind]] = list()
            segment_pred_prob[aln_names[aln_ind]].append(['',(segment_pred, segment_prob,'')])
        prob_to_stats = mass_test(segment_pred_prob, step=probStep)
        for prob, stats in prob_to_stats.items():
            tprs[i].append(stats[0])
            fprs[i].append(1-stats[1])
        aucs.append(auc(fprs[i], tprs[i]))

    mean_fpr = np.mean(np.array(fprs), axis=0)
    mean_tpr = np.mean(np.array(tprs), axis=0)
    std_tpr = np.std(tprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    return mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr

def write_stats_csv(penalty_to_stats, range_of_probs, csv_loc):
    keys = sorted(penalty_to_stats.keys())
    penalty_row = list(np.repeat(keys,3))
    distance_row = ['TPR','FPR','STD-tpr']*len(keys)
    penalty_row.insert(0, 'Number estimators')
    distance_row.insert(0, 'Probability')
    with open(csv_loc, mode='w') as output_csv:
        csv_writer = csv.writer(output_csv, delimiter=',')
        csv_writer.writerow(penalty_row)
        csv_writer.writerow(distance_row)
        for i, prob in enumerate(range_of_probs):
            temprow = [prob]
            for penalty in keys:
                temprow.extend([penalty_to_stats[penalty][0][i], penalty_to_stats[penalty][1][i], penalty_to_stats[penalty][4][i]])
            csv_writer.writerow(temprow)

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    if comm_args.output_path == None:
        comm_args.output_path = comm_args.csv_path.replace('.csv', '')+'_crossval'
    X, y, sample_weight, aln_names = load_data(comm_args.csv_path, top_segments=1, abs_length=comm_args.absolute_length)
    number_folds = comm_args.number_folds
    n_estimators = comm_args.n_estimators
    probStep = float(comm_args.probability_step)
    modelName = comm_args.model_type[0]

    nEst_to_stats = dict()
    for numberEstimator in n_estimators:
        mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr = calc_stats_by_folds(aln_names, 
                                                                            number_folds, 
                                                                            X, y, sample_weight,
                                                                            getModel(modelName, numberEstimator),
                                                                            probStep)
        nEst_to_stats[numberEstimator] = (mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr)

    write_stats_csv(nEst_to_stats, np.arange(0, 1+2*probStep, probStep), comm_args.output_path+'.csv')
    
    fig, ax = plt.subplots()
    color_indexes = np.linspace(0, 1, len(n_estimators))
    
    for i, (nEstimators, stats) in enumerate(nEst_to_stats.items()):
        (mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr) = stats
        plot_roc_curve(ax, mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr,
                         label = f'N estimators: {nEstimators}', 
                         color = color_indexes[i])

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
           title="Receiver operating characteristic")
    ax.legend(loc=4)
    plt.savefig(comm_args.output_path+'.svg', dpi=600)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



