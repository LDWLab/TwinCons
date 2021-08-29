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
from sklearn.calibration import calibration_curve, CalibratedClassifierCV
from twincons.twcSVMtest import normalize_features
from twincons.twcCrossValidate import load_data, plot_roc_curve, make_idx
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
    parser.add_argument('-mt','--model_type', help='Select the model of the classifier. List multiple for comparison. (Deafult: RandomForest)', 
                        nargs='+', choices=['RandomForest','ExtraTrees', 'AdaBoost'], default=['RandomForest'])
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def calibrate(clf, X_test, y_test, X_valid, y_valid, weights_test):
    from sklearn.metrics import log_loss, brier_score_loss
    # Gaussian Naive-Bayes with no calibration
    prob_pos_clf = clf.predict_proba(X_test)[:, 1]

    # Tree with Isotonic calibration
    calClassifierIso = CalibratedClassifierCV(clf, method="isotonic", cv="prefit")
    calClassifierIso.fit(X_valid, y_valid)
    iso_clf_probs = calClassifierIso.predict_proba(X_test)
    iso_score = log_loss(y_test, iso_clf_probs)


    # Gaussian Naive-Bayes with sigmoid calibration
    calClassifierSig = CalibratedClassifierCV(clf, method="sigmoid", cv="prefit")
    calClassifierSig.fit(X_valid, y_valid)
    sig_clf_probs = calClassifierSig.predict_proba(X_test)
    sig_score = log_loss(y_test, iso_clf_probs)

    print("Log loss scores: (the smaller the better)")
    print("With isotonic calibration: %1.3f" % iso_score)
    print("With sigmoid calibration: %1.3f" % sig_score)

    print("Brier scores: (the smaller the better)")
    clf_score = brier_score_loss(y_test, prob_pos_clf, weights_test)
    print("No calibration: %1.3f" % clf_score)
    clf_isotonic_score = brier_score_loss(y_test, iso_clf_probs[:, 1], weights_test)
    print("With isotonic calibration: %1.3f" % clf_isotonic_score)
    clf_sigmoid_score = brier_score_loss(y_test, sig_clf_probs[:, 1], weights_test)
    print("With sigmoid calibration: %1.3f" % clf_sigmoid_score)

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

    write_stats_csv(nEst_to_stats, np.arange(0, 1+probStep, probStep), comm_args.output_path+'.csv')
    
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



