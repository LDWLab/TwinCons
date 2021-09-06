#!/usr/bin/env python3
'''Train a calibrated classifier.
'''

import numpy as np
import os, sys, argparse, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sklearn.calibration import CalibratedClassifierCV
from twincons.twcSVMtrain import train_classifier
from twincons.twcSVMtest import load_csv_data, csv_iterator, trim_data_by_top_segments
from twincons.twcCrossValidate import make_idx
import pickle as cPickle

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file with segment data for training and calibration', type=str)
    parser.add_argument('pickle', help='Path to output trained and calibrated classifer', type=str)
    parser.add_argument('-p',"--penalty", help='Penalty to train the classifier. Default: 1', type=float, default=1)
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to the top segments that cover\
        \nthis percentage of the total normalized length and weight. (Default = 0.5)', type=float, default=0.5)

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
    clf_score = brier_score_loss(y_test, prob_pos_clf, sample_weight=weights_test)
    print("No calibration: %1.3f" % clf_score)
    clf_isotonic_score = brier_score_loss(y_test, iso_clf_probs[:, 1], sample_weight=weights_test)
    print("With isotonic calibration: %1.3f" % clf_isotonic_score)
    clf_sigmoid_score = brier_score_loss(y_test, sig_clf_probs[:, 1], sample_weight=weights_test)
    print("With sigmoid calibration: %1.3f" % clf_sigmoid_score)

    if iso_score > sig_score:
        return calClassifierSig
    return calClassifierIso

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    #if comm_args.output_path == None:
    #    comm_args.output_path = comm_args.csv_path.replace('.csv', '')+'_calibrate'
    csv_list = csv_iterator(comm_args.csv_path)
    csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    X, y, sample_weight, maxX, maxY, minX, minY, aln_names = load_csv_data(csv_list)

    train_ind_alns, validate_ind_alns, test_ind_alns = make_idx(aln_names, 3)
    train_ind = [item for sublist in train_ind_alns for item in sublist]
    validate_ind = [item for sublist in validate_ind_alns for item in sublist]
    test_ind = [item for sublist in test_ind_alns for item in sublist]
    X_train, y_train, X_test, y_test, X_valid, y_valid = X[train_ind], y[train_ind], X[test_ind], y[test_ind], X[validate_ind], y[validate_ind]
    classifier = train_classifier(X_train, y_train, comm_args.penalty, 'auto', 'rbf', sample_weight=np.array(sample_weight)[train_ind])
    calibratedClassifier = calibrate(classifier, X_test, y_test, X_valid, y_valid, np.array(sample_weight)[test_ind])

    with open(comm_args.pickle, 'wb') as classifier_output:
        cPickle.dump(calibratedClassifier, classifier_output)

    min_max_features = {"maxX":maxX, "maxY": maxY, "minX":minX, "minY":minY}
    data = [min_max_features, commandline_arguments]
    with open(str(comm_args.pickle.replace('.pkl',''))+".json", 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



