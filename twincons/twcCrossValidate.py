#!/usr/bin/env python3
'''Cross validate penalty selection for a given file.
'''
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import os, sys, random, csv, argparse
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from twincons.twcSVMtest import load_and_assign_data, \
                                trim_data_by_top_segments, \
                                csv_iterator, mass_test, \
                                use_absolute_length_of_segments, \
                                normalize_features
from twincons.twcSVMtrain import train_classifier

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('-o', '--output_path', help='Path and name for output files. (Default = csv_path_crossval)', type=str, default=None)
    parser.add_argument('-abs', '--absolute_length', help='Use the segment absolute length as X coordinate. Does not weight segments.', 
                                    action="store_true", default=False)
    parser.add_argument('-nf', '--number_folds', help='Number of folds to split the training dataset into. (Default = 3)', type=int, default=3)
    parser.add_argument('-p', '--penalties', nargs='+', help='List of penalty options to evaluate. (Default = 0.05, 0.1, 1, 2, 5, 10, 20, 50, 100',
                                    default=[0.05, 0.1, 1, 2, 5, 10, 20, 50, 100], type=float)
    parser.add_argument('-dt', '--range_distances', nargs=3, metavar=('Start', 'Stop', 'Step'), 
                                    help='Range of distances from the decision boundary to evaluate.\
                                    \nDefault (-20, 20, 0.05).', default=[-20, 20, 0.05], type = float)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def chunks(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def shuffle2d(arr2d, rand=random):
    """Shuffes entries of 2-d array arr2d, preserving shape."""
    reshape = []
    data = []
    iend = 0
    for row in arr2d:
        data.extend(row)
        istart, iend = iend, iend+len(row)
        reshape.append((istart, iend))
    rand.shuffle(data)
    return [data[istart:iend] for (istart,iend) in reshape]

def make_idx(aln_names, num_of_splits):
    aln_ixs = dict()
    i = 0
    for entry in aln_names:
        if entry not in aln_ixs.keys():
            aln_ixs[entry] = list()
        aln_ixs[entry].append(i)
        i += 1
    mixed_list = list()
    for aln,pos in aln_ixs.items():
        num_chunks = len(pos)//num_of_splits
        mixed_list.append(shuffle2d(chunks(pos, num_of_splits)))
    return list(zip(*mixed_list))

def cv_by_alns(aln_names, num_of_splits):
    i = 0
    mixed_indexes = make_idx(aln_names, num_of_splits)
    flattened = list()
    for index in mixed_indexes:
        flattened.append(sorted([item for sublist in index for item in sublist]))
    for i,ele in enumerate(flattened):
        ind_test = ele
        ind_train = flattened[i+1:] + flattened[:i]
        ind_train = [item for sublist in ind_train for item in sublist]
        yield ind_train, ind_test

def predict_test_set(classifier, test_segment):
    segment_pred = classifier.predict(test_segment.reshape(1,-1))[0]
    segment_dist = classifier.decision_function(test_segment.reshape(1,-1))[0]
    segment_prob = classifier.predict_proba(test_segment.reshape(1,-1))[0][1]
    return segment_pred, segment_dist, segment_prob

def load_data(csv_location, top_segments=1, abs_length=False):
    csv_list = csv_iterator(csv_location)
    csv_list = trim_data_by_top_segments(csv_list, top_segments)
    if abs_length:
        csv_list = use_absolute_length_of_segments(csv_list)
    X, y, sample_weight, aln_names = load_and_assign_data(csv_list)
    return np.asarray(X), np.asarray(y), sample_weight, aln_names

def set_SVM_train_params(penalty=1, gamma='auto', kernel='rbf'):
    return penalty, gamma, kernel

def calc_stats_by_folds(aln_names, number_folds, X, y, sample_weight, penalty, gamma, kernel, start, stop, step):
    tprs, fprs, aucs = list(), list(), list()
    for i, (train_ind, test_ind) in enumerate(cv_by_alns(aln_names, number_folds)):
        tprs.append(list())
        fprs.append(list())
        segment_pred_dist = dict()
        X_test, y_test, X_train, y_train, = X[test_ind], y[test_ind], X[train_ind], y[train_ind]
        sample_weight_train = list(np.asarray(sample_weight)[train_ind])
        maxX, maxY, minX, minY = max(X_train[:, 0]), max(X_train[:, 1]), min(X_train[:, 0]), min(X_train[:, 1])
        X_train_norm = np.asarray(normalize_features(list(zip(X_train[:,0], X_train[:,1])), maxX, maxY, minX, minY))

        classifier = train_classifier(X_train_norm, y_train, penalty, gamma, kernel, sample_weight=sample_weight_train)

        for segment, aln_ind in zip(X_test, test_ind):
            test_segment = np.array([float(segment[0])/float(maxX),float(segment[1])/float(maxY)])
            segment_pred, segment_dist, segment_prob = predict_test_set(classifier, test_segment)
            if aln_names[aln_ind] not in segment_pred_dist.keys():
                segment_pred_dist[aln_names[aln_ind]] = list()
            segment_pred_dist[aln_names[aln_ind]].append(['',(segment_pred, segment_dist,'')])
        dist_to_stats = mass_test(segment_pred_dist, min_threshold=start, max_threshold=stop, step=step)
        for dist, stats in dist_to_stats.items():
            tprs[i].append(stats[0])
            fprs[i].append(1-stats[1])
        aucs.append(auc(fprs[i], tprs[i]))

    mean_fpr = np.mean(np.array(fprs), axis=0)
    mean_tpr = np.mean(np.array(tprs), axis=0)
    std_tpr = np.std(tprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    return mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr

def plot_roc_curve(axis, mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr, label='', color=''):
    if label == '':
        label = "Mean ROC"
    if color == '':
        color_line = 'blue'
        color_std = 'grey'
    else:
        color_line = plt.cm.viridis(color)
        color_std = plt.cm.viridis(color)
    
    axis.plot(mean_fpr, mean_tpr, 
            label=f'{label} (AUC = {mean_auc:.2f} $\pm$ {std_auc:.2f})',
            lw=2, alpha=.8, color=color_line)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    axis.fill_between(mean_fpr, tprs_lower, tprs_upper, color=color_std, alpha=.1)

def write_stats_csv(penalty_to_stats, range_of_dists, csv_loc):
    keys = sorted(penalty_to_stats.keys())
    penalty_row = list(np.repeat(keys,3))
    distance_row = ['TPR','FPR','STD-tpr']*len(keys)
    penalty_row.insert(0, 'Penalty')
    distance_row.insert(0, 'Distance')
    with open(csv_loc, mode='w') as output_csv:
        csv_writer = csv.writer(output_csv, delimiter=',')
        csv_writer.writerow(penalty_row)
        csv_writer.writerow(distance_row)
        for i, dist in enumerate(range_of_dists):
            temprow = [dist]
            for penalty in keys:
                temprow.extend([penalty_to_stats[penalty][0][i], penalty_to_stats[penalty][1][i], penalty_to_stats[penalty][4][i]])
            csv_writer.writerow(temprow)

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    if comm_args.output_path == None:
        comm_args.output_path = comm_args.csv_path.replace('.csv', '')+'_crossval'
    X, y, sample_weight, aln_names = load_data(comm_args.csv_path, top_segments=1, abs_length=comm_args.absolute_length)
    number_folds = comm_args.number_folds
    penalties = comm_args.penalties
    start_dist, stop_dist, step_dist = comm_args.range_distances[0], comm_args.range_distances[1], comm_args.range_distances[2],

    penalty_to_stats = dict()
    for var_penalty in penalties:
        penalty, gamma, kernel = set_SVM_train_params(var_penalty, 'auto', 'rbf')
        mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr = calc_stats_by_folds(aln_names, 
                                                                            number_folds, 
                                                                            X, y, sample_weight,
                                                                            penalty, gamma, kernel,
                                                                            start_dist, stop_dist, step_dist)
        penalty_to_stats[var_penalty] = (mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr)

    write_stats_csv(penalty_to_stats, np.arange(start_dist, stop_dist, step_dist), comm_args.output_path+'.csv')
    
    fig, ax = plt.subplots()
    color_indexes = np.linspace(0, 1, len(penalties))
    
    for i, (penalty, stats) in enumerate(penalty_to_stats.items()):
        (mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr) = stats
        plot_roc_curve(ax, mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr,
                         label = f'Penalty: {penalty}', 
                         color = color_indexes[i])

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
           title="Receiver operating characteristic")
    ax.legend(loc=4)
    plt.savefig(comm_args.output_path+'.svg', dpi=600)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



