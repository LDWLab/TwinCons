import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import auc, plot_roc_curve
from sklearn import svm
import os, sys, random
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bin.SVM_test import load_and_assign_data, trim_data_by_top_segments, csv_iterator, use_absolute_length_of_segments, mass_test
from bin.SVM_train import train_classifier

#######################
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
    return segment_pred, segment_dist

def normalize_train_set(X, Y, maxX, maxY):
    norm_xy = list()
    for x, y in zip(X, Y):
        norm_xy.append([float(x)/float(maxX), float(y)/float(maxY)])
    return np.asarray(norm_xy)

def load_data(csv_location, top_segments=1, abs_length=False):
    csv_list = csv_iterator(csv_location)
    csv_list = trim_data_by_top_segments(csv_list, top_segments)
    if abs_length:
        csv_list = use_absolute_length_of_segments(csv_list)
    X, y, sample_weight, aln_names = load_and_assign_data(csv_list)
    return np.asarray(X), np.asarray(y), sample_weight, aln_names

def set_SVM_train_params(penalty=1, gamma='auto', kernel='rbf'):
    return penalty, gamma, kernel

def calc_stats_by_folds(aln_names, number_folds, X, y, penalty, gamma, kernel):
    tprs, fprs, aucs = list(), list(), list()
    for i, (train_ind, test_ind) in enumerate(cv_by_alns(aln_names, number_folds)):
        tprs.append(list())
        fprs.append(list())
        segment_pred_dist = dict()
        X_test, y_test, X_train, y_train, = X[test_ind], y[test_ind], X[train_ind], y[train_ind]
        maxX, maxY = max(X_train[:, 0]), max(X_train[:, 1])
        X_train_norm = normalize_train_set(X_train[:,0], X_train[:,1], maxX, maxY)

        classifier = train_classifier(X_train_norm, y_train, penalty, gamma, kernel)

        for segment, aln_ind in zip(X_test, test_ind):
            test_segment = np.array([float(segment[0])/float(maxX),float(segment[1])/float(maxY)])
            segment_pred, segment_dist = predict_test_set(classifier, test_segment)
            if aln_names[aln_ind] not in segment_pred_dist.keys():
                segment_pred_dist[aln_names[aln_ind]] = list()
            segment_pred_dist[aln_names[aln_ind]].append(['',(segment_pred, segment_dist,'')])
        dist_to_stats = mass_test(segment_pred_dist, min_threshold=-2, max_threshold=2, step=0.01)
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

def main(commandline_arguments):
    X, y, sample_weight, aln_names = load_data("./data/CSV/PRST_cg09_it1_lt3.csv", top_segments=0.5, abs_length=False)
    number_folds = 3
    #penalties = [0.1, 1, 2, 5, 10, 20, 50, 100]
    penalties = [0.1, 1, 2]

    penalty_to_stats = dict()
    for var_penalty in penalties:
        penalty, gamma, kernel = set_SVM_train_params(var_penalty, 'auto', 'rbf')
        mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr = calc_stats_by_folds(aln_names, 
                                                                            number_folds, 
                                                                            X, y, 
                                                                            penalty, gamma, kernel)
        penalty_to_stats[var_penalty] = (mean_tpr, mean_fpr, mean_auc, std_auc, std_tpr)
    
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
    plt.savefig("./PRST_ROC_3folds.png", dpi=600)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



