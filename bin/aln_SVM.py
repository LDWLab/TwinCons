#!/usr/bin/env python3
"""
Generate SVM from alignment segments.
Two modes:
    Train: Computes a decision function from csv generated with MultiPhyMeas
    Test: Evaluates alignment entries in csv with a provided decision function
"""
#print(__doc__)
import re, sys, csv, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import svm
from statistics import mean 
import _pickle as cPickle

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('-tr','--train',help='Provide path to save classifier as a pickle binary file.', type=str)
    mode.add_argument('-te','--test',help='Provide path to classifier pickle binary file.', type=str)
    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-mf','--max_features',help='Provide path to file with comma separated maximal values \
                                                    \nfor each feature class used to generate the classifier.', type=str)
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to longest N segments.', type=int)
    parser.add_argument('-tcp','--test_classifier_precision', help='Provided csv is annotated for testing the classifier. \
                                        \nChoose between two ways of determining a positive result:\
                                        \n- by average distance of segments fom the decision boundary (dist);\
                                        \n- by the presence of at least one positive segment (id).', choices=['dist', 'id'])
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

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

def trim_data_by_top_segments(csv_list, number_top_segments):
    '''Returns trimmed down list by number of longest segments'''
    from heapq import nlargest
    aln_id_to_data = dict()
    for row in csv_list:
        if row[0] not in aln_id_to_data.keys():
            aln_id_to_data[row[0]] = list()
        aln_id_to_data[row[0]].append((row[2], row[3], row[1], row[4]))
    output_list = list()
    for aln_id, data_items in aln_id_to_data.items():
        for segment_data in nlargest(number_top_segments, data_items):
            output_list.append([aln_id, segment_data[2], segment_data[0], segment_data[1], segment_data[3]])
    return output_list

def load_csv_data(comm_args, csv_list, max_features=''):
    '''
    Reads a csv outputted from MultiPhyMeas.py, entries starting with A_ and B_ are considered as + data
    and entries starting with C_ and D_ as - data. Normalizes all the data for x and y in the range 0,1
    using the max(x) and max(y).
    '''
    data_xy = []
    data_identity = []
    data_weights = []
    aln_names = []
    for row in csv_list:
        aln_name = row[0]
        aln_names.append(aln_name)
        if re.match(r'^A_|^B_',aln_name):
            data_xy.append([float(row[2]),float(row[3])])
            data_identity.append(1)
            data_weights.append(row[1])
        elif re.match(r'^C_|^D_',aln_name):
            data_xy.append([float(row[2]),float(row[3])])
            data_identity.append(0)
            data_weights.append(float(row[1]))
    data_xy_normx = []
    if max_features == '':
        maxX = max(np.asarray(data_xy)[:,0])
        maxY = max(np.asarray(data_xy)[:,1])
    else:
        maxX = float(max_features[0])
        maxY = float(max_features[1])
    for tups in data_xy:
        if comm_args.train:
            data_xy_normx.append([float(tups[0]/maxX),float(tups[1]/maxY)])
        elif comm_args.test:
            data_xy_normx.append([float(tups[0])/maxX,float(tups[1])/maxY])
    return np.asarray(data_xy_normx), np.asarray(data_identity), data_weights, maxX, maxY, aln_names

def test_function(csv_list, decision_function, max_features):
    '''
    Executes prediction and distance calculation on each
    segment from the input for a given decision function.
    '''
    segment_pred_dist = {}
    for entry in csv_list:
        test_segment = np.array([float(entry[2])/float(max_features[0]),float(entry[3])/float(max_features[1])])
        segment_pred = decision_function.predict(test_segment.reshape(1,-1))[0]
        segment_dist = decision_function.decision_function(test_segment.reshape(1,-1))[0]
        if str(entry[0]) not in segment_pred_dist.keys():
            segment_pred_dist[str(entry[0])] = []
        segment_pred_dist[str(entry[0])].append([entry[4],(segment_pred,segment_dist)])
    return segment_pred_dist

def flatten_alignment_segment_stats_to_groups(segment_pred_dist):
    grouped_data={}
    for aln_name in segment_pred_dist:
        dist_sum, pos_sum, tot_segm = 0,0,0
        if aln_name.split("_")[0] not in grouped_data.keys():
            grouped_data[aln_name.split("_")[0]] = []
        for segm in segment_pred_dist[aln_name]:
            tot_segm += 1
            dist_sum += segm[1][1]
            if segm[1][0] == 1:
                pos_sum += 1
        print(aln_name, dist_sum, pos_sum, tot_segm)
        grouped_data[aln_name.split("_")[0]].append((dist_sum, pos_sum, tot_segm))
    return grouped_data

def bypass_zero_division(x,y):
    try:
        return x/y
    except ZeroDivisionError:
        return 0

def calc_avedist_stats(grouped_data, thr, tp, tn, fp, fn):
    for group, segments in grouped_data.items():
        if group == 'A' or group == 'B':
            fn += sum(i[0]/i[2] <  thr for i in segments)
            tp += sum(i[0]/i[2] >= thr for i in segments)
        if group == 'C' or group == 'D':
            tn += sum(i[0]/i[2] <  thr for i in segments)
            fp += sum(i[0]/i[2] >= thr for i in segments)
    return tp, tn, fp, fn

def calc_identity_stats(grouped_data, tp, tn, fp, fn):
    for group, segments in grouped_data.items():
        if group == 'A' or group == 'B':
            for i in segments:
                if i[1] == 0:
                    fn += 1
                else:
                    tp += 1
        if group == 'C' or group == 'D':
            for i in segments:
                if i[1] == 0:
                    tn += 1
                else:
                    fp += 1
    print(tp, tn, fp, fn)
    return tp, tn, fp, fn

def mass_test(segment_pred_dist, grouped_data, distance_or_identity, min_threshold=0,max_threshold=2, step=0.1):
    '''For evaluating the decision function with 
    great number of alignments in separate groups.
    '''
    ###   Calculate Sensitivity(tpr), Specificity(tnr) and Precision   ###
    ###   over a range of distances from the decision boundary         ###
    dist_to_stats = {}
    for thr in np.arange(min_threshold,max_threshold,step):
        tp, tn, fp, fn = 0,0,0,0
        if distance_or_identity == 'dist':
            tp, tn, fp, fn = calc_avedist_stats(grouped_data, thr, tp, tn, fp, fn)
        elif distance_or_identity == 'id':
            tp, tn, fp, fn = calc_identity_stats(grouped_data, tp, tn, fp, fn)
        tpr = bypass_zero_division(tp,tp+fn)
        tnr = bypass_zero_division(tn,tn+fp)
        precision = bypass_zero_division(tp,tp+fp)
        dist_to_stats[thr] = (tpr, tnr, precision)
        print ("Threshold "+str(thr),'\n',"tpr", tpr,"tnr", tnr,'\n',"precision",precision)
    return dist_to_stats

def plot_decision_function(classifier, X, y, sample_weight, axis, fig, title, aln_names, distance_or_identity='',  decision_levels=''):
    '''
    Plots the data and the decision function. Besides the classifier function, takes in the sample weights
    and plots it as a size of the datapoints. If they are different than 1.
    '''
    # plot the decision function
    xx, yy = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    vir_cmap = plt.cm.get_cmap('viridis')

    ###   Draws each decision level (if present) with color from the viridis colormap   ###
    if decision_levels is not '' and distance_or_identity == 'dist':
        # Specifies under and over values to first and last color of the colormap
        vir_cmap.set_under(vir_cmap(0))
        vir_cmap.set_over(vir_cmap(1))
        CS1 = axis.contourf(xx, yy, Z, sorted(decision_levels.keys()), 
                            colors=vir_cmap(np.linspace(0, 1, len(decision_levels))),
                            extend='both')
        cbar = fig.colorbar(CS1, ticks=sorted(decision_levels.keys()))
        levels_labels = ["%.2f, %.2f, %.2f" %decision_levels[thr] for thr in sorted(decision_levels.keys())]
        cbar.set_ticklabels(["%.1f %s" % (thr,lev) for lev,thr in zip(levels_labels, sorted(decision_levels.keys()))])
        ###   Draws a line for each decision level   ###
        if (sorted(decision_levels.keys())[1] - sorted(decision_levels.keys())[0]) >= 0.5:
            CS2 = axis.contour(xx, yy, Z, colors='black', levels=sorted(decision_levels.keys()),
                              alpha=0.75, linestyles=['-'], linewidths=0.5)
            axis.clabel(CS2, fmt='%2.1f', colors='black', fontsize=4)

    ###   Draws the decision function as a red line   ###
    else:
        axis.contour(xx, yy, Z, colors='r', levels=[0], linestyles=['-'], linewidths=0.5)

    abs_length = [float(n)**2 for n in sample_weight]
    if distance_or_identity == 'id':
        title +=" with tpr "+"{:2.2f}".format(decision_levels[next(iter(decision_levels))][0])\
                +" and tnr "+"{:2.2f}".format(decision_levels[next(iter(decision_levels))][1])\
                +" per alignment"
        mappedDict = dict(zip(sorted(set(aln_names)), range(len(aln_names))))
        aln_identities = list(map(lambda x: mappedDict[x], aln_names))
        scatter = axis.scatter(X[:, 0], X[:, 1], c=aln_identities, alpha=0.75, s=abs_length,
                 cmap=plt.cm.viridis, edgecolors='black')
        inv_map = {v: k for k, v in mappedDict.items()}
        aln_labels = list()
        for label in scatter.legend_elements()[1]:
            aln_labels.append(inv_map[int(re.findall(r'\d+', label)[0])])
        plt.legend(scatter.legend_elements()[0], aln_labels,
                    loc="upper right", title="Alignments")

    if distance_or_identity == '':
        axis.scatter(X[:, 0], X[:, 1], c=y, alpha=0.75, s=abs_length,
                 cmap=plt.cm.bone, edgecolors='black')
    #axis.axis('off')
    axis.set_title(title)

def plot_test_histograms(grouped_data):
    '''Some plots to visualize different group results'''
    fig, axes = plt.subplots(2, 2, sharex=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for ax, group in zip(axes.ravel(), sorted(grouped_data.keys())):
       ax.hist(grouped_data[group])
       ax.set_title(group)
    plt.savefig('./data/outputs/CSV/2ftesting.png', dpi=600)
    plt.clf()

    list_for_hist = []
    for group in sorted(grouped_data.keys()):
       list_for_hist.append(grouped_data[group])
    fig, axes = plt.subplots(1, 1)
    axes.hist(list_for_hist, label=sorted(grouped_data.keys()),histtype='barstacked', bins=30)
    axes.legend(prop={'size': 10})
    plt.savefig('./data/outputs/CSV/2ftesting_hists.png', dpi=600)
    plt.clf()

def read_max_features(max_features_path):
    with open(max_features_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for first_row in csv_reader:
            return first_row

def train_classifier(comm_args, X, y, maxX, maxY, penalty=10, gamma='auto', sample_weight=''):

    ###   Fit the classifier   ###
    decision_function = svm.SVC(C=penalty, gamma=gamma)
    if sample_weight != '':
        decision_function.fit(X, y, sample_weight=sample_weight)
    else:
        decision_function.fit(X, y)
    
    ###   Save max values of features   ###
    with open(str(comm_args.train)+".maxvals", 'w') as max_features:
        csv_writer = csv.writer(max_features, delimiter=',')
        csv_writer.writerow([maxX, maxY])
    print("Max on X axis:", maxX, "\nMax on Y axis:", maxY)
    
    ###   Save the classifier   ###
    with open(comm_args.train, 'wb') as fid:
        cPickle.dump(decision_function, fid)
    
    return decision_function

def main(commandline_arguments):
    '''Main entry point'''
    comm_args = create_and_parse_argument_options(commandline_arguments)

    ###   Load alignment segment data   ###
    csv_list = csv_iterator(comm_args.csv_path)
    if comm_args.top_segments:
        csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)

    ###   Testing mode   ###
    if comm_args.test:
        if not comm_args.max_features:
            parser.print_help()
            raise ValueError("In test mode providing maximal values for the classifier features is required!")
        max_features = read_max_features(comm_args.max_features)
        decision_function = cPickle.load(open(comm_args.test, 'rb'))
        segment_pred_dist = test_function(csv_list, decision_function, max_features)
        if comm_args.test_classifier_precision:
            grouped_data = flatten_alignment_segment_stats_to_groups(segment_pred_dist)
            dist_to_se_sp_pr = mass_test(segment_pred_dist, grouped_data, 
                comm_args.test_classifier_precision, min_threshold=-2.5, max_threshold=5, step=0.5)
            #plot_test_histograms(grouped_data)

    ###   Training mode  ###
    if comm_args.train:
        X, y, sample_weight, maxX, maxY, aln_names = load_csv_data(comm_args, csv_list)
        decision_function = train_classifier(comm_args, X, y, maxX, maxY, sample_weight=sample_weight)

    ###   Plot the classifier   ###
    if comm_args.plot_df:
        fig, axes = plt.subplots(1, 1, )
        if comm_args.test:
            X, y, sample_weight, maxX, maxY, aln_names = load_csv_data(comm_args, csv_list, max_features=max_features)
            if comm_args.test_classifier_precision:
                plot_decision_function(decision_function,X,y, sample_weight, axes, fig,
                                    "Decision function",aln_names, comm_args.test_classifier_precision, decision_levels=dist_to_se_sp_pr)
            else:
                plot_decision_function(decision_function,X,y, sample_weight, axes, fig,
                                    "Decision function",aln_names, decision_levels=dist_to_se_sp_pr)
        if comm_args.train:
            plot_decision_function(decision_function,X,y, sample_weight, axes, fig, aln_names,
                                    "Decision function")
        plt.tight_layout()
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))