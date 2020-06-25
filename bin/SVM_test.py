#!/usr/bin/env python3
"""
Evaluates alignment entries in csv generated from MultiPhyMeas.
Requires a decision function and max features generated from SVM_train.
Train and test only with the same parameters!
Such parameters can be % cutting gaps, center mass segments, top segments.
"""
#print(__doc__)
import re, sys, csv, math, argparse, json
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
    parser.add_argument('pickle',help='Provide path to classifier pickle binary file. The script will also search for an identically\
                                    \nnamed file with extension ".maxvals" containing comma separated maximal values, for example:\
                                    \npickle file:\trandom_test.pkl\
                                    \nmaximal values:\trandom_test.pkl.maxvals', type=str)
    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-et','--evalue_threshold', help='Threshold to consider an alignment significant. (Default 0.05)',
                                                                                            type=float, default=0.05)
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to the top segments that cover\
                                    \n this percentage of the total normalized length and weight. (Default = 0.5)', 
                                                                                            type=float, default=0.5)
    parser.add_argument('-cms','--center_mass_segments', help='Use the average position (Center of mass) \
                                    \nfrom all segments per alignment. No sample weighing.', action="store_true")
    parser.add_argument('-tcp','--test_classifier_precision', help='Provided csv is annotated for testing the classifier. \
                                    \nChoose between two ways of determining a positive result:\
                                    \n- by average distance of segments fom the decision boundary (dist);\
                                    \n- by the presence of at least one positive segment (id).', choices=['dist', 'id'])
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
        aln_id_to_data[row[0]].append((float(row[2]), float(row[3]), float(row[1]), float(row[4])))
    return aln_id_to_data

def trim_data_by_top_segments(csv_list, perc_top_segments):
    '''Returns trimmed down list by percentage coverage of the 
    total segment weights and total normalized segment lengths '''
    from heapq import nlargest
    from operator import itemgetter
    aln_id_to_data = csv_to_segments_annotated_by_aln(csv_list)
    output_list = list()
    for aln_id, data_items in aln_id_to_data.items():
        limit_weight = sum(x[1] for x in data_items)*perc_top_segments
        limit_nlength = sum(x[0] for x in data_items)*perc_top_segments
        curr_weight, curr_nlength = 0, 0
        for segment_data in sorted(data_items, key=lambda item: item[1:2], reverse=True):
            output_list.append([aln_id, segment_data[2], segment_data[0], segment_data[1], segment_data[3]])
            curr_weight  += segment_data[1]
            curr_nlength += segment_data[0]
            if curr_weight > limit_weight and curr_nlength > limit_nlength:
                break
    return output_list

def average(xs):
    N = float(len(xs))
    return tuple(sum(col)/N for col in zip(*xs))

def recalculate_data_by_averaging_segments(csv_list):
    aln_id_to_data = csv_to_segments_annotated_by_aln(csv_list)
    output_list = list()
    for aln_id, data_items in aln_id_to_data.items():
        averaged_results = average(data_items)
        output_list.append([aln_id, averaged_results[0]*100, averaged_results[2], averaged_results[1], 'NA'])
    return output_list

def load_csv_data(csv_list, max_features=''):
    '''
    Reads a csv outputted from MultiPhyMeas.py, entries starting with A_ and B_ are considered as + data
    and entries starting with C_ and D_ as - data. Normalizes all the data for x and y in the range 0,1
    using the max(x) and max(y).
    '''
    data_xy = []
    data_identity = []
    data_weights = []
    aln_names = []
    i = 0
    curr_name = ''
    for row in csv_list:
        aln_name = row[0].split(".")[0]
        aln_names.append(aln_name)
        data_xy.append([float(row[2]),float(row[3])])
        data_weights.append(float(row[1]))
        if re.match(r'^A_|^B_',aln_name):
            data_identity.append(1)
        elif re.match(r'^C_|^D_',aln_name):
            data_identity.append(0)
        else:
            data_identity.append(i)
            if curr_name != row[0]:
                curr_name = row[0]
                i+=1
    data_xy_normx = []
    if max_features == '':
        maxX = max(np.asarray(data_xy)[:,0])
        maxY = max(np.asarray(data_xy)[:,1])
    else:
        maxX = float(max_features[0])
        maxY = float(max_features[1])
    for tups in data_xy:
        data_xy_normx.append([float(tups[0])/maxX,float(tups[1])/maxY])
    return np.asarray(data_xy_normx), np.asarray(data_identity), data_weights, maxX, maxY, aln_names

def test_function(csv_list, classifier, max_features):
    '''
    Executes prediction and distance calculation on each
    segment from the input for a given decision function.
    '''
    segment_pred_dist = {}
    for entry in csv_list:
        test_segment = np.array([float(entry[2])/float(max_features[0]),float(entry[3])/float(max_features[1])])
        segment_pred = classifier.predict(test_segment.reshape(1,-1))[0]
        segment_dist = classifier.decision_function(test_segment.reshape(1,-1))[0]
        if str(entry[0]) not in segment_pred_dist.keys():
            segment_pred_dist[str(entry[0])] = []
        segment_pred_dist[str(entry[0])].append([entry[4],(segment_pred,segment_dist,entry[1])])
    return segment_pred_dist

def flatten_alignment_segment_stats_to_groups(segment_pred_dist, by_group=False):
    grouped_data={}
    for aln_name in segment_pred_dist:
        if by_group:
            outname = aln_name.split("_")[0]
        else:
            outname = aln_name.split(".")[0]
        dist_sum, pos_sum, tot_segm = 0,0,0
        if outname not in grouped_data.keys():
            grouped_data[outname] = []
        for segm in segment_pred_dist[aln_name]:
            tot_segm += 1
            dist_sum += segm[1][1]*segm[1][2]
            if segm[1][0] == 1:
                pos_sum += 1
        grouped_data[outname].append((dist_sum, pos_sum, tot_segm))
    return grouped_data

def bypass_zero_division(x,y):
    if y != 0:
        return x/y
    else:
        return 0

def calc_avedist_stats(grouped_data, thr, tp, tn, fp, fn):
    for group, alignments in grouped_data.items():
        if group == 'A' or group == 'B':
            fn += sum(aln[0]/aln[2] <  thr for aln in alignments)
            tp += sum(aln[0]/aln[2] >= thr for aln in alignments)
        if group == 'C' or group == 'D':
            tn += sum(aln[0]/aln[2] <  thr for aln in alignments)
            fp += sum(aln[0]/aln[2] >= thr for aln in alignments)
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

def draw_thresholds(axis, fig, X, xx, yy, Z, decision_levels, clean=False):

    vir_cmap = plt.cm.get_cmap('viridis')
    # Specifies under and over values to first and last color of the colormap
    vir_cmap.set_under(vir_cmap(0))
    vir_cmap.set_over(vir_cmap(1))
    CS1 = axis.contourf(xx, yy, Z, sorted(decision_levels.keys()), 
                        colors=vir_cmap(np.linspace(0, 1, len(decision_levels))),
                        extend='both')
    cbar = fig.colorbar(CS1, ticks=sorted(decision_levels.keys()))
    if clean:
        cbar.set_ticklabels(["%.1f" %thr for thr in sorted(decision_levels.keys())])
        cbar.ax.set_title('Distance')
    else:
        levels_labels = ["%.2f, %.2f, %.2f" %decision_levels[thr] for thr in sorted(decision_levels.keys())]
        cbar.set_ticklabels(["%.1f %s" % (thr,lev) for lev,thr in zip(levels_labels, sorted(decision_levels.keys()))])
        cbar.ax.set_title('Dist  tpr, tnr, prec', x=5)
    ###   Draws a line for each decision level   ###
    if math.ceil(max(X[:, 1])) == 1 and math.ceil(max(X[:, 0])) == 1:
        CS2 = axis.contour(xx, yy, Z, colors='magenta', levels=sorted(decision_levels.keys()),
                          alpha=0.75, linestyles=[':'], linewidths=0.5)
        axis.clabel(CS2, fmt='%2.1f', colors='magenta', fontsize=5)

def plot_decision_function(classifier, X, y, sample_weight, axis, fig, title, aln_names, label_order_tups=None, distance_or_identity='', decision_levels=''):
    '''
    Plots the data and the decision function. Besides the classifier function, takes in the sample weights
    and plots it as a size of the datapoints. If they are different than 1.
    '''
    # plot the decision function
    xx, yy = np.meshgrid(np.linspace(0, math.ceil(max(X[:, 0])), 100), np.linspace(0, math.ceil(max(X[:, 1])), 100))
    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    vir_cmap = plt.cm.get_cmap('viridis')

    ###   Draws each decision level (if present) with color from the viridis colormap   ###
    if decision_levels is not '' and distance_or_identity == 'dist':
        draw_thresholds(axis, fig, X, xx, yy, Z, decision_levels)

    ###   Draws the decision function as a red line   ###
    else:
        axis.contour(xx, yy, Z, levels=[0],colors='r', linestyles=['-'], linewidths=0.5)

    abs_length = [float(n)**2 for n in sample_weight]
    if distance_or_identity == 'id':
        if decision_levels is not '':
            title +=" with tpr "+"{:2.2f}".format(decision_levels[next(iter(decision_levels))][0])\
                    +" and tnr "+"{:2.2f}".format(decision_levels[next(iter(decision_levels))][1])\
                    +" per alignment"
        mappedDict = dict(zip(sorted(set(aln_names)), range(len(aln_names))))
        aln_identities = list(map(lambda x: mappedDict[x], aln_names))
        scatter = axis.scatter(X[:, 0], X[:, 1], c=aln_identities, alpha=0.75, s=abs_length,
                 cmap=plt.cm.tab20, edgecolors='black')
        inv_map = {v: k for k, v in mappedDict.items()}
        aln_labels = list()
        for label in scatter.legend_elements()[1]:
            aln_labels.append(inv_map[int(re.findall(r'\d+', label)[0])])
        plt.legend(scatter.legend_elements()[0], aln_labels,
                    loc="upper right", title="Alignments")

    ###   Plot scatter of segments   ###
    if decision_levels is not '':
        scatter = axis.scatter(X[:, 0], X[:, 1], c=y, alpha=0.5, 
        s=abs_length, cmap=plt.cm.bone, edgecolors='black')
    else:
        import seaborn as sns
        from operator import itemgetter
        dummy_levels = dict()
        for thr in np.arange(-0.5,1,0.1):
            dummy_levels[thr]=0
        draw_thresholds(axis, fig, X, xx, yy, Z, dummy_levels, clean=True)
        label_order = []
        scatter = sns.scatterplot(X[:, 0], X[:, 1], hue=aln_names, s=abs_length, 
            palette="tab20", edgecolor='black')
        ###   Legend labels ordering   ###
        handles, labels = axis.get_legend_handles_labels()

        count_bellow_1=0
        for tup in sorted(label_order_tups, key = itemgetter(1)):
            if tup[1] < 0.05:
                count_bellow_1+=1
            label_order.append(tup[2])
        ordered_labels = [labels[idx] for idx in label_order]
        ordered_handles = [handles[idx] for idx in label_order]
        lgnd = plt.legend(ordered_handles[:count_bellow_1],
                            ordered_labels[:count_bellow_1], 
                            title="Alignment E value")

    #Size legend needs work
    # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    # size_labels = list()
    # for label in labels:
    #     size_labels.append(round(math.sqrt(int(re.findall(r'\d+', label)[0]))))
    # legend2 = axis.legend(handles, size_labels, loc="lower right", title="Segment length")
    
    
    plt.xlim(0, math.ceil(max(X[:, 0])))
    plt.ylim(0, math.ceil(max(X[:, 1])))
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

def read_features(features_path):
    with open(features_path) as f:
        data = json.load(f)
    return data[1], [data[0]["maxX"],data[0]["maxY"]]

def main(commandline_arguments):
    '''Main entry point'''
    comm_args, parser = create_and_parse_argument_options(commandline_arguments)

    ###   Load alignment segment data   ###
    csv_list = csv_iterator(comm_args.csv_path)
    if comm_args.top_segments:
        csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    if comm_args.center_mass_segments:
        csv_list = recalculate_data_by_averaging_segments(csv_list)

    ###   Test data   ###
    #max_features = read_max_features(comm_args.pickle+".maxvals")
    train_args, max_features = read_features(comm_args.pickle+".json")
    classifier = cPickle.load(open(comm_args.pickle, 'rb'))
    segment_pred_dist = test_function(csv_list, classifier, max_features)

    if comm_args.test_classifier_precision:
        grouped_data = flatten_alignment_segment_stats_to_groups(segment_pred_dist, by_group=True)
        dist_to_se_sp_pr = mass_test(segment_pred_dist, grouped_data, 
            comm_args.test_classifier_precision, min_threshold=-1, max_threshold=0.5, step=0.1)
        #plot_test_histograms(grouped_data)
    else:
        grouped_data = flatten_alignment_segment_stats_to_groups(segment_pred_dist)
        
        #In the case of large segments there will be few segments and they will be
        #far away from the boundary => E value nearing 0.
        #In the case of many small segments their distance to the boundary will be 
        #accumulated resulting in big negative number (larger than any segment can
        # attain on it's own) => E value nearing infinity.
        alnid_with_evalue = list()
        alnid_to_eval = dict()
        i=0
        for aln in sorted(grouped_data.keys()):
            eval = math.exp((grouped_data[aln][0][0]/grouped_data[aln][0][2])*-1)
            alnid_to_eval[aln] = eval
            alnid_with_evalue.append((aln, eval, i))
            i+=1


    ###   Plot the classifier   ###
    if comm_args.plot_df:
        plot_title = re.sub('.csv','',str(re.sub(r'.*/','',comm_args.csv_path))).replace("_", " ")
        X, y, sample_weight, maxX, maxY, aln_names = load_csv_data(csv_list, max_features=max_features)
        fig, axes = plt.subplots(1, 1, )
        if comm_args.test_classifier_precision:
            plot_decision_function(classifier,X,y, sample_weight, axes, fig,
                                plot_title, aln_names, 
                                decision_levels=dist_to_se_sp_pr, 
                                distance_or_identity=comm_args.test_classifier_precision)
        else:
            aln_names_eval = list()
            for alnid in aln_names:
                aln_names_eval.append(("{:2.1e}".format(alnid_to_eval[alnid])+" "+alnid[:15]))
            plot_decision_function(classifier,X,y, sample_weight, axes, fig,
                                plot_title,aln_names_eval,label_order_tups=alnid_with_evalue)
        plt.tight_layout()
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))