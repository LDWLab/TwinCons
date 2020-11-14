#!/usr/bin/env python3
"""
Evaluates alignment entries in csv generated from MultiPhyMeas.
Requires a decision function and max features generated from SVM_train.
Train and test only with the same parameters!
Such parameters can be % cutting gaps, center mass segments, top segments.
"""
import re, sys, csv, math, argparse, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from operator import itemgetter
import _pickle as cPickle

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data', type=str)
    parser.add_argument('output_path', help='Path to output significant segment results', type=str)
    parser.add_argument('pickle',help='Provide path to classifier pickle binary file. The script will also search for an identically\
                                    \nnamed file with extension ".maxvals" containing comma separated maximal values, for example:\
                                    \npickle file:\trandom_test.pkl\
                                    \nmaximal values:\trandom_test.pkl.maxvals', type=str)
    parser.add_argument('-pd','--plot_df', help='Path to output plot for the decision function.', type=str)
    parser.add_argument('-et','--evalue_threshold', help='Use confidence calculation for segment determination. Uses the range distance thresholds Start.',
                                                                            action="store_true")
    parser.add_argument('-ts','--top_segments', help='Limit input for each alignment to the top segments that cover\
                                    \n this percentage of the total normalized length and weight. (Default = 0.5)', 
                                                                                            type=float, default=0.5)
    parser.add_argument('-l','--length_type_calculation', help='Choose what type of segment calculation should be used.\
        \n\t absolute:   absolute length of the segments.\
        \n\t normalized: length of segments is normalized with the total alignment length.\
        \n\t cms:        average position (center of mass) from all segments per alignment.', choices=['absolute', 'normalized', 'cms'], default='normalized')
    parser.add_argument('-cms','--center_mass_segments', help='Use the average position (Center of mass) \
                                    \nfrom all segments per alignment.', action="store_true")
    calculate_positive = parser.add_mutually_exclusive_group(required=True)
    calculate_positive.add_argument('-tcp','--test_classifier_precision', help='Provided csv is annotated for testing the classifier.', action="store_true")
    calculate_positive.add_argument('-tqa','--test_query_alignments', help='Provided csv is a query and is not annotated for testing the classifier.', action="store_true")
    parser.add_argument('-dt', '--range_distance_thresholds', nargs=3, metavar=('Start', 'End', 'Step'), 
                                    help='Range of distances from the decision boundary to evaluate. Also works with --evalue_threshold\
                                    \nDefault for non evalue (-20, 20, 0.05).\
                                    \nDefault for evalue distance (0, 1, 0.0001)', default=[-20, 20, 0.05], type = float)
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
        aln_id_to_data[row[0]].append((float(row[2]), float(row[3]), float(row[1]), str(row[4])))
    return aln_id_to_data

def trim_data_by_top_segments(csv_list, perc_top_segments):
    '''Returns trimmed down list by percentage coverage of the 
    total segment weights and total normalized segment lengths '''
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
        averaged_results = average([x[:-1] for x in data_items])
        output_list.append([aln_id, averaged_results[2], averaged_results[0], averaged_results[1], 'NA'])
    return output_list

def use_absolute_length_of_segments(csv_list):
    '''Moves absolute segment length to accurate position.
    Assigns weights of 1 to all segments.'''
    return [[x[0], 1, x[1], x[3], x[4]] for x in csv_list]

def load_and_assign_data(csv_list):
    '''Reads a csv outputted from CalculateSegments.py, 
    entries starting with A_ and B_ are considered as +
    data and entries starting with C_ and D_ as - data.
    '''
    data_xy, data_identity, data_weights, aln_names = list(), list(), list(), list()
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
    return data_xy, data_identity, data_weights, aln_names

def load_csv_data(csv_list, max_features=''):
    '''
    Normalizes all the data for x and y in the range 0,1
    using the max(x) and max(y).
    '''
    data_xy, data_identity, data_weights, aln_names = load_and_assign_data(csv_list)
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
            dist_sum += segm[1][1]
            if segm[1][0] == 1:
                pos_sum += 1
        grouped_data[outname].append((dist_sum, pos_sum, tot_segm))
    return grouped_data

def bypass_zero_division(x,y):
    if y != 0:
        return x/y
    else:
        return 0

def compare_thr(thr, segments, eval=False):
    negative, positive = 0, 0
    for i in segments:
        if eval:
            positive += (math.exp(i[1][1]*i[1][2]*-1) < thr)
            negative += (math.exp(i[1][1]*i[1][2]*-1) >= thr)
        else:
            positive += (i[1][1] > thr)
            negative += (i[1][1] <= thr)
    return (negative, positive)

def calc_identity_stats(segment_pred_dist, thr, tp, tn, fp, fn, eval=False):
    for aln, segments in segment_pred_dist.items():
        if re.match('^A_|^B_', aln) is not None:
            if compare_thr(thr, segments, eval=eval)[1] == 0:
                fn += 1
            if compare_thr(thr, segments, eval=eval)[1] > 0:
                tp += 1
        if re.match('^C_|^D_', aln) is not None:
            if compare_thr(thr, segments, eval=eval)[1] == 0:
                tn += 1
            if compare_thr(thr, segments, eval=eval)[1] > 0:
                fp += 1
    #print(tp, tn, fp, fn)
    return tp, tn, fp, fn

def mass_test(segment_pred_dist, min_threshold=0, max_threshold=2, step=0.1, eval=False):
    '''For evaluating the decision function with 
    great number of alignments in separate groups.
    '''
    ###   Calculate Sensitivity(tpr), Specificity(tnr) and Precision   ###
    ###   over a range of distances from the decision boundary         ###
    dist_to_stats = {}
    for thr in np.arange(min_threshold,max_threshold,step):
        tp, tn, fp, fn = 0,0,0,0
        tp, tn, fp, fn = calc_identity_stats(segment_pred_dist, thr, tp, tn, fp, fn, eval=eval)
        tpr = bypass_zero_division(tp,tp+fn)
        tnr = bypass_zero_division(tn,tn+fp)
        precision = bypass_zero_division(tp,tp+fp)
        dist_to_stats[thr] = (tpr, tnr, precision)
        #print ("Threshold "+str(thr),'\n',"tpr", tpr,"tnr", tnr,'\n',"precision",precision)
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
    
    ###   Draws a line for each decision level   ###
    if math.ceil(max(X[:, 1])) == 1 and math.ceil(max(X[:, 0])) == 1:
        CS2 = axis.contour(xx, yy, Z, colors='magenta', levels=sorted(decision_levels.keys()),
                          alpha=0.75, linestyles=[':'], linewidths=0.5)
        axis.clabel(CS2, fmt='%2.1f', colors='magenta', fontsize=5)
    if clean:
        cbar.set_ticklabels(["%.1f" %thr for thr in sorted(decision_levels.keys())])
        cbar.ax.set_title('Distance')
    else:
        levels_labels = ["%.2f, %.2f, %.2f" %decision_levels[thr] for thr in sorted(decision_levels.keys())]
        cbar.set_ticklabels(["%.1f %s" % (thr,lev) for lev,thr in zip(levels_labels, sorted(decision_levels.keys()))])
        cbar.ax.set_title('Dist  tpr, tnr, prec', x=5)
    return True

def plot_decision_function(classifier, X, y, sample_weight, axis, fig, title, aln_names, label_order_tups=None, cms_or_identity=False, decision_levels='', thresholds=''):
    '''
    Plots the data and the decision function. Besides the classifier function, takes in the sample weights
    and plots it as a size of the datapoints. If they are different than 1.
    '''
    # plot the decision function
    xx, yy = np.meshgrid(np.linspace(0, math.ceil(max(X[:, 0])), 100), np.linspace(0, math.ceil(max(X[:, 1])), 100))
    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)


    ###   Draws each decision level (if present) with color from the viridis colormap   ###
    if decision_levels is not '':
        draw_thresholds(axis, fig, X, xx, yy, Z, decision_levels)

    ###   Draws the decision function as a red line   ###
    #else:
    axis.contour(xx, yy, Z, levels=[0],colors='r', linestyles=['-'], linewidths=0.5)

    abs_length = [float(n)**2 for n in sample_weight]
    if abs_length == sample_weight:
        abs_length = [float(n)*20 for n in sample_weight]
        edgecolor = None
    else:
        edgecolor = "black"

    ###   Plot scatter of segments   ###
    if decision_levels is not '':
        scatter = axis.scatter(X[:, 0], X[:, 1], c=y, alpha=0.5, 
        s=abs_length, cmap=plt.cm.bone, edgecolor=edgecolor)
    else:
        import seaborn as sns
        from operator import itemgetter
        dummy_levels = dict()
        for thr in np.arange(float(thresholds[0]),float(thresholds[1])+float(thresholds[2]),float(thresholds[2])):
            dummy_levels[thr]=0
        draw_thresholds(axis, fig, X, xx, yy, Z, dummy_levels, clean=True)
        label_order = []
        scatter = sns.scatterplot(X[:, 0], X[:, 1], hue=aln_names, 
                palette="tab20", edgecolor=edgecolor, s=abs_length)

        ###   Legend labels ordering   ###
        handles, labels = axis.get_legend_handles_labels()

        count_bellow_1=0
        if cms_or_identity:
            for tup in sorted(label_order_tups, key = itemgetter(1)):
                if tup[1] < 1:
                    count_bellow_1+=1
                label_order.append(tup[2])
        else:
            for tup in sorted(label_order_tups, key = itemgetter(3), reverse=True):
                if tup[1] > 0:
                    count_bellow_1+=1
                label_order.append(tup[2])
        ordered_labels = [labels[idx] for idx in label_order]
        ordered_handles = [handles[idx] for idx in label_order]
        lgnd = plt.legend(ordered_handles[:count_bellow_1],
                            ordered_labels[:count_bellow_1], 
                            title=f"Alignments with segments\nabove threshold {thresholds[1]}")

    #Size legend needs work
    # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    # size_labels = list()
    # for label in labels:
    #     size_labels.append(round(math.sqrt(int(re.findall(r'\d+', label)[0]))))
    # legend2 = axis.legend(handles, size_labels, loc="lower right", title="Segment length")
    
    
    plt.xlim(0, math.ceil(max(X[:, 0])))
    plt.ylim(0, math.ceil(max(X[:, 1])))
    axis.set_title(title)
    return True

def read_features(features_path):
    with open(features_path) as f:
        data = json.load(f)
    return data[1], [data[0]["maxX"],data[0]["maxY"]]

def write_aln_rows(segments, csv_writer, aln):
    for segment in segments:
        if segment[1][0] != 0:
            csv_writer.writerow([aln, segment[1][1], segment[0], math.exp(segment[1][1]*segment[1][2]*-1)])
    return True

def main(commandline_arguments):
    '''Main entry point'''
    comm_args, parser = create_and_parse_argument_options(commandline_arguments)

    ###   Load alignment segment data   ###
    csv_list = csv_iterator(comm_args.csv_path)
    if comm_args.top_segments:
        csv_list = trim_data_by_top_segments(csv_list, comm_args.top_segments)
    if comm_args.length_type_calculation == 'cms':
        csv_list = recalculate_data_by_averaging_segments(csv_list)
    if comm_args.length_type_calculation == 'absolute':
        csv_list = use_absolute_length_of_segments(csv_list)

    if comm_args.evalue_threshold and (comm_args.range_distance_thresholds[0] == -20 or comm_args.range_distance_thresholds[1] == 20):
        comm_args.range_distance_thresholds[0] = 0
        comm_args.range_distance_thresholds[1] = 2
        comm_args.range_distance_thresholds[2] = 0.01

    ###   Test data   ###
    train_args, max_features = read_features(comm_args.pickle+".json")
    classifier = cPickle.load(open(comm_args.pickle, 'rb'))
    segment_pred_dist = test_function(csv_list, classifier, max_features)

    if comm_args.test_classifier_precision:
        dist_to_se_sp_pr = mass_test(segment_pred_dist,
            min_threshold=float(comm_args.range_distance_thresholds[0]), 
            max_threshold=float(comm_args.range_distance_thresholds[1])+float(comm_args.range_distance_thresholds[2]), 
            step=float(comm_args.range_distance_thresholds[2]), eval=comm_args.evalue_threshold)
        levels_labels_csv = [dist_to_se_sp_pr[thr] for thr in sorted(dist_to_se_sp_pr.keys())]
        roc_stats = [(thr,lev) for lev,thr in zip(levels_labels_csv, sorted(dist_to_se_sp_pr.keys()))]
        with open(comm_args.output_path, mode ="w") as output_csv:
            csv_writer = csv.writer(output_csv)
            csv_writer.writerow(['Alignment name', 'Distance from boundary', 'Alignment position', 'E-value'])
            for aln, segments in segment_pred_dist.items():
                write_aln_rows(segments, csv_writer, aln)
            csv_writer.writerow(['ROC stats'])
            csv_writer.writerow(['Distance from boundary', 'TPR', 'TNR', 'PREC'])
            for row in roc_stats:
                csv_writer.writerow([row[0],row[1][0],row[1][1],row[1][2]])
    if comm_args.test_query_alignments:
        alnid_with_evalue = list()
        alnid_to_eval = dict()
        i = 0
        if comm_args.length_type_calculation == 'cms':
            grouped_data = flatten_alignment_segment_stats_to_groups(segment_pred_dist)
            for aln in sorted(grouped_data.keys()):
                eval = math.exp(grouped_data[aln][0][0]*-1)
                alnid_to_eval[aln] = eval
                alnid_with_evalue.append((aln, eval, i))
                i+=1
        with open(comm_args.output_path, mode ="w") as output_csv:
            csv_writer = csv.writer(output_csv)
            csv_writer.writerow(['Alignment name', 'Distance from boundary', 'Alignment position', 'E-value'])
            for aln, segments in segment_pred_dist.items():
                segment_id = compare_thr(float(comm_args.range_distance_thresholds[1]), segments)
                alnid_with_evalue.append((aln, segment_id[1], i, max(segments, key = itemgetter(1))[1][1] ))
                alnid_to_eval[aln] = segment_id[1]
                i+=1
                write_aln_rows(segments, csv_writer, aln)

    ###   Plot the classifier   ###
    if comm_args.plot_df:
        plot_title = re.sub('.csv','',str(re.sub(r'.*/','',comm_args.csv_path))).replace("_", " ")
        X, y, sample_weight, maxX, maxY, aln_names = load_csv_data(csv_list, max_features=max_features)
        fig, axes = plt.subplots(1, 1, )
        if comm_args.test_classifier_precision:
            plot_decision_function(classifier,X,y, sample_weight, axes, fig,
                                plot_title, aln_names, 
                                decision_levels=dist_to_se_sp_pr, 
                                cms_or_identity=comm_args.center_mass_segments)
        elif comm_args.test_query_alignments:
            aln_names_eval = list()
            for alnid in aln_names:
                if comm_args.center_mass_segments:
                    aln_names_eval.append(("{:2.1e}".format(alnid_to_eval[alnid])+" "+alnid[:15]))
                else:
                    aln_names_eval.append((str(alnid_to_eval[alnid])+" "+alnid[:15]))
            plot_decision_function(classifier,X,y, sample_weight, axes, fig,
                                plot_title,aln_names_eval,label_order_tups=alnid_with_evalue,
                                cms_or_identity=comm_args.center_mass_segments,
                                thresholds=comm_args.range_distance_thresholds)
        plt.tight_layout()
        plt.savefig(comm_args.plot_df, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))