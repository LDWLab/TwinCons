#!/usr/bin/env python3
"""Calculates segments for multiple or single alignments"""
import TwinCons
import re, os, sys, csv, getopt, argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
#from labellines import labelLine, labelLines
from operator import itemgetter

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_path', help='Path to folder with alignment files.')
    parser.add_argument('output_path', help='Path to image for output.')
    parser.add_argument('-co', '--calculation_options', help='Options for score calculation. See README for details.', required=True, nargs='+')
    parser.add_argument('-t','--length_threshold', help='Threshold for number of allowed bad scores when calculating length of positive sections.', type=int, default=1, required=False)
    parser.add_argument('-avew','--average_weight', help='Instead of plotting total weight per segment, plot average weight', action="store_true", required=False, default=False)
    parser.add_argument('-it','--intensity_threshold', help='Threshold for intensity over which a score is considered truly positive.', type=float, default=1, required=False)
    parser.add_argument('-s','--structure_path', help='Path to folder with structure files; names should match alignment groups within files.')
    parser.add_argument('-c','--csv', help='Output length and weight distributions in a csv file. Uses the output file name specified by appending .csv', required=False, action="store_true")
    parser.add_argument('-l','--legend', help='Draw a legend', default=False, action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def scatter_plot(alns_to_segment_stats, legend=False, average_weight=False):
    '''Outputs scatter plot image with colors depending on the number of the input alignments.
    '''
    ###   Defining Color scheme   ###
    fig, ax = plt.subplots(1, 1, )
    if len(alns_to_segment_stats) == 20:
        colors = matplotlib.cm.tab20(np.linspace(0, 1, len(alns_to_segment_stats)))
    elif len(alns_to_segment_stats) == 10:
        #In cases of just 10 alignments it assumes two sets of 5 each and uses a seismic/divergent gradient.
        #colors = matplotlib.cm.seismic(np.linspace(0, 1, len(weight_distr)))
        colors = matplotlib.cm.PRGn(np.linspace(0, 1, len(alns_to_segment_stats)))
    elif len(alns_to_segment_stats) > 21 and re.match('A_', sorted(alns_to_segment_stats.keys())[0]):
        #Creates a color mapping for groups of alignments defined with prepended A_, B_, C_ and so on.
        #Used in cases where there are too many alignments to properly discriminate with colors.
        colors=[]
        letters=[]
        for letter in sorted(alns_to_segment_stats.keys()):
            letters.append(letter.split('_')[0])

        color_set = matplotlib.cm.viridis(np.linspace(0, 1, len(set(letters))))
        letter_color_map = dict(zip(sorted(set(letters)),color_set))

        for letter in letters:
            colors.append(letter_color_map[letter])
    else:
        colors = matplotlib.cm.tab20(np.linspace(0, 1, len(alns_to_segment_stats)))

    ###   Plotting   ###
    sorted_names = sorted(alns_to_segment_stats.keys())
    file_number = 0
    label_order_tups,label_order = [], []
    for file, color in zip(sorted_names,colors):
        scaled_lengths = [n[1] for n in alns_to_segment_stats[file]]
        if average_weight:
            segment_weights = [n[2]/n[0] for n in alns_to_segment_stats[file]]
        else:
            segment_weights = [n[2] for n in alns_to_segment_stats[file]]
        segment_lengths = [n[0] for n in alns_to_segment_stats[file]]
        abs_length = [n**2 for n in segment_lengths]
        plt.scatter(scaled_lengths, segment_weights, label=re.sub(r'\.fas.*','',file),marker='.',s=abs_length,color=color)
        label_order_tups.append((file_number,max(segment_weights)+max(segment_lengths)))
        file_number += 1
    
    ###   Legend labels ordering   ###
    handles, labels = ax.get_legend_handles_labels()
    for tup in sorted(label_order_tups, key = itemgetter(1), reverse=True):
        label_order.append(tup[0])
    ordered_labels = [labels[idx] for idx in label_order]
    ordered_handles = [handles[idx] for idx in label_order]

    ###   Legend block   ###
    if legend:
        lgnd = plt.legend(ordered_handles,ordered_labels,bbox_to_anchor=(1.04,1), borderaxespad=0)
        for n in range(len(alns_to_segment_stats)):
            lgnd.legendHandles[n]._sizes = [30]
    return fig, ax, plt

def csv_output(output_path, alns_to_segment_stats):
    '''Writes out data used for generating the plot in a csv file.
    '''
    if output_path.endswith('.png'):
        file_for_writing = re.sub(r'.png','.csv',output_path)
    elif output_path.endswith('.csv'):
        file_for_writing = output_path
    else:
        file_for_writing = output_path+'.csv'

    with open(file_for_writing, mode ="w") as output_csv:
        csv_writer = csv.writer(output_csv)
        csv_writer.writerow(['File name', 'Segment length', 'Normalized segment length','Total segment weight', 'Alignment position'])
        for file, segments in sorted(alns_to_segment_stats.items()):
            for one_segment in sorted(segments, key=lambda x: x[0]):
                aln_pos = str(one_segment[3][0])+'-'+str(one_segment[3][1])
                #print([file, one_segment[0], one_segment[1],one_segment[2], one_segment[3]])
                csv_writer.writerow([file, one_segment[0], one_segment[1],one_segment[2], aln_pos])
    return True

def split_by_thresholds(score_list, intensity_thr, length_low):
    '''Takes in list with alignment positions and scores 
    returns a list of list with segments split by thresholds
    '''
    low_count, high_count = 0, 0
    split_by = list()
    for i in range(len(score_list)):
        if score_list[i][1] < intensity_thr:
            low_count += 1
        if score_list[i][1] >= intensity_thr:
            high_count += 1
            low_count = 0
        if low_count >= length_low and high_count != 0:
            split_by.append(i)
            high_count, low_count = 0, 0
    segments = [score_list[i : j] for i, j in zip([0] + split_by, split_by + [None])] 
    return segments

def calculate_segment_stats(unfiltered_segment_list, thr, aln_length):
    segment_stats = list()
    for unf_segment in unfiltered_segment_list:
        weight = 0
        res = [idx for idx, val in enumerate(unf_segment) if val[1] > thr]
        if len(res) == 0:
            continue
        segment = unf_segment[res[0]:res[len(res)-1]+1]
        for pos in segment:
           if pos[1] > thr:
               weight += pos[1]
        segment_stats.append((len(segment), len(segment)/aln_length, weight, (segment[0][0],segment[len(segment)-1][0])))
    return segment_stats

def main(commandline_args):
    comm_args = create_and_parse_argument_options(commandline_args)
    alns_to_segment_stats = dict()
    for file in os.listdir(comm_args.alignment_path):
        if not re.findall(r'(.*\/)(.*)(\.fasta|\.fas|\.fa)',comm_args.alignment_path+file):
            print("Skipping non-alignment file "+file)
            continue
        
        ###   Constructing arguments for TwinCons   ###
        out_dict={}
        list_for_twincons = ['-a',comm_args.alignment_path+file, '-r']
        for opt in comm_args.calculation_options:
            if re.search(r'_', opt):
                list_for_twincons.append('-'+opt.split('_')[0])
                list_for_twincons.append(opt.split('_')[1])
            else:
                list_for_twincons.append('-'+opt)
        print(list_for_twincons)
        
        ###   Executing TwinCons   ###
        alnindex_score, sliced_alns, number_of_aligned_positions, gap_mapping = TwinCons.main(list_for_twincons)
        inv_map = dict()
        for k, v in gap_mapping.items():
            if v in inv_map.keys():
                continue
            inv_map[v] = k
        
        ###   Calculating segment stats  ###
        score_list = list()
        for x in alnindex_score.keys():
            pos = x
            if len(inv_map) > 0:
                pos = inv_map[x]
            score_list.append((pos,alnindex_score[x][0]))
            out_dict[pos] = alnindex_score[x][0]
        unfiltered_segment_list = split_by_thresholds(score_list, comm_args.intensity_threshold, comm_args.length_threshold+1)
        segment_stats = calculate_segment_stats(unfiltered_segment_list, comm_args.intensity_threshold, number_of_aligned_positions)
        alns_to_segment_stats[file.split('.')[0]] = segment_stats

    fig, ax, plt = scatter_plot(alns_to_segment_stats, legend=comm_args.legend, average_weight=comm_args.average_weight)
    if comm_args.output_path.endswith('.png'):
        plt.savefig(comm_args.output_path, dpi=600, bbox_inches='tight')
    else:
        plt.savefig(comm_args.output_path+'.png', dpi=600, bbox_inches='tight')
    if comm_args.csv:
        csv_output(comm_args.output_path, alns_to_segment_stats)

if __name__ == "__main__":
    main(sys.argv[1:])