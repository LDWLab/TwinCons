#!/usr/bin/env python3
"""Calculates segments for multiple or single alignments"""
import re, os, sys, csv, argparse, matplotlib
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

import twincons.TwinCons as TwinCons

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    input_path = parser.add_mutually_exclusive_group(required=True)
    input_path.add_argument('-a', '--alignment_path', help='Path to folder with alignment files.')
    input_path.add_argument('-twc', '--twincons_path', help='Path to folder with csv output files from TwinCons.py')
    parser.add_argument('output_path', help='Path to image for output.')
    parser.add_argument('-t','--length_threshold', help='Threshold for consecutive low scores that split positive segments.\
                                                \nDefault: 3', type=int, default=3)
    parser.add_argument('-it','--intensity_threshold', help='Threshold for intensity over which a score is considered truly positive.\
                                                \nDefault: 1', type=float, default=1)
    parser.add_argument('-avew','--average_weight', help='Use average weight for segments, instead of using their total weight.', action="store_true", default=False)
    parser.add_argument('-np','--treat_highly_negative_as_conserved', help='Treat low scoring positions as conserved for segment calculation. \
                                                \nConsiders the absolute for negative positions when comparing with intensity threshold.', action="store_true", default=False)
    parser.add_argument('-c','--csv', help='Output length and weight distributions in a csv file. \
                                                \nUses the output file name specified by appending .csv', action="store_true")
    parser.add_argument('-p','--plot', help='Plot a scatter of the segments.', default=False, action="store_true")
    parser.add_argument('-l','--legend', help='Draw a legend.', default=False, action="store_true")
    parser.add_argument('-co', '--calculation_options', help='Options for TwinCons calculation. See README for details.', nargs='+')
    #parser.add_argument('-s','--structure_path', help='Path to folder with structure files; names should match alignment groups within files.')
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
        if len(segment_lengths) == 0:
            segment_lengths.append(0)
        if len(segment_weights) == 0:
            segment_weights.append(0)
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

def split_by_thresholds(score_list, intensity_thr, length_low, treat_highly_negative_as_conserved=False):
    '''Takes in list with alignment positions and scores 
    returns a list of list with segments split by thresholds
    '''
    low_count, high_count = 0, 0
    split_by = list()
    for i in range(len(score_list)):
        if treat_highly_negative_as_conserved:
            test_score = abs(score_list[i][1])
        else:
            test_score = score_list[i][1]
        if test_score < intensity_thr:
            low_count += 1
        if test_score >= intensity_thr:
            high_count += 1
            low_count = 0
        if low_count >= length_low and high_count != 0:
            split_by.append(i)
            high_count, low_count = 0, 0
    segments = [score_list[i : j] for i, j in zip([0] + split_by, split_by + [None])] 
    return segments

def calculate_segment_stats(unfiltered_segment_list, thr, aln_length, treat_highly_negative_as_conserved=False):
    segment_stats = list()
    for unf_segment in unfiltered_segment_list:
        weight = 0
        res = [idx for idx, val in enumerate(unf_segment) if val[1] > thr]
        if len(res) == 0:
            continue
        segment = unf_segment[res[0]:res[len(res)-1]+1]
        for pos in segment:
            if treat_highly_negative_as_conserved:
                pos_score = abs(pos[1])
            else:
                pos_score = pos[1]
            if pos_score > thr:
               weight += pos_score
        segment_stats.append((len(segment), len(segment)/aln_length, weight, (segment[0][0],segment[len(segment)-1][0])))
    return segment_stats

def parse_arguments_for_twc(calculation_options, list_for_twincons):
    for opt in calculation_options:
        if re.search(r'_', opt):
            list_for_twincons.append('-'+opt.split('_')[0])
            list_for_twincons.append(opt.split('_')[1])
        else:
            list_for_twincons.append('-'+opt)
    return list_for_twincons

def run_twincons(file_dir, calculation_options):
    ###   Constructing arguments for TwinCons   ###
    list_for_twincons = ['-a',file_dir, '-r']
    list_for_twincons = parse_arguments_for_twc(calculation_options, list_for_twincons)
    print(list_for_twincons)
    
    ###   Executing TwinCons   ###
    return TwinCons.main(list_for_twincons)

def load_twincons_csv(file_dir):
    first_line = True
    alnindex_score = dict()
    with open(file_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if first_line:
                first_line = False
                continue
            alnindex_score[row[0]] = float(row[1])
    return alnindex_score, len(alnindex_score), dict()

def calc_segments_for_aln(alnindex_score, num_alned_pos, int_thr, length_thr, highly_neg_as_pos=False, path_type='aln_path'):
    '''Calculating segment stats'''
    out_dict = dict()
    score_list = list()
    for x in alnindex_score.keys():
        position = x
        if path_type == 'aln_path':
            pos_score = alnindex_score[x][0]
        elif path_type == 'twc_path':
            pos_score = alnindex_score[x]
        else:
            raise IOError("Missing input for TwinCons data!")
        score_list.append((position, pos_score))
        out_dict[position] = pos_score
    unfiltered_segment_list = split_by_thresholds(score_list, 
                                                  int_thr, 
                                                  length_thr, 
                                                  treat_highly_negative_as_conserved=highly_neg_as_pos)
    segment_stats = calculate_segment_stats(unfiltered_segment_list, 
                                            int_thr, 
                                            num_alned_pos, 
                                            treat_highly_negative_as_conserved=highly_neg_as_pos)
    return segment_stats

def main(commandline_args):
    comm_args = create_and_parse_argument_options(commandline_args)
    alns_to_segment_stats = dict()
    if comm_args.alignment_path:
        if not comm_args.calculation_options:
            raise IOError("Must specify twincons options when running on alignment files!")
        file_dir = comm_args.alignment_path
        regex = r'(.*\/)(.*)(\.fasta|\.fas|\.fa)'
    elif comm_args.twincons_path:
        file_dir = comm_args.twincons_path
        regex = r'(.*\/)(.*)(\.csv)'
    else:
        raise IOError("Missing input for TwinCons data!")
    
    for file in os.listdir(file_dir):
        if not re.findall(regex, file_dir+file):
            continue
        if comm_args.alignment_path:
            path_type = 'aln_path'
            try:
                alnindex_score, sliced_alns, number_of_aligned_positions, gap_mapping = run_twincons(file_dir+file, comm_args.calculation_options)
            except Exception:
                raise Exception("TwinCons failed to run!")
        elif comm_args.twincons_path:
            path_type = 'twc_path'
            alnindex_score, number_of_aligned_positions, gap_mapping = load_twincons_csv(file_dir+file)
        else:
            raise IOError("Missing input for TwinCons data!")
        segment_stats = calc_segments_for_aln(alnindex_score, 
                                            number_of_aligned_positions, 
                                            comm_args.intensity_threshold,
                                            comm_args.length_threshold,
                                            highly_neg_as_pos=comm_args.treat_highly_negative_as_conserved,
                                            path_type=path_type)
        alns_to_segment_stats[file.split('.')[0]] = segment_stats
    if comm_args.plot:
        fig, ax, plt = scatter_plot(alns_to_segment_stats, legend=comm_args.legend, average_weight=comm_args.average_weight)
        if comm_args.output_path.endswith('.png'):
            plt.savefig(comm_args.output_path, dpi=600, bbox_inches='tight')
        else:
            plt.savefig(comm_args.output_path+'.png', dpi=600, bbox_inches='tight')
    if comm_args.csv:
        csv_output(comm_args.output_path, alns_to_segment_stats)

if __name__ == "__main__":
    main(sys.argv[1:])