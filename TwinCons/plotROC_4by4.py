#!/usr/bin/env python3
import re, os, csv, sys, pprint
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc

def read_csv(csv_location):
    '''Put csv in list'''
    output_list = list()
    segment_data = True
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if segment_data:
                if row[0] == 'ROC stats':
                    segment_data = False
                continue
            output_list.append(row)
    return output_list[1:]

def construct_param_struc(split_label, datastruc, tpr, fpr):
    
    if (split_label[1],split_label[6]) not in datastruc.keys():
        datastruc[(split_label[1],split_label[6])] = dict()

    if (split_label[2],split_label[5]) not in datastruc[(split_label[1],split_label[6])].keys():
        datastruc[(split_label[1],split_label[6])][(split_label[2],split_label[5])] = list()

    datastruc[(split_label[1],split_label[6])][(split_label[2],split_label[5])].append(((split_label[3],split_label[7],split_label[9],split_label[10]), (tpr, fpr)))
    return datastruc

def plot_five_by_two(datastruc, file_name, titleplot):
    rows = len(datastruc) 
    cols = list()
    for key in datastruc.keys():
        cols.append(len(datastruc[key]))
    fig, axs = plt.subplots(nrows=rows, ncols=max(cols),
                            sharey=True, sharex=True, 
                            figsize=(16, 24))
    fig.suptitle(titleplot, fontsize=16)
    for row, lt in zip(list(range(rows)), datastruc.keys()):
        #print (row, lt, end=" ")
        for col, it in zip(list(range(max(cols))), datastruc[lt].keys()):
            #print(col, it)
            subplot_title = f'Matrix: {lt[0]}, Length thr: {lt[1]},\n\
Weighted: {it[0]}, Intensity thr: {it[1]}'
            #print(subplot_title)
            axs[row,col].set_title(subplot_title)
            datatruc_idx = range(0, len(datastruc[lt][it]))
            colormap_ix = [item for sublist in [[x]*8 for x in [0.2,0.5,0.8]] for item in sublist]
            markers = ['.','.','1','1']*6
            colormaps = [plt.cm.Purples, plt.cm.Greens]*12
            linestyles = ['-','-','-','-',':',':',':',':']*3
            for label_to_data, i in zip(datastruc[lt][it], datatruc_idx):
                accuracy = auc(label_to_data[1][1],label_to_data[1][0])
                gap_label = "{:.0%} ".format(float(label_to_data[0][0].replace('p','.')[2:]))
                ts_label = label_to_data[0][2].replace('p','.')
                legend_label = gap_label+f"{label_to_data[0][1]} {ts_label} {label_to_data[0][3]} "+"{:.1%}".format(float(accuracy))
                axs[row,col].plot(label_to_data[1][1], 
                                  label_to_data[1][0], 
                                  label=legend_label, 
                                  color=colormaps[i](colormap_ix[i]),
                                  linestyle=linestyles[i],
                                  marker=markers[i],
                                  linewidth=1, markersize=4)
                axs[row,col].legend(loc=4, title="% cut gaps   AUC", prop={'size': 6})
                axs[row,col].set_xlabel("False positive rate")
                axs[row,col].set_ylabel("True positive rate")
    plt.savefig(file_name, dpi=600)

directory = "./data/test_twc_parameters/out_stats/eval/"
bbs,rprot,indeli = dict(),dict(),dict()
for file in os.listdir(directory):
    if not re.findall(r'(.*)(\.csv)',file):
        #print("Skipping non-csv file "+file)
        continue
    csv_data = read_csv(directory+file)
    tpr = [float(x[1]) for x in csv_data]
    fpr = [1-float(x[2]) for x in csv_data]
    split_label = file.replace('.csv','').split("_")
    #print(split_label)
    if split_label[0] == "BBSvBBS":
        bbs = construct_param_struc(split_label, bbs, tpr, fpr)
    if split_label[0] == "BBSvrProt":
        rprot = construct_param_struc(split_label, rprot, tpr, fpr)
    if split_label[0] == "BBSvIND":
        indeli = construct_param_struc(split_label, indeli, tpr, fpr)

plot_five_by_two(bbs, "./data/outputs/PNG/param_test/BBS-BBS_ev.png","BBS vs BBS")
plot_five_by_two(rprot, "./data/outputs/PNG/param_test/BBS-rProt_ev.png","BBS vs rProt")
plot_five_by_two(indeli, "./data/outputs/PNG/param_test/BBS-INDELI_ev.png","BBS vs INDELI")

