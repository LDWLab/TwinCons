#!/usr/bin/env python3
import re, os, csv, sys
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
    if split_label[2] not in datastruc.keys():
        datastruc[split_label[2]] = dict()
    if split_label[3] not in datastruc[split_label[2]].keys():
        datastruc[split_label[2]][split_label[3]] = list()
    datastruc[split_label[2]][split_label[3]].append((split_label[1], (tpr, fpr)))
    return datastruc

def plot_five_by_two(datastruc, file_name, titleplot):
    rows = len(datastruc) 
    cols = list()
    for key in datastruc.keys():
        cols.append(len(datastruc[key]))
    fig, axs = plt.subplots(nrows=rows, ncols=max(cols),
                            sharey=True, sharex=True, 
                            figsize=(20, 12))
    fig.suptitle(titleplot, fontsize=16)
    for row, lt in zip(list(range(rows)), datastruc.keys()):
        #print (row, lt)
        for col, it in zip(list(range(max(cols))), datastruc[lt].keys()):
            subplot_title = "Segments split by "+str(int(lt[2])+1)+" low positions\n\
TWC low position threshold: "+it.replace('p','.')[2:]
            axs[row,col].set_title(subplot_title)
            color_idx = np.linspace(0, 1, len(datastruc[lt][it]))
            for label_to_data, i in zip(datastruc[lt][it], color_idx):
                accuracy = auc(label_to_data[1][1],label_to_data[1][0])
                legend_label = "{:.0%}".format(float(label_to_data[0].replace('p','.')[2:]))+"    {:.1%}".format(float(accuracy))
                axs[row,col].plot(label_to_data[1][1], 
                                  label_to_data[1][0], 
                                  label=legend_label, 
                                  color=plt.cm.viridis(i))
                axs[row,col].legend(loc=4, title="% cut gaps   AUC")
                axs[row,col].set_xlabel("False positive rate")
                axs[row,col].set_ylabel("True positive rate")
    plt.savefig(file_name, dpi=600)

directory = "./data/test_twc_parameters/out_stats/BL_unw/SeparateSegments/tcp_id_C1_ts1/"
bbs,rprot,indeli = dict(),dict(),dict()
for file in os.listdir(directory):
    if not re.findall(r'(.*)(\.csv)',file):
        #print("Skipping non-csv file "+file)
        continue
    csv_data = read_csv(directory+file)
    tpr = [float(x[1]) for x in csv_data]
    fpr = [1-float(x[2]) for x in csv_data]
    split_label = file.replace('.csv','').split("_")
    label = split_label[1]+"-"+split_label[2]+"-"+split_label[3]
    if split_label[0] == "BBSMlvBBSl":
        bbs = construct_param_struc(split_label, bbs, tpr, fpr)
    if split_label[0] == "BBSlvrProt":
        rprot = construct_param_struc(split_label, rprot, tpr, fpr)
    if split_label[0] == "BBSlvIND":
        indeli = construct_param_struc(split_label, indeli, tpr, fpr)

plot_five_by_two(bbs, "./data/outputs/PNG/param_test/BBS-BBS_C2_ts1.png","BBS vs BBS all segments penalty 2")
plot_five_by_two(rprot, "./data/outputs/PNG/param_test/BBS-rProt_C2_ts1.png","BBS vs rProt all segments penalty 2")
plot_five_by_two(indeli, "./data/outputs/PNG/param_test/BBS-INDELI_C2_ts1.png","BBS vs INDELI all segments penalty 2")

