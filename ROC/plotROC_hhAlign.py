#!/usr/bin/env python3
import re, os, csv, sys, pprint
sys.path.append(os.path.dirname(os.path.abspath(__name__)))
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc

def read_csv(csv_location):
    '''Put csv in list'''
    output_list = list()
    first_line = True
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if first_line:
                first_line = False
                continue
            output_list.append((row[0], min(float(row[4]),float(row[12]))))
    return output_list

def calculate_tpr_fpr(csv_data, eval_thr):
    from bin.SVM_test import bypass_zero_division
    tp, tn, fp, fn = 0, 0, 0, 0
    for aln in csv_data:
        if re.match('^A_|^B_', aln[0]) is not None:
            if aln[1] > eval_thr:
                fn += 1
            if aln[1] <= eval_thr:
                tp += 1
        if re.match('^C_|^D_', aln[0]) is not None:
            if aln[1] > eval_thr:
                tn += 1
            if aln[1] <= eval_thr:
                fp += 1
    tpr = bypass_zero_division(tp,tp+fn)
    tnr = bypass_zero_division(tn,tn+fp)
    return (tpr, 1-tnr)


def construct_param_struc(split_label, datastruc, tpr, fpr):
    if split_label[1] not in datastruc.keys():
        datastruc[split_label[1]] = list()
    datastruc[split_label[1]].append((tpr, fpr))
    return datastruc

def plot_subplot(datastruc, axs, titlesubplot):
    axs.set_title(titlesubplot)
    color_idx = np.linspace(0, 1, len(datastruc.keys()))
    for key, i in zip(datastruc.keys(), color_idx):
        fpr = [x[1] for x in datastruc[key]]
        tpr = [x[0] for x in datastruc[key]]
        accuracy = auc(fpr, tpr)
        legend_label = key[1:]+"%    {:.1%}".format(float(accuracy))
        axs.plot(fpr, tpr, 
                  label=legend_label, 
                  color=plt.cm.viridis(i))
        axs.legend(loc=4, title="% cut gaps   AUC")
        axs.set_xlabel("False positive rate")
        axs.set_ylabel("True positive rate")
    return True

directory = "/home/ppenev/Dropbox-Gatech/_2019_Petar_Methods-Score/DATA/hhalign/"
bbs,rprot,indeli = dict(),dict(),dict()
for file in os.listdir(directory):
    if not re.findall(r'(.*)(\.csv)',file):
        #print("Skipping non-csv file "+file)
        continue
    csv_data = read_csv(directory+file)
    split_label = file.replace('.csv','').split(".")
    for eval_thr in np.linspace(0, 1, 1000):
        tpr, fpr = calculate_tpr_fpr(csv_data, float(eval_thr))
        if split_label[0] == "BALI":
            bbs = construct_param_struc(split_label, bbs, tpr, fpr)
        if split_label[0] == "rProt":
            rprot = construct_param_struc(split_label, rprot, tpr, fpr)
        if split_label[0] == "IND":
            indeli = construct_param_struc(split_label, indeli, tpr, fpr)

fig, axes = plt.subplots(1, 3,
                            sharey=True, 
                            figsize=(20, 12))
fig.suptitle("ROC HHalign", fontsize=16)

for datastruc, axs, title in zip([bbs,rprot,indeli], axes, ["BALI","rProtein","INDELIBLE"]):
    plot_subplot(datastruc, axs, title)

plt.savefig("./data/outputs/PNG/param_test/ROC_HHalign.png", dpi=600)