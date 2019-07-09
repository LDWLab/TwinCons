#!/usr/bin/env python3
'''
SVM from alignment segments
Two modes:
    Train: Computes a decision function from csv generated with MultiPhyMeas
    Test: Evaluates alignment entries in csv with a provided decision function
'''
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
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('csv_path', help='Path to csv file storing alignment segment data')
    mode.add_argument('-tr','--train',help='Train on provided alignment segment data.', action="store_true")
    mode.add_argument('-te','--test',help='Test provided alignment segment data.', action="store_true")
    parser.add_argument('-dp','--decision_function_path', help='Path to pre-computed decision function. Required if running in test mode.', required=False)
    parser.add_argument('-mf','--max_features',nargs='+',help='If running in test mode, the program requires the maximal values for each feature class used to generate the decision funciton.')
    commandline_args = parser.parse_args(argument_list)
    return commandline_args, parser

def load_our_data(csv_location):
    '''
    Reads a csv outputted from MultiPhyMeas.py, entries starting with A_ and B_ are considered as + data
    and entries starting with C_ and D_ as - data. Normalizes all the data for x and y in the range 0,1
    using the max(x) and max(y).
    '''
    data_xyz = []
    data_identity = []
    data_weights = []
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                if re.match(r'^A_|^B_',row[0]):
                    data_xyz.append([float(row[2]),float(row[3]),float(row[1])])
                    data_identity.append(1)
                    data_weights.append(row[1])
                elif re.match(r'^C_|^D_',row[0]):
                    data_xyz.append([float(row[2]),float(row[3]),float(row[1])])
                    data_identity.append(0)
                    data_weights.append(float(row[1]))
            line_count+=1
    data_xyz_normx = []
    for tups in data_xyz:
        data_xyz_normx.append([float(tups[0]/max(np.asarray(data_xyz)[:,0])),float(tups[1]/max(np.asarray(data_xyz)[:,1])),float(tups[2]/max(np.asarray(data_xyz)[:,2]))])
    return np.asarray(data_xyz_normx), np.asarray(data_identity), data_weights, max(np.asarray(data_xyz)[:,0]), max(np.asarray(data_xyz)[:,1]), max(np.asarray(data_xyz)[:,2])

def plot_decision_function(classifier,X,y, sample_weight, axis, title):
    '''
    Plots the data and the decision function. Besides the classifier function, takes in the sample weights
    and plots it as a size of the datapoints. If they are different than 1.
    '''
    # plot the decision function
    xx, yy, zz = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    #xx, yy = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel(), zz.ravel()])
    print(classifier.score(X, y))
    Z = Z.reshape(xx.shape)
    Z = Z[:,:,0]
    xx = xx[:,:,0]
    yy = yy[:,:,0]
    
    # plot the line, the points, and the nearest vectors to the plane
    axis.contour(xx, yy, Z, colors='r', levels=[0], alpha=0.5, linestyles=['-'])
    #axis.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone)
    if sample_weight[0] == 1:           #Needs fixing, what if the first is one but the rest are something else?
        axis.scatter(X[:, 0], X[:, 1], c=y, alpha=0.9,
                    cmap=plt.cm.viridis, edgecolors='black')
    else:
        abs_length = [float(n)**2 for n in sample_weight]
        axis.scatter(X[:, 0], X[:, 1], c=y, s=abs_length, alpha=0.9,
                    cmap=plt.cm.bone, edgecolors='black')

    #axis.axis('off')
    axis.set_title(title)

def main(commandline_arguments):
    '''Main entry point'''
    comm_args, parser = create_and_parse_argument_options(commandline_arguments)
    if comm_args.train:
        #Load our data
        X = []
        y = []
        X, y, sample_weight, maxX, maxY, maxZ = load_our_data(comm_args.csv_path)
        sample_weight_constant = np.ones(len(X))

        print("Max on X axis:", maxX, "\nMax on Y axis:", maxY, "\nMax on Z axis:", maxZ)
        
        ### Fit the model   ###
        clf_weights = svm.SVC(C=10, gamma='auto')
        clf_weights.fit(X, y, sample_weight=sample_weight)

        ###   Save the classifier   ###
        with open('longtest_classifier.pkl', 'wb') as fid:
            cPickle.dump(clf_weights, fid)

    if comm_args.test:
        if not comm_args.decision_function_path or not comm_args.max_features:
            parser.print_help()
            raise ValueError("In test mode decision_function_path and max_features are required arguments!")
        

#fig, axes = plt.subplots(1, 1, figsize=(14, 6))
#plot_decision_function(clf_weights,X,y, sample_weight, axes[1],
#                       "Modified weights")
#plt.savefig(sys.argv[2], dpi=600)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))