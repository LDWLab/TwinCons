#!/usr/bin/env python3
import os
import numpy as np

from twincons import MatrixInfo
from twincons.TwinCons import subs_matrix, nucl_matrix
from twincons.MatrixLoad import PAMLmatrix

def parseQIJmx(mxpath):
    f = open ( mxpath , 'r')
    triangular_mx = [[num for num in line.rstrip('\n').split(' ') ] for line in f if line.strip() != "" ]

subMatrices = MatrixInfo.available_matrices

print("Name MaxBefore MinBefore MaxAfter MinAfter")
for mxName in subMatrices:
    mx = subs_matrix(mxName)

mxfiles = [f for f in os.listdir(str(os.getcwd())+'/matrices/') if os.path.isfile(os.path.join(str(os.getcwd())+'/matrices/', f))]
for mxFile in mxfiles:
    mxName = mxFile.replace('.dat','')
    mx = np.array(PAMLmatrix(str(os.getcwd())+'/matrices/'+mxFile).lodd)

for nucMxName in ["identity", "blastn", "trans"]:
    mx = nucl_matrix(nucMxName)