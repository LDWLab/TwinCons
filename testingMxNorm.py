#!/usr/bin/env python3
import os
import numpy as np

from twincons import MatrixInfo
from twincons.TwinCons import subs_matrix, nucl_matrix
from twincons.MatrixLoad import PAMLmatrix


def normalizeMx(mx):
    denom = mx.max()-mx.min()
    zeroAfterNorm = (1-mx.min())/(mx.max()-mx.min())
    normMX = np.divide(np.add(mx,-1*mx.min()),denom)
    return np.add(normMX,-1*zeroAfterNorm)*10

subMatrices = MatrixInfo.available_matrices

print("Name MaxBefore MinBefore MaxAfter MinAfter")
for mxName in subMatrices:
    mx = subs_matrix(mxName)
    mxNormAroundZero = normalizeMx(mx)
    print(mxName, mx.max(), mx.min(), mxNormAroundZero.max(), mxNormAroundZero.min())

mxfiles = [f for f in os.listdir(str(os.getcwd())+'/matrices/') if os.path.isfile(os.path.join(str(os.getcwd())+'/matrices/', f))]
for mxFile in mxfiles:
    mxName = mxFile.replace('.dat','')
    mx = np.array(PAMLmatrix(str(os.getcwd())+'/matrices/'+mxFile).lodd)
    mxNormAroundZero = normalizeMx(mx)
    print(mxName, mx.max(), mx.min(), mxNormAroundZero.max(), mxNormAroundZero.min())

for nucMxName in ["identity", "blastn", "trans"]:
    mx = nucl_matrix(nucMxName)
    mxNormAroundZero = normalizeMx(mx)
    print(nucMxName, mx.max(), mx.min(), mxNormAroundZero.max(), mxNormAroundZero.min())