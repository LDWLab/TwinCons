#!/usr/bin/env python3
import os, sys
import numpy as np
from subprocess import STDOUT, check_output, TimeoutExpired, CalledProcessError

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from twincons import MatrixInfo
from twincons.TwinCons import subs_matrix, nucl_matrix
from twincons.MatrixLoad import PAMLmatrix

subMatrices = MatrixInfo.available_matrices

def write_temp_matrix(subsMatrix):
    with open("./matrixAdjustments/tempMx", "w") as f:
        for row in subsMatrix:
            numRow = [str(x) for x in row]
            f.write(' '.join(numRow))
            f.write('\n')

def runJointProbEstimator(mxName, secs, mxSize):
    '''Runs SJP_bioinfo to calculate the validity of a substitution matrix and discover its joint probaibilities.
    Returns True or False depending on success. Saves the output in a file named after the matrix in the folder matrixCA.'''
    cmd = f'/mnt/e/Programs/TwinCons/matrixAdjustments/SJP_bioinfo {mxSize}\
            /mnt/e/Programs/TwinCons/matrixAdjustments/tempMx\
            /mnt/e/Programs/TwinCons/matrixAdjustments/matrixCA/{mxName}'
    success = False
    try:
        output = check_output(cmd, stderr=STDOUT, shell=True, timeout=secs)
        print(f"Joint probabilitites of matrix {mxName} were discovered successfully.")
        success = True
    except TimeoutExpired:
        print(f"Timeout of {secs} seconds has expired. Skipping matrix {mxName}!")
    except CalledProcessError:
        print(f"Matrix {mxName} caused error, but its joint probabilitites were calculated successfully.")
        success = True
    except:
        print(f"Matrix {mxName} caused an unknown error. Skipping matrix {mxName}!")
    return success

mxfiles = [f for f in os.listdir(str(os.getcwd())+'/matrices/structureDerived/') if os.path.isfile(os.path.join(str(os.getcwd())+'/matrices/structureDerived/', f))]
for mxFile in mxfiles:
    mxName = mxFile.replace('.dat','')
    mx = np.array(PAMLmatrix(str(os.getcwd())+'/matrices/structureDerived/'+mxFile).lodd)
    write_temp_matrix(mx)
    succeded = runJointProbEstimator(mxName, 10, 20)
    print(succeded)

for nucMxName in ["identity", "blastn", "trans"]:
    mx = nucl_matrix(nucMxName)
    write_temp_matrix(mx)
    succeded = runJointProbEstimator(nucMxName, 10, 4)
    print(succeded)

for mxName in subMatrices:
    mx = subs_matrix(mxName)
    write_temp_matrix(mx)
    succeded = runJointProbEstimator(mxName, 10, 20)
    print(succeded)