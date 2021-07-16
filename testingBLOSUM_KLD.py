#!/usr/bin/env python3
import os
import numpy as np

from twincons import MatrixInfo
from twincons.TwinCons import subs_matrix, nucl_matrix, baseline_matrix
from twincons.MatrixLoad import PAMLmatrix

def parseQIJorSIJmx(mxpath):
    '''Given a path to .qij or .sij file with triangular BLOSUM mx returns a symetric repreentation of the matrix.'''
    mxType = mxpath.split('.')[2]
    with open(mxpath, "r") as f:
        if mxType == 'qij':
            rows = f.readlines()[4:]
        elif mxType == 'sij':
            rows = f.readlines()[5:]
        else:
            raise IOError(f"Unsupported matrix {mxpath}")
    triangular_mx = [[float(num) for num in line.rstrip('\n').split(' ') if num != '' ] for line in rows if line.strip() != "" ]
    
    #Construct symetric matrix
    sym_mx=np.zeros((20,20))
    for i in range(len(triangular_mx)):
        for j in range(len(triangular_mx[i])):
            if triangular_mx[i][j] != '':
                sym_mx[j][i]=triangular_mx[i][j]
    i_lower = np.tril_indices(len(sym_mx), -1)
    sym_mx[i_lower] = sym_mx.T[i_lower]

    #Check if qij sums to near 1
    if mxType == 'qij':
        if np.isclose([sum(sum(sym_mx))],[1], rtol=1e-02) == False:
            raise ValueError(f"The target frequencies of matrix {mxpath} do not sum to 1!")

    return sym_mx

def parseBLOSUMoutFile(outFilePath):
    '''Given a path to BLOSUM .out file returns the marginal probabilities as numpy array of floats.'''
    with open(outFilePath, "r") as f:
        pis = f.readlines()[37]
    return np.array([float(pi) for pi in pis.split()])

def calculateMXentropy(qij, pipj):
    '''Given two matrices of equal size (qij and pipj) returns the sij MX entropy.'''
    qijOverpij = qij/pipj
    qijTimesSij = np.log2(qijOverpij)*qij
    return sum(sum(qijTimesSij))

def calculateMXrelativeEntropy(qij, pipj):
    '''Given two matrices of equal size (qij and pipj) returns the sij MX relative entropy.'''
    pijOverQij = pipj/qij
    pijTimesReverseSij = np.log2(pijOverQij)*pipj
    return sum(sum(pijTimesReverseSij))

def calculateSymmetricDivergence(qij, pipj):
    '''Given two matrices of equal size (qij and pipj) returns the sij MX symmetric divergence.'''
    qijMinuspipj = qij-pipj
    symmetricDivergence = qijMinuspipj*np.log2(qij/pipj)
    return sum(sum(symmetricDivergence))

def main ():
    print("Matrix name,Uniform baseline,pi baseline,Entropy,Relative entropy,Symmetric divergence,Half bit adjustment,Is KL equal to symmetric diverence")
    for mxPath in os.listdir('./matrices/BLOSUM/'):
        if mxPath.endswith('.qij'):
            mxName = mxPath.split('.')[0]
            qij = parseQIJorSIJmx(f'./matrices/BLOSUM/{mxName}.qij')
            sij = parseQIJorSIJmx(f'./matrices/BLOSUM/{mxName}.sij')
            pi = parseBLOSUMoutFile(f'./matrices/BLOSUM/{mxName}.out')
            
            #Get original BL matrix and calculate its baseline from the uniform and from the pi
            halfBit = subs_matrix(mxName)
            halfBitBase = baseline_matrix(halfBit)
            halfBitBasePi = baseline_matrix(halfBit, pi)
            mxBaseLine = halfBitBase[0][0]-halfBit[0][0]
            mxBaseLinePi = halfBitBasePi[0][0]-halfBit[0][0]

            #Check the sij rounding after scaling by 2
            isHalfBitCorrect = np.isclose(halfBit,(sij*2).round())

            #Calculate entropy and divergencies of original matrices
            mxentropy = calculateMXentropy(qij, np.outer(pi,pi.T))
            mxrelativeEntropy = calculateMXrelativeEntropy(qij, np.outer(pi,pi.T))
            mxsymmetricDivergence = calculateSymmetricDivergence(qij, np.outer(pi,pi.T))

            #Check if the sum of the two KL divergencies are equal to the symmetric divergence
            areKLequalToSymmDivergence = np.isclose([mxentropy + mxrelativeEntropy],[mxsymmetricDivergence])[0]
            outData = [str(round(x, 3)) for x in [mxBaseLine, mxBaseLinePi, mxentropy, mxrelativeEntropy, mxsymmetricDivergence]]
            print(f'{mxName},'+','.join(outData)+f',{isHalfBitCorrect.all()},{areKLequalToSymmDivergence}')
            '''
            If q is not too far from the product of its marginals, then
            the two K-L divergences will be similar in magnitude and
            approximately half of the symmetric divergence
            '''

main()

# subMatrices = MatrixInfo.available_matrices

# print("Name MaxBefore MinBefore MaxAfter MinAfter")
# for mxName in subMatrices:
#     mx = subs_matrix(mxName)

# mxfiles = [f for f in os.listdir(str(os.getcwd())+'/matrices/') if os.path.isfile(os.path.join(str(os.getcwd())+'/matrices/', f))]
# for mxFile in mxfiles:
#     mxName = mxFile.replace('.dat','')
#     mx = np.array(PAMLmatrix(str(os.getcwd())+'/matrices/'+mxFile).lodd)

# for nucMxName in ["identity", "blastn", "trans"]:
#     mx = nucl_matrix(nucMxName)