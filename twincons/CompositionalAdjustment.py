'''
Compositional adjustment of substitution matrices (Yu, Wootton & Altschul 2003,
doi.org/10.1073/pnas.2533904100; Yu & Altschul 2005, doi.org/10.1093/bioinformatics/bti070).

Pure numpy replacement for the newton_direct_solve program previously bundled with
TwinCons, which wrapped the public-domain NCBI composition_adjustment library.
It reproduces that program's results (validated to its 9-decimal output precision) in
the mode TwinCons used: the adjusted matrix keeps the relative entropy of the original
matrix evaluated in the new compositional context.
'''
import numpy as np

DEFAULT_PSEUDOCOUNTS = 20


class CompositionalAdjustmentError(ArithmeticError):
    '''Raised when the target frequency optimization does not converge.'''


def _constraint_matrix(alphsize):
    '''A such that A @ x.ravel() gives the column sums of x followed by the sums of rows 1..n-1.'''
    indices = np.arange(alphsize * alphsize).reshape(alphsize, alphsize)
    A = np.zeros((2 * alphsize - 1, alphsize * alphsize))
    for j in range(alphsize):
        A[j, indices[:, j]] = 1.0
    for i in range(1, alphsize):
        A[alphsize + i - 1, indices[i, :]] = 1.0
    return A


def lambda_for_composition(scores, row_probs, col_probs):
    '''Solves sum_ij r_i c_j exp(lambda s_ij) = 1 as the original program did: double lambda
    from 1 until the slope is positive, then take half Newton steps until a step is <= 1e-7.'''
    weights = np.outer(row_probs, col_probs)
    lam, slope_sum = 1.0, 0.0
    while 1e-7 >= slope_sum:
        lam *= 2
        slope_sum += (scores * np.exp(lam * scores) * weights).sum()
    delta, iterations = 1.0, 0
    while abs(delta) > 1e-7 and iterations <= 299:
        terms = np.exp(lam * scores) * weights
        delta = (1.0 - terms.sum()) / (scores * terms).sum()
        lam += 0.5 * delta
        iterations += 1
    return lam


def optimize_target_frequencies(joint_probs, row_sums, col_sums, relative_entropy=None, tol=1e-8, maxits=2000):
    '''
    Port of NCBI Blast_OptimizeTargetFrequencies: finds target frequencies x closest (in
    relative entropy) to joint_probs with the given row and column sums and, unless
    relative_entropy is None, with that relative entropy. Returns (x, converged, iterations).
    '''
    alphsize = len(row_sums)
    mA = 2 * alphsize - 1
    constrain = relative_entropy is not None
    m = mA + 1 if constrain else mA
    A = _constraint_matrix(alphsize)
    b = np.concatenate([col_sums, row_sums[1:]])
    q = np.asarray(joint_probs, dtype=float).ravel()
    x = q.copy()
    z = np.zeros(mA + 1)
    iterations = 0
    rnorm = np.inf
    with np.errstate(all='ignore'):
        old_scores = np.log(q / np.outer(row_sums, col_sums).ravel())
        while iterations <= maxits:
            log_ratio = np.log(x / q)
            grad_objective = log_ratio + 1
            resids_x = -grad_objective + A.T @ z[:mA]
            resids_z = np.empty(m)
            resids_z[:mA] = b - A @ x
            if constrain:
                eta = z[mA]
                re_terms = log_ratio + old_scores
                grad_re = re_terms + 1
                resids_x += eta * grad_re
                resids_z[mA] = relative_entropy - np.dot(x, re_terms)
            rnorm = np.sqrt(np.dot(resids_x, resids_x) + np.dot(resids_z, resids_z))
            if not rnorm > tol:
                break
            iterations += 1
            if iterations > maxits:
                break
            # Newton step on the KKT system, block-reduced to J D^-1 J^T as in the NCBI code.
            Dinv = x / (1 - eta) if constrain else x.copy()
            W = np.empty((m, m))
            W[:mA, :mA] = (A * Dinv) @ A.T
            work = resids_x * Dinv
            rhs = resids_z.copy()
            rhs[:mA] -= A @ work
            if constrain:
                scaled_grad = Dinv * grad_re
                W[mA, mA] = np.dot(grad_re, scaled_grad)
                W[mA, :mA] = W[:mA, mA] = A @ scaled_grad
                rhs[mA] -= np.dot(grad_re, work)
            try:
                L = np.linalg.cholesky(W)
            except np.linalg.LinAlgError:
                break
            step_z = np.linalg.solve(L.T, np.linalg.solve(L, rhs))
            step_x = resids_x + (grad_re * step_z[mA] if constrain else 0)
            step_x = (step_x + A.T @ step_z[:mA]) * Dinv
            # Largest step (at most 1/0.95) that keeps x positive, scaled by 0.95.
            to_boundary = -x / step_x
            to_boundary = to_boundary[(to_boundary >= 0) & (to_boundary < 1 / .95)]
            alpha = 0.95 * (to_boundary.min() if to_boundary.size else 1 / .95)
            x = x + alpha * step_x
            z[:m] += alpha * step_z
    converged = iterations <= maxits and rnorm <= tol and (not constrain or z[m - 1] < 1)
    return x.reshape(alphsize, alphsize), converged, iterations


def adjust_matrix(joint_probs, freqs1, freqs2, length1, length2, pseudocounts=DEFAULT_PSEUDOCOUNTS):
    '''
    Returns the compositionally adjusted score matrix (natural-log units) for two residue
    compositions, given the joint probabilities that define a substitution matrix.
    Each composition is mixed with the matrix background using pseudocounts/(length+pseudocounts).
    '''
    jp = np.asarray(joint_probs, dtype=float)
    jp = jp / jp.sum()
    implicit_rows, implicit_cols = jp.sum(axis=1), jp.sum(axis=0)
    old_scores = np.log(jp / np.outer(implicit_rows, implicit_cols))
    f1 = np.asarray(freqs1, dtype=float)
    f2 = np.asarray(freqs2, dtype=float)
    f1, f2 = f1 / f1.sum(), f2 / f2.sum()
    w1 = pseudocounts / (length1 + pseudocounts)
    w2 = pseudocounts / (length2 + pseudocounts)
    rows = implicit_rows * w1 + (1.0 - w1) * f1
    cols = implicit_cols * w2 + (1.0 - w2) * f2
    relative_entropy = None
    # Keep the relative entropy of the original matrix in the new context. A matrix with a
    # non-negative expected score has no lambda, so it is adjusted without that constraint.
    if (old_scores * np.outer(f1, f2)).sum() < -1e-10:
        lam = lambda_for_composition(old_scores, f1, f2)
        relative_entropy = float((np.outer(f1, f2) * np.exp(lam * old_scores) * lam * old_scores).sum())
    target_freqs, converged, iterations = optimize_target_frequencies(jp, rows, cols, relative_entropy)
    if not converged:
        raise CompositionalAdjustmentError(
            f"Compositional adjustment did not converge after {iterations} iterations for these group compositions.")
    return np.log(target_freqs / np.outer(rows, cols))
