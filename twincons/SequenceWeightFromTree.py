"""Sequence weights and group definitions from phylogenetic trees of an alignment."""
import os, Bio.Align
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from twincons.twcSupportFunctions import alignment_array

def pairwise_distances(calculator, seqs1, seqs2):
    '''
    Vectorized DistanceCalculator._pairwise: returns an array whose [p, q] entry is the
    distance between seqs1[p] and seqs2[q], given as alignment_array() byte arrays.
    Uses the calculator's own scoring matrix and skip letters.
    '''
    skip_letters = [letter.encode('ascii') for letter in calculator.skip_letters]
    skip1 = np.isin(seqs1, skip_letters)
    skip2 = np.isin(seqs2, skip_letters)
    distances = np.empty((len(seqs1), len(seqs2)))
    if calculator.scoring_matrix is None:
        length = seqs1.shape[1]
        if length == 0:
            distances[:] = 1
            return distances
        # Compare against blocks of seqs2 at once, keeping each block's boolean array under ~2e7 elements.
        block = max(1, 20_000_000 // max(1, seqs1.size))
        for start in range(0, len(seqs2), block):
            stop = start + block
            valid = ~(skip1[:, None, :] | skip2[None, start:stop, :])
            score = ((seqs1[:, None, :] == seqs2[None, start:stop, :]) & valid).sum(axis=2)
            distances[:, start:stop] = 1 - score / length
        return distances

    scoring = np.asarray(calculator.scoring_matrix, dtype=float)
    lookup = np.full(256, -1)
    for index, letter in enumerate(calculator.scoring_matrix.alphabet):
        lookup[ord(letter)] = index
    codes1 = lookup[seqs1.view(np.uint8)]
    codes2 = lookup[seqs2.view(np.uint8)]
    diagonal = np.diag(scoring)
    for q in range(len(seqs2)):
        valid = ~(skip1 | skip2[q])
        bad = valid & ((codes1 < 0) | (codes2[q] < 0))
        if bad.any():
            p, position = np.argwhere(bad)[0]
            letter = seqs1[p, position] if codes1[p, position] < 0 else seqs2[q, position]
            raise ValueError(f"Bad letter '{letter.decode()}' at position '{position}'")
        c1, c2 = np.where(valid, codes1, 0), np.where(valid, codes2[q], 0)
        score = np.where(valid, scoring[c1, c2], 0).sum(axis=1)
        max_score = np.maximum(np.where(valid, diagonal[c1], 0).sum(axis=1),
                               np.where(valid, diagonal[c2], 0).sum(axis=1))
        with np.errstate(divide='ignore', invalid='ignore'):
            distances[:, q] = np.where(max_score == 0, 1, 1 - score / max_score)
    return distances

def distance_matrix(aln_obj, calculator):
    '''Same result as calculator.get_distance(aln_obj), computed with numpy.'''
    names = [record.id for record in aln_obj]
    seqs = alignment_array(aln_obj)
    lower_triangle = [[0]]
    for i in range(1, len(names)):
        lower_triangle.append(pairwise_distances(calculator, seqs[:i], seqs[i:i+1])[:, 0].tolist() + [0])
    return DistanceMatrix(names, lower_triangle)

def tree_construct(aln_obj, nj=False, nucl=False, ladderize=True, calc_mx='blosum62'):
    '''
    Constructs and returns a tree from an alignment object.
    '''
    if nucl:
        if calc_mx == 'blosum62':
            calc_mx = 'blastn'
        calculator = DistanceCalculator(calc_mx)
    else:
        calculator = DistanceCalculator(calc_mx)
    dist_mx = distance_matrix(aln_obj, calculator)
    constructor = DistanceTreeConstructor()
    if nj:
        tree = constructor.nj(dist_mx)
    else:
        tree = constructor.upgma(dist_mx)
    if ladderize:
        tree.ladderize()
    return tree

def voronoi_convergence(distances):
    '''Given sample x sequence distances, splits one vote per sample among its closest sequences.'''
    closest = distances == distances.min(axis=1, keepdims=True)
    return (closest / closest.sum(axis=1, keepdims=True)).sum(axis=0)

def leaf_distance_sums(tree, names):
    '''For each named leaf, the sum of tree distances to all leaves in names.
    Equivalent to summing tree.distance() over all pairs, in O(n * depth).'''
    depths = tree.depths()
    parents = {child: clade for clade in tree.find_clades() for child in clade.clades}
    leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
    targets = {leaves[name] for name in names}
    below = {}
    for clade in tree.find_clades(order='postorder'):
        below[clade] = (clade in targets) + sum(below[child] for child in clade.clades)
    total_depth = sum(depths[leaves[name]] for name in names)
    sums = list()
    for name in names:
        leaf = leaves[name]
        lca_depth_sum = depths[leaf] * below[leaf]
        child, clade = leaf, parents.get(leaf)
        while clade is not None:
            lca_depth_sum += depths[clade] * (below[clade] - below[child])
            child, clade = clade, parents.get(clade)
        sums.append(len(names) * depths[leaf] + total_depth - 2 * lca_depth_sum)
    return sums

def branch_sharing_weights(tree, names):
    '''
    Thompson, Higgins & Gibson (1994) weights, as used by ClustalW: each branch length is shared
    equally by the leaves below it, and a leaf's weight is the sum of its shares along the path
    from the root. Closely related sequences split their common branches, so redundant
    sequences do not accumulate weight.
    '''
    leaves_below = {}
    for clade in tree.find_clades(order='postorder'):
        leaves_below[clade] = 1 if clade.is_terminal() else sum(leaves_below[child] for child in clade.clades)
    share = {tree.root: 0.0}
    for clade in tree.find_clades(order='preorder'):
        for child in clade.clades:
            share[child] = share[clade] + (child.branch_length or 0.0) / leaves_below[child]
    leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
    return [share[leaves[name]] for name in names]

DEFAULT_VORONOI_SAMPLES = 100000

def _one_hot(codes, symbols):
    '''(rows, columns) symbol codes -> float32 (rows, columns * symbols) one-hot matrix.'''
    rows, columns = codes.shape
    matrix = np.zeros((rows, columns * symbols), dtype=np.float32)
    matrix[np.repeat(np.arange(rows), columns), (np.arange(columns) * symbols + codes).ravel()] = 1
    return matrix

def voronoi_weights(aln_obj, repeat=DEFAULT_VORONOI_SAMPLES, calc_mx='identity'):
    '''
    Voronoi sequence weights (Sibbald & Argos 1990): random sequences pick, at every column,
    one of the characters present there with equal probability; each random sequence gives one
    vote, split between the alignment sequences closest to it. Uses numpy's global random state.
    '''
    calculator = DistanceCalculator(calc_mx)
    seqs = alignment_array(aln_obj)
    n, length = seqs.shape
    symbols, codes = np.unique(seqs, return_inverse=True)
    codes = codes.reshape(n, length)
    column_codes = [np.unique(codes[:, column]) for column in range(length)]
    choices_per_column = np.array([len(column) for column in column_codes])
    padded_codes = np.zeros((length, choices_per_column.max()), dtype=codes.dtype)
    for column, column_symbols in enumerate(column_codes):
        padded_codes[column, :len(column_symbols)] = column_symbols
    # Identity distances depend only on the number of identical columns, which a one-hot matrix
    # product counts exactly; other distance models go through pairwise_distances.
    fast_identity = calculator.scoring_matrix is None and not calculator.skip_letters
    if fast_identity:
        seq_one_hot = _one_hot(codes, len(symbols))
    chunk_size = max(1, 25_000_000 // max(1, length * len(symbols)))
    convergence_vr = np.zeros(n)
    for start in range(0, repeat, chunk_size):
        samples = min(chunk_size, repeat - start)
        picks = (np.random.random_sample((samples, length)) * choices_per_column).astype(int)
        sample_codes = padded_codes[np.arange(length), picks]
        if fast_identity:
            matches = _one_hot(sample_codes, len(symbols)) @ seq_one_hot.T
            distances = 1 - matches.astype(float) / length
        else:
            distances = pairwise_distances(calculator, seqs, symbols[sample_codes]).T
        convergence_vr += voronoi_convergence(distances)
    return (convergence_vr / convergence_vr.sum()).tolist()

WEIGHTING_ALGORITHMS = ['pairwise', 'voronoi', 'clustalw']

def calculate_weight_vector(aln_obj, algorithm='pairwise', calc_mx='identity', repeat=DEFAULT_VORONOI_SAMPLES, nucl=False):
    '''Returns one weight per sequence of aln_obj, summing to 1.'''
    if algorithm not in WEIGHTING_ALGORITHMS:
        raise ValueError("Invalid algorithm type. Expected one of: %s" % WEIGHTING_ALGORITHMS)
    if algorithm == 'voronoi':
        return voronoi_weights(aln_obj, repeat=repeat, calc_mx=calc_mx)
    if len(aln_obj) == 1:
        return [1.0]
    tree = tree_construct(aln_obj, nucl=nucl, calc_mx=calc_mx)
    names = [seq_obj.id for seq_obj in aln_obj]
    if algorithm == 'pairwise':
        raw_weights = leaf_distance_sums(tree, names)
    else:
        raw_weights = branch_sharing_weights(tree, names)
    total = sum(raw_weights)
    # Identical sequences give a tree without branch lengths; they are then equally weighted.
    if total == 0:
        return [1 / len(names)] * len(names)
    return [weight / total for weight in raw_weights]

def find_deepest_ancestors(tree):
    '''
    Finds and returns a dictionary with keys the deepest ancestors and values their children.
    '''
    treedepths_int = tree.depths(unit_branch_lengths=True)
    deepestanc_to_child={}
    for anc in treedepths_int:
        if treedepths_int[anc] == 1:
            for child in tree.get_terminals():
                if anc.is_parent_of(child):
                    if anc not in deepestanc_to_child:
                        deepestanc_to_child[anc]=[]
                    deepestanc_to_child[anc].append(child.name)

    return deepestanc_to_child

def slice_by_anc(unsliced_aln_obj, deepestanc_to_child):
    '''
    Slices an alignment into different alignments
    by the groupings defined in deepestanc_to_child.
    '''
    anc_names={}
    sliced_dict={}
    i = 1
    for anc in deepestanc_to_child:
        #print(anc, os.path.commonprefix(deepestanc_to_child[anc]).replace("_",""))
        anc_names[anc] = os.path.commonprefix(deepestanc_to_child[anc]).replace("_","")
        #In the case of no common name found between sequences from 1 group
        if os.path.commonprefix(deepestanc_to_child[anc])[:-1] == '':
            anc_names[anc] = i
            i += 1
    
    if len(anc_names) != 2:
        raise ValueError("For now does not support more than two groups! Offending groups are "+str(', '.join(anc_names.values())))
    
    for anc in deepestanc_to_child:
        what = Bio.Align.MultipleSeqAlignment([])
        for entry in unsliced_aln_obj:
            if entry.id in deepestanc_to_child[anc]:
                what.append(entry)
        sliced_dict[anc_names[anc]]=what
    return sliced_dict
