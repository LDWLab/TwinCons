#!/usr/bin/env python3
"""Calculate importance of sequences based on phylogenetic tree."""

import os, sys, Bio.Align, argparse
import numpy as np
from Bio import Phylo
from io import StringIO
from numpy.random import choice
from collections import Counter
from statistics import mean, stdev
from itertools import combinations, product
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from twincons.twcSupportFunctions import read_align, slice_by_name, alignment_array

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', help='Path to alignment file')
    parser.add_argument('-ind','--indelible_tree_file', help='Path to INDELible trees log file. If defined alignment_file, should lead to folder with alignments.', type=str)
    parser.add_argument('-split','--split_by_tree', help='Split the provided alignment in two files by the deepest branching point of a constructed tree. Path for output.', action="store_true")
    parser.add_argument('-save','--save_split_alignments', help='Path for output of split trees', type=str)
    parser.add_argument('-nc','--nucleotide_alignments', help='Alignments are nucleotidic', action="store_true")
    parser.add_argument('-d','--dhat', help='Calculate D average.', action="store_true")
    parser.add_argument('-w','--weight_by_tree', help='Weight sequences by tree topology.', action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

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
        for q in range(len(seqs2)):
            valid = ~(skip1 | skip2[q])
            score = ((seqs1 == seqs2[q]) & valid).sum(axis=1)
            distances[:, q] = 1 - score / length if length else 1
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

def generate_sequence_sampled_from_alignment(aln_obj):
    outseq = ''
    i = 0
    while i < len(aln_obj[0]):
        aa_list = list(set(aln_obj[:, i]))
        distribution = Counter(aln_obj[:, i])
        choice_distr = []
        for aa in aa_list:
            choice_distr.append(distribution[aa]/len(aln_obj[:, i]))
        outseq += choice(aa_list)
        i+=1
    return outseq

def voronoi_convergence(distances):
    '''Given sample x sequence distances, splits one vote per sample among its closest sequences.'''
    convergence_vr = np.zeros(distances.shape[1])
    for sample_distances in distances:
        closest = sample_distances == sample_distances.min()
        convergence_vr[closest] += 1/closest.sum()
    return convergence_vr

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

def calculate_weight_vector(aln_obj, algorithm='pairwise', calc_mx='identity', repeat=1000, nucl=False, chunk_size=100):
    alg_types = ['voronoi', 'pairwise']
    if algorithm not in alg_types:
        raise ValueError("Invalid algorithm type. Expected one of: %s" % alg_types)
    if algorithm == 'voronoi':
        # Random sequences pick, at every column, one of the characters present there with equal probability.
        calculator = DistanceCalculator(calc_mx)
        seqs = alignment_array(aln_obj)
        column_chars = [np.unique(seqs[:, column]) for column in range(seqs.shape[1])]
        choices_per_column = np.array([len(chars) for chars in column_chars])
        padded_chars = np.full((len(column_chars), choices_per_column.max()), b'-', dtype='S1')
        for column, chars in enumerate(column_chars):
            padded_chars[column, :len(chars)] = chars
        convergence_vr = np.zeros(len(aln_obj))
        for start in range(0, repeat, chunk_size):
            samples = min(chunk_size, repeat - start)
            picks = (np.random.random_sample((samples, len(column_chars))) * choices_per_column).astype(int)
            test_seqs = padded_chars[np.arange(len(column_chars)), picks]
            convergence_vr += voronoi_convergence(pairwise_distances(calculator, seqs, test_seqs).T)
        return (convergence_vr / convergence_vr.sum()).tolist()
    if algorithm == 'pairwise':
        tree = tree_construct(aln_obj, nucl=nucl, calc_mx=calc_mx)
        distance_sums = leaf_distance_sums(tree, [seq_obj.id for seq_obj in aln_obj])
        return [i/sum(distance_sums) for i in distance_sums]

def read_indeli_trees_file(file_path):
    '''Reads indelible trees.txt and outputs a dictionary with keys the REPs
    and values Bio.Phylo tree objects based on the tree string.
    '''
    with open (file_path) as trees_file:
        lines = trees_file.read().splitlines()
    tree_data = lines[6:]
    output_dict = {}
    for tree_entry in tree_data:
        tree = Phylo.read(StringIO(tree_entry.split('\t')[8]), "newick")
        output_dict[tree_entry.split('\t')[3]] = tree
    return output_dict

def find_deepest_ancestors(tree):
    '''
    Finds and returns a dictionary with keys the deepest ancestors and values their children.
    '''
    treedepths_int = tree.depths(unit_branch_lengths=True)
    treedepths = tree.depths()
    common_ancestors=[]
    deepestanc_to_child={}
    for anc in treedepths:
        if treedepths_int[anc] == 1:
            #print(anc,treedepths[anc],treedepths_int[anc])
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

def generate_weight_vectors(tree):
    treedepths_int = tree.depths(unit_branch_lengths=True)
    treedepths = tree.depths()
    dict_clades = {}
    for terminal_leaf in tree.get_terminals():
        dict_clades[terminal_leaf.name] = terminal_leaf
    wei_vr = []
    for leaf in sorted(dict_clades):
        #wei_vr.append((1/2**treedepths_int[dict_clades[leaf]]))        #Accounting for topology of tree
        #print (leaf, 1/2**treedepths_int[dict_clades[leaf]])
        wei_vr.append(treedepths[dict_clades[leaf]])
        print(leaf, treedepths[dict_clades[leaf]])
    return wei_vr

def calc_d_hat(tree, seq_names1, seq_names2):

    named_intra1, intra1 = pairwise_dist(tree, seq_names1)
    named_intra2, intra2 = pairwise_dist(tree, seq_names2)
    named_inter_dist, inter_dist = dict(), list()
    for p1, p2 in product(seq_names1, seq_names2):
        d = tree.distance(p1, p2)
        inter_dist.append(d)
        named_inter_dist[(p1, p2)] = d

    dhat = ((mean(intra1)+mean(intra2))/2)/mean(inter_dist)
    return dhat

def pairwise_dist(tree, seq_names):
    dmat, dlist = {}, []
    for p1, p2 in combinations(seq_names, 2):
        d = tree.distance(p1, p2)
        dmat[(p1, p2)] = d
        dlist.append(d)
    return dmat, dlist

def write_group_tagged_alns(file_handle, alngroup_name, alignment):
    with open(file_handle, "w") as tagged_aln_file:
        for entry in alignment:
            tagged_aln_file.write(str(">"+str(alngroup_name)+"_"+str(entry.id)+"\n"+entry.seq+"\n"))
    return True

def output_split_alignments(sliced_aln_dict, output_path, aln_file):
    alignment_file_name = aln_file.rsplit('/',1)[1]
    for alngroup_name, alignment in sliced_aln_dict.items():
        output_handle = output_path+str(alngroup_name)+"_"+alignment_file_name
        write_group_tagged_alns(output_handle, alngroup_name, alignment)
    return True

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    aln_name = str(comm_args.alignment_file).split("/")[-1].replace(".fasta", "").replace(".fas", "")
    if comm_args.indelible_tree_file:
        indeli_trees = read_indeli_trees_file(comm_args.indelible_tree_file)
        for name,tree in indeli_trees.items():
            alignIO_out = read_align(comm_args.alignment_file+'INDELI_TRUE_'+name+'.fas')
            deepest_anc = find_deepest_ancestors(tree)
            gapped_sliced_alns = slice_by_anc (alignIO_out, deepest_anc)
            print(name, gapped_sliced_alns)
            with open (comm_args.alignment_file+'tg_'+name+'.fas', "w") as tagged_aln_file:
                for aln_group_name, alignment in gapped_sliced_alns.items():
                    for entry in alignment:
                        tagged_aln_file.write(str(">"+str(aln_group_name)+"_"+str(entry.id)+"\n"+entry.seq+"\n"))
        sys.exit()
    alignIO_out = read_align(comm_args.alignment_file)
    if comm_args.nucleotide_alignments:
        for sequence in alignIO_out:
            sequence.seq = sequence.seq.back_transcribe()
        tree = tree_construct(alignIO_out, nucl=True)
    else:
        tree = tree_construct(alignIO_out)

    if comm_args.split_by_tree:
        deepestanc_to_child = find_deepest_ancestors(tree)
        sliced_dict = slice_by_anc(alignIO_out, deepestanc_to_child)
        if comm_args.save_split_alignments:
            output_split_alignments(sliced_dict, comm_args.save_split_alignments, comm_args.alignment_file)
            sys.exit()
    else:
        sliced_dict = slice_by_name(alignIO_out)

    if comm_args.dhat:
        # #D hat here figure out why it always ranges between 0.4-0.7
        seq_names = list()
        group_names = list()
        for group in sliced_dict:
            if group[-1] == 'b':
                continue
            group_names.append(group)
            seq_names.append([x.id for x in sliced_dict[group]])
        dhat = calc_d_hat(tree, seq_names[0], seq_names[1])
        named_intra1, intra1 = pairwise_dist(tree, seq_names[0])
        named_intra2, intra2 = pairwise_dist(tree, seq_names[1])
        print(aln_name, "\t", ' '.join(group_names), dhat)
        print(named_intra1, intra1)
        print(named_intra2, intra2)
        sys.exit()

    wei_vr_dict, pairwise_dict, pairwise_lists = {}, {}, {}
    for alngroup_name in sliced_dict:
        tree = tree_construct(sliced_dict[alngroup_name])
        seq_names = tree.get_terminals()
        pairwise_dict[alngroup_name], pairwise_lists[alngroup_name] = pairwise_dist(tree, seq_names)
        wei_vr_dict[alngroup_name] = generate_weight_vectors(tree)
    
    print (aln_name, end="\t")
    for group in pairwise_lists:
        print(group, mean(pairwise_lists[group]),stdev(pairwise_lists[group]), end="\t")
        #print(group, [x for x in pairwise_lists[group]])
    print()
    # for group in wei_vr_dict:
    #     print(group, sum(wei_vr_dict[group])/len(wei_vr_dict[group]))
    #     print(group, [x for x in wei_vr_dict[group]])
    #     #print(group, [x / max(wei_vr_dict[group]) for x in wei_vr_dict[group]])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))