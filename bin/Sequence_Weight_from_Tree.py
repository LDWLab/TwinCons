#!/usr/bin/env python3
"""Calculate importance of sequences based on phylogenetic tree."""

import os, re, sys, Bio.Align, argparse
import numpy as np
from Bio import Phylo
from Bio import AlignIO
from io import StringIO
from AlignmentGroup import AlignmentGroup
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import PhyMeas

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', help='Path to alignment file')
    parser.add_argument('-ind','--indelible_tree_file', help='Path to INDELible trees log file. If defined alignment_file, should lead to folder with alignments.', type=str)
    parser.add_argument('-split','--split_by_tree', help='Split the provided alignment in two files by the deepest branching point of a constructed tree. Path for output.', type=str)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def tree_contruct(aln_obj, nj=False):
    '''
    Constructs and returns a tree from an alignment object.
    '''
    calculator = DistanceCalculator('blosum62')
    dist_mx = calculator.get_distance(aln_obj)
    constructor = DistanceTreeConstructor()
    if nj:
        tree = constructor.nj(dist_mx)
    else:
        tree = constructor.upgma(dist_mx)
    tree.ladderize()
    return tree

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
        anc_names[anc] = os.path.commonprefix(deepestanc_to_child[anc])[:-1]
        #In the case of no common name found between sequences from 1 group
        if os.path.commonprefix(deepestanc_to_child[anc])[:-1] is '':
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
        wei_vr.append(treedepths[dict_clades[leaf]])
        #print(leaf, treedepths[dict_clades[leaf]])
        #print (leaf, 1/2**treedepths_int[dict_clades[leaf]])
    return wei_vr

def write_group_tagged_alns(file_handle, alngroup_name, alignment):
    with open(file_handle, "w") as tagged_aln_file:
        for entry in alignment:
            tagged_aln_file.write(str(">"+str(alngroup_name)+"_"+str(entry.id)+"\n"+entry.seq+"\n"))
    return True

def output_split_alignments(sliced_aln_dict, comm_args):
    alignment_file_name = comm_args.alignment_file.rsplit('/',1)[1]
    for alngroup_name, alignment in sliced_aln_dict.items():
        output_handle = comm_args.split_by_tree+str(alngroup_name)+"_"+alignment_file_name
        write_group_tagged_alns(output_handle, alngroup_name, alignment)
    return True

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    if comm_args.indelible_tree_file:
        indeli_trees = read_indeli_trees_file(comm_args.indelible_tree_file)
        for name,tree in indeli_trees.items():
            alignIO_out = PhyMeas.read_align(comm_args.alignment_file+'INDELI_TRUE_'+name+'.fas')
            deepest_anc = find_deepest_ancestors(tree)
            gapped_sliced_alns = slice_by_anc (alignIO_out, deepest_anc)
            print(name, gapped_sliced_alns)
            with open (comm_args.alignment_file+'tg_'+name+'.fas', "w") as tagged_aln_file:
                for aln_group_name, alignment in gapped_sliced_alns.items():
                    for entry in alignment:
                        tagged_aln_file.write(str(">"+str(aln_group_name)+"_"+str(entry.id)+"\n"+entry.seq+"\n"))
        sys.exit()
    alignIO_out = PhyMeas.read_align(comm_args.alignment_file)
    #sliced_alns = PhyMeas.slice_by_name(alignIO_out)

    tree = tree_contruct(alignIO_out)
    deepestanc_to_child = find_deepest_ancestors(tree)
    sliced_dict = slice_by_anc(alignIO_out, deepestanc_to_child)
    
    if comm_args.split_by_tree:
        if output_split_alignments(sliced_dict, comm_args):
            sys.exit()

    wei_vr_dict = {}
    for alngroup_name in sliced_dict:
        tree = tree_contruct(sliced_dict[alngroup_name])
        wei_vr_dict[alngroup_name] = generate_weight_vectors(tree)
    
    # for group in wei_vr_dict:
    #     print(group, sum(wei_vr_dict[group]))
    #     print(group, [x / max(wei_vr_dict[group]) for x in wei_vr_dict[group]])

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))