#!/usr/bin/env python3
import re
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

def read_align(aln_path):
    '''Reads the fasta file and gets the sequences.
    '''
    alignment = AlignIO.read(open(aln_path), "fasta")
    for record in alignment:
        record.seq = record.seq.upper()
    return alignment

def slice_by_name(unsliced_aln_obj):
    '''
    Slices an alignment into different alignments
    by first string (should hold protein name)
    Returns a dictionary with keys being the alignment
    names pointing to the alignIO objects
    '''
    prot_list=[]
    sliced_dict={}
    for entry in unsliced_aln_obj:
        prot_list.append(entry.id.split("_")[0])
    uniq_prot_list=set(prot_list)
    for prot in uniq_prot_list:
        what = MultipleSeqAlignment([])
        for entry in unsliced_aln_obj:
            if re.match(prot,entry.id.split("_")[0]):
                what.append(entry)
        sliced_dict[prot]=what
    return sliced_dict            #Iterate over the dict and create instances of AlignmentGroup