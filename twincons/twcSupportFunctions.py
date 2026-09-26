#!/usr/bin/env python3
import shutil
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

def find_executable(name, purpose):
    '''Returns the full path of an external program, or raises a helpful error.'''
    executable = shutil.which(name)
    if executable is None:
        raise OSError(f"{name} was not found on PATH. It is required {purpose}.")
    return executable

def read_align(aln_path):
    '''Reads the fasta file and gets the sequences.
    '''
    with open(aln_path) as aln_handle:
        alignment = AlignIO.read(aln_handle, "fasta")
    for record in alignment:
        record.seq = record.seq.upper()
    return alignment

def slice_by_name(unsliced_aln_obj):
    '''
    Slices an alignment into different alignments by the sequence id prefix
    before the first underscore. Returns a dictionary, ordered by first
    appearance in the alignment, of group name -> MultipleSeqAlignment.
    '''
    sliced_dict = {}
    for entry in unsliced_aln_obj:
        group_name = entry.id.split("_")[0]
        if group_name not in sliced_dict:
            sliced_dict[group_name] = MultipleSeqAlignment([])
        sliced_dict[group_name].append(entry)
    return sliced_dict