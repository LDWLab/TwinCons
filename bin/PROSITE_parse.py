#!/usr/bin/env python3
"""Crete mafft combination commands from directory with alignments."""
import os, re, sys, argparse, itertools, glob, random
sys.path.append(os.path.dirname(os.path.abspath(__name__)))
import numpy as np
from Bio import AlignIO


def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Merge alignments from a given directory')
	parser.add_argument('output_directory', help='Path for outputting mergers.')
	parser.add_argument('-p','--prosite_path', help='Path to PROSITE folder.')
	parser.add_argument('-prst_choice','--prosite_choice', help='Good, bad, or both merges.', choices=['good', 'bad', 'both'])
	parser.add_argument('-i','--indeli_path', help='Path to INDELI folder.')
	parser.add_argument('-e','--ecod_path', help='Path to ECOD families folder.')
	parser.add_argument('-el','--ecod_level', help='Architecture level to merge', type=int, default=4)
	parser.add_argument('-rp','--rprot_path', help='Path to rProteins folder.')
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def read_parse_doc_file(comm_args):
	'''
	Reads the doc file from PROSITE and creates a dictionary of keys PDCs and values a list of associated PS families.
	'''
	docfile = comm_args.prosite_path+"/prosite.doc"
	pdc_to_ps={}
	with open (docfile, 'r', encoding="latin-1") as in_file:
		data = in_file.read().replace('\n', ' ')
	indexes = re.split('{END} | {BEGIN}',data)
	del indexes[1::2]
	indexes = list(filter(None, indexes))
	for fam in indexes:
		elements_in_fam = list(filter(None, re.split('{|}',fam)))
		elements_in_fam = list(filter(lambda element: element.strip(), elements_in_fam))
		elements_in_fam = list(filter(lambda element: not re.compile(r'//').search(element), elements_in_fam))
		for index, element in enumerate(elements_in_fam):
			elements_in_fam[index] = re.sub(r';.*', '', elements_in_fam[index])
		pdc_to_ps[elements_in_fam.pop(0)] = elements_in_fam
	return pdc_to_ps

def parse_indeli_folder(comm_args):
	'''
	Creates a dictionary of sequences from a folder with sequences, for mafft merging.
	Makes sure alignments with the same number of sequences are grouped together. (ConSurf becnhmarking files)
	'''
	alnfiles_list = os.listdir(comm_args.indeli_path)
	pdc_to_ps={}
	for alnfile in alnfiles_list:
		if alnfile.split('.')[0] not in pdc_to_ps:
			pdc_to_ps[alnfile.split('.')[0]]=[]
		pdc_to_ps[alnfile.split('.')[0]].append(comm_args.indeli_path+alnfile)
	return(pdc_to_ps)

def parse_indeli_folder_mine(comm_args):
	'''
	Creates a dictionary of sequences from a folder with sequences, for mafft merging.
	Makes sure alignments with the same number of sequences are grouped together. (ConSurf becnhmarking files)
	'''
	alnfiles_list = os.listdir(comm_args.indeli_path)
	pdc_to_ps={}
	pdc_to_ps[1] = alnfiles_list
	return(pdc_to_ps)

def parse_rprot_folder(comm_args):
	'''
	Combines rProtein alignments and INDELI results into data structure for mafft merging module.
	'''
	indeli_list = glob.glob(comm_args.indeli_path+"*.fas")
	rprot_list = glob.glob(comm_args.rprot_path+"*.fas")
	pdc_to_ps = {}
	for pair in itertools.product(indeli_list, rprot_list):
		pdc_to_ps[pair] = [pair[0],pair[1]]
	return pdc_to_ps

def hash_constructor(data, aln):
	if len(data) == 0:	#trivial case, we have no element therefore we return empty list
		return aln 
	else: 				#if we have elements
		first_value = data[0] #we take the first value
		data = {first_value : hash_constructor(data[1:],aln)} #data[1:] will return a list with every value but the first value
		return data #this is called after the last recursion is called

def build_nested(segs, text, container):
	'''
	Recursive constructor of a hash of hashes meant to organize the X -> H -> T -> F -> file_path.
	Container gets augmented to hold that structure.
	In cases of overlapping and similar Family levels (e.g 1.1.1.2 and 1.1.1.2-1 - two alignments of one family)
	it uses a tuple to store the file paths
	'''
	head = segs[0]
	tail = segs[1:]
	if not tail:
		if head not in container:
			container[head] = text
		else:
			container[head] = (container[head], text)
	else:
		if head not in container:
			container[head] = {}
		build_nested(tail, text, container[head])

def parse_ecod_folder(comm_args):
	'''
	Produces combinations of family alignments for specified ECOD level.
	'''
	aln_list = glob.glob(comm_args.ecod_path+"*.fas")
	fnamelist_to_alnpath = {}
	for alignment in aln_list:
		#parses through the aln filename to grab only the architecture levels in a fname list
		fname = re.split('\.|-', alignment.split("/")[int(len(alignment.split("/"))-1)])[:-1]
		build_nested(fname, alignment, fnamelist_to_alnpath)
	
	for xlevel in sorted(fnamelist_to_alnpath):
		print('X level: ',xlevel)
		for hlevel in fnamelist_to_alnpath[xlevel]:
			print('\tH level ',hlevel)
			for tlevel in fnamelist_to_alnpath[xlevel][hlevel]:
				print ('\t\tT level ',tlevel, fnamelist_to_alnpath[xlevel][hlevel][tlevel])


def merger_commands(sequence_groups,output_path):
	'''
	Executes mafft merging with a provided groups of sequences.
	Writes out the resulting mergers to a specified directory.
	'''
	for group, seq_list in sequence_groups.items():
			print(group)
			for fam_comb in itertools.combinations(seq_list, 2):
				#This one is unique
				alns_for_merging = "concat_"+re.sub(r'.msa|.fas|.*/','',fam_comb[0])+"_"+re.sub(r'.*/','',fam_comb[1])
				os.system("cat "+fam_comb[0]+" "+fam_comb[1]+" > "+alns_for_merging)
				os.system("ruby /usr/local/bin/makemergetable.rb "+fam_comb[0]+" "+fam_comb[1]+" > subMSAtable")
				os.system("mafft --merge subMSAtable "+alns_for_merging+" > "+output_path+re.sub(r'concat_','',alns_for_merging))
				os.system("rm "+alns_for_merging)
	return True

def prosite_cleanup(list_with_ids, prst_loc):
	from bin.TwinCons import read_align
	from os import listdir
	from os.path import isfile, join
	prst_aln_dir = f"{prst_loc}/prosite_alignments_handle/"
	prst_alns = [f.replace(".msa", '') for f in listdir(prst_aln_dir) if isfile(join(prst_aln_dir, f))]
	present_alns = list(set(prst_alns).intersection(set(list_with_ids)))
	filtered_ids = list()
	for prst_id in present_alns:
		aln = read_align(f"{prst_loc}/prosite_alignments_handle/{prst_id}.msa")
		if (20 < len(aln) < 100) and (50 < aln.get_alignment_length() < 500):
			filtered_ids.append(prst_id)
	return filtered_ids

def prosite_merge(list_with_bad, prst_aln_loc, out_dir):
	'''Specify whether its good or bad in the out_dir variable.'''
	family_combinations = list(set(list(itertools.combinations(list_with_bad, 2))))
	for fam_comb in family_combinations:
		alns_for_merging = f"concat_{fam_comb[0]}_{fam_comb[1]}.fas"
		file1 = f"{prst_aln_loc}/{fam_comb[0]}.msa"
		file2 = f"{prst_aln_loc}/{fam_comb[1]}.msa"
		os.system(f"cat {file1} {file2} > {out_dir}/{alns_for_merging}")
		os.system(f"ruby /usr/local/bin/makemergetable.rb {file1} {file2} > subMSAtable")
		os.system(f"mafft --quiet --merge subMSAtable {out_dir}/{alns_for_merging} > {out_dir}/{fam_comb[0]}_{fam_comb[1]}.fas")
		os.system(f"rm {out_dir}/{alns_for_merging}")
	os.system("rm subMSAtable")
	return True

def main(commandline_arguments):
	'''Main entry point'''
	comm_args = create_and_parse_argument_options(commandline_arguments)
	
	if comm_args.prosite_path:
		if not comm_args.prosite_choice:
			raise IOError("Choose what type of merge to do with -prst_choice")
		pdc_to_ps = read_parse_doc_file(comm_args)
		prosite_bad_list = list()
		prosite_good_list = list()
		for pdc, ps_list in pdc_to_ps.items():
			if len(ps_list) == 1:
				prosite_bad_list.append(ps_list[0])
			elif len(ps_list) > 1:
				prosite_good_list.append(ps_list)
		
		if comm_args.prosite_choice == 'bad' or comm_args.prosite_choice == 'both':
			bad_filtered_ids = prosite_cleanup(prosite_bad_list, comm_args.prosite_path)
			prosite_merge(bad_filtered_ids, 
						f"{comm_args.prosite_path}/prosite_alignments_handle", 
						f"{comm_args.output_directory}/BAD")
		
		if comm_args.prosite_choice == 'good' or comm_args.prosite_choice == 'both':
			good_filtered_ids = list()
			for unfiltered_good in prosite_good_list:
				temp_filtered = list()
				temp_filtered = prosite_cleanup(unfiltered_good, comm_args.prosite_path)
				if len(temp_filtered) > 1:
					good_filtered_ids.append(temp_filtered)

			for filtered_good in good_filtered_ids:
				prosite_merge(filtered_good, 
							f"{comm_args.prosite_path}/prosite_alignments_handle", 
							f"{comm_args.output_directory}/GOOD")
	
	elif comm_args.rprot_path and comm_args.indeli_path:
		pdc_to_ps = parse_rprot_folder(comm_args)
		merger_commands(pdc_to_ps, comm_args.output_directory)
	elif comm_args.indeli_path:
		pdc_to_ps = parse_indeli_folder(comm_args)
		merger_commands(pdc_to_ps, comm_args.output_directory)
	elif comm_args.ecod_path:
		parse_ecod_folder(comm_args)
	


if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))