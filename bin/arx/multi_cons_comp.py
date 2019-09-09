#!/usr/bin/env python3
import re, os, sys, getopt, plotly, single_cons_comp
import plotly.graph_objs as go
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Name of the pdb file should start with 5 letters of the anchor sequence organism
def usage():
	print (\
	"USAGE:\n./multi_cons_comp.py -a [alignment_file_path] -p [protein_struc_path] -o [output_file_path] -h\n\
	-a: defines path to folder with alignment files. Works only on fasta type of alignments.\tREQUIRED\n\
	-o: defines output image path. Default is ./test_set/output/alignment_file_name.png\n\
	-g: columns with gaps greater than this percentage will be removed; default 1(decimal between 0 and 1)\n\
	    Setting this to 1 will substitute all gaps will random selection from the 20 amino-acid residues.\n\
	-t: threshold for number of allowed bad scores when calculating length of positive sections; default 1 (integer)\n\
	-m: what substitution matrix to use (blosum45, pam50, etc.) Default is blosum62.\n\
	-h: prints this\
")

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:o:g:t:m:h', ['alignment=', 'output=', 'percentage=', 'threshold=', 'matrix=', 'help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

output_path = 'test.png'
my_mx = 'blosum62'
gap_perc = 1
thr_distr = 1
for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-a', '--alignment'):
		input_path = arg
	elif opt in ('-o', '--output'):
		output_path = arg
	elif opt in ('-g'):
		gap_perc = float(arg)
	elif opt in ('-t'):
		thr_distr = int(arg)
	elif opt in ('-m'):
		my_mx = arg
	else:
		usage()
		sys.exit(2)

'''
def plotting()
	sns.set(style="darkgrid")
	rs = np.random.RandomState(8)
'''
def lookahead(iterable):
	"""Pass through all values from the given iterable, augmented by the
	information if there are more values to come after the current one
	(True), or if it is the last value (False).
	"""
	# Get an iterator and pull the first value.
	it = iter(iterable)
	last = next(it)
	# Run the iterator to exhaustion (starting from the second value).
	for val in it:
	    # Report the *previous* value (more to come).
	    yield last, True
	    last = val
	# Report the last value.
	yield last, False

def make_length_distr(df):
	'''
	Takes in dataframe with values per file and returns
	a length distribution dictionary with keys files and
	values lengths of uninterrupted positive scoring positions.
	(Can be interupted by 1 low scoring position (or more set with -t)
	This means that a sequence of ++-+-+-+-- will have a length of 5)
	'''
	length_distr={}
	weight_distr={}
	for file in df:
		i=0
		k=0
		w=0
		for pos,has_more in lookahead(df[file]):
			if pos is None:
				pass
			else:
				if pos > 1:
					i+=1
					w+=pos
				elif pos <= 1:
					if k == thr_distr:
						if file in length_distr.keys():
							length_distr[file].append(i)
							weight_distr[file].append(w)
							w=0
							i=0
							k=0
						else:
							weight_distr[file]=[]
							weight_distr[file].append(w)
							length_distr[file]=[]
							length_distr[file].append(i)
							w=0
							i=0
							k=0
					elif k < thr_distr:
						w+=pos
						#i+=1	#Think about that one
						k+=1
			if has_more is False:
				if file in length_distr.keys():
					weight_distr[file].append(w)
					length_distr[file].append(i)
					w=0
					i=0
				else:
					weight_distr[file]=[]
					length_distr[file]=[]
					weight_distr[file].append(w)
					length_distr[file].append(i)
					w=0
					i=0
	return length_distr, weight_distr

def ribbon_plot(newdict, bin_edges):
	traces = []
	xtickvals = []
	y_raw = bin_edges
	samples = sorted(list(newdict.keys()))
	sample_labels = [re.compile(r"\..*").sub("", m) for m in samples]
	for i in range(0, len(samples)):
		print(samples[i])
		xtickvals.append((i+1)*2+0.5)
		z_raw = newdict[samples[i]]
		#print(samples[i],z_raw)
		x = []
		y = []
		z = []
		for j in range(1, len(z_raw)):
			z.append([z_raw[j], z_raw[j]])
			y.append([y_raw[j], y_raw[j]])
			x.append([(i+1)*2, (i+1)*2+1])
		traces.append(dict(z=z, x=x, y=y, showscale=False, type='surface'))
	layout = go.Layout(
						title='Segment length distributions',
						scene = dict(
						xaxis=dict(
							tickmode="array", ticktext=sample_labels, tickvals=xtickvals),
						yaxis=dict(
							title="Length of uninterrupted positive scores"),
						zaxis=dict(
							title="Number of segments"))
						)
	fig = go.Figure(data=traces, layout=layout)
	plotly.offline.plot(fig, filename=output_path)

def make_hist(input_dict):
	maxlength=0
	for x in input_dict:
		if maxlength < max(input_dict[x]):
			maxlength = max(input_dict[x])
	bins = list(range(0,int(maxlength)+1))
	newdict={}
	for x in input_dict:
		hist, bin_edges = np.histogram(input_dict[x], bins=bins)
		newdict[x]=hist
	return newdict, bin_edges

def main():
	group_dict={}
	for file in os.listdir(input_path):
		#print(file)
		if re.findall(r'(.*\/)(.*)(\.fasta|\.fas)',input_path+file):
			out_dict,consnam_list=single_cons_comp.main(input_path+file,my_mx,gap_perc)
			group_dict[file]=out_dict
		else:
			raise ValueError("Directory must have only .fas or .fasta alignment files!")
	df = pd.DataFrame.from_dict(group_dict)
	length_distr, weight_distr = make_length_distr(df)
	lendict, len_bin_edges = make_hist (length_distr)
	'''
	for x in sorted(weight_distr):
		print(x, weight_distr)
	weidict, wei_bin_edges = make_hist (weight_distr)
	print(wei_bin_edges)
	for x in sorted(weidict):
		print(x, weidict[x])
	'''
	ribbon_plot(lendict, len_bin_edges)
	ax = plt.subplot()
	for k, v in lendict.items():
		plt.plot(range(1, len(v) + 1), v, '.-', label=k)
	ax.set_ylim(0,20)
	plt.legend()
	plt.savefig(output_path, dpi=600)


if __name__ == "__main__":
	main()
