#!/usr/bin/env python3
"""Calculates conservation score for multiple alignments"""
import PhyMeas, SlidingWindow
import re, os, sys, csv, getopt, plotly, single_cons_comp, argparse
import plotly.graph_objs as go
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats
#from labellines import labelLine, labelLines
from operator import itemgetter

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Calculate and visualize conservation between two groups of sequences from multiple alignments.')
	parser.add_argument('alignment_path', help='Path to folder with alignment files.')
	parser.add_argument('output_path', help='Path to image for output.')
	parser.add_argument('-t','--threshold', help='Threshold for number of allowed bad scores when calculating length of positive sections.', type=int, default=1, required=False)
	parser.add_argument('-s','--structure_path', help='Path to folder with structure files; names should match alignment groups within files.')
	parser.add_argument('-w','--window', help='Window for sliding the two groups', type=int, required=False)
	parser.add_argument('-c','--csv', help='Output length and weight distributions in a csv file', required=False, action="store_true")
	parser.add_argument('-l','--leg', help='Do not write out a legend', default=False, action="store_true")
	output_group = parser.add_mutually_exclusive_group(required=True)
	output_group.add_argument('-ps', '--scatter_plot', help='Plots a scatter of the length for positive stretches and their total score.', action="store_true")
	output_group.add_argument('-pr', '--ribbon_plot', help='Plots a 3D ribbon of the length for positive stretches and their total score.', action="store_true")
	output_group.add_argument('-pm', '--multi_plot', help='Plots a scatter of the representation for each alignment. For many (>20) alignments.', action="store_true")
	output_group.add_argument('-mm', '--multi_max', help='Plots a scatter of the max values for each alignment. Used for determining separation threshold.', action="store_true")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def make_length_distr(df,comm_args,group_dict,aln_total_lengths):
	'''Takes in dataframe with values per file and returns
	a length distribution dictionary with keys files and
	values lengths of uninterrupted positive scoring positions.
	(Can be interupted by 1 low scoring position (or more set with -t)
	This means that a sequence of ++-+-+-+-- will have a length of 5)
	'''
	length_distr={}
	weight_distr={}
	length_to_weight={}
	for file in df:
		l,i,k,w=0,0,0,0
		alignment_length = len(group_dict[file])
		for pos,has_more in PhyMeas.lookahead(df[file]):
			l+=1
			#print(l,file)
			if has_more is False:
				scaled_length = i/aln_total_lengths[file][1]
				#print(i, 'last')
				if i > 0 and w > 1:
					if file not in length_distr.keys():
						length_distr[file]=[]
						weight_distr[file]=[]
						###
						length_to_weight[file]={}
						length_to_weight[file][i]=[]
						###
					length_distr[file].append(i)
					weight_distr[file].append((scaled_length,w))
					###
					if i not in length_to_weight[file].keys():
						length_to_weight[file][i]=[]
					length_to_weight[file][i].append((w,l,scaled_length))
					###
			if pos > 1:
				#print(i, 'greater')
				k=0
				i+=1
				w+=pos
			elif pos <= 1:
				scaled_length = i/aln_total_lengths[file][1]
				if k == comm_args.threshold:
					#print(i, 'smaller k is threshold')
					if i > 0 and w > 1:
						if file not in length_distr.keys():
							length_distr[file]=[]
							weight_distr[file]=[]
							###
							length_to_weight[file]={}
							length_to_weight[file][i]=[]
							###
						length_distr[file].append(i)
						weight_distr[file].append((scaled_length,w))
						###
						if i not in length_to_weight[file].keys():
							length_to_weight[file][i]=[]
						length_to_weight[file][i].append((w,l,scaled_length))
						###
						i,k,w=0,0,0
				elif k < comm_args.threshold:
					#print(i, 'smaller k is not threshold')
					#w+=pos	#Think about these two
					#i+=1	
					k+=1
	return length_distr, weight_distr, length_to_weight

def slope(point1, point2):
	return np.arctan2(point2[1]-point1[1],point2[0]-point1[0])

def make_hist(input_dict):
	maxlength=0
	for file in input_dict:
		if maxlength < max(input_dict[file]):
			maxlength = max(input_dict[file])
	bins = list(range(0,int(maxlength)+1))
	newdict={}
	for file in input_dict:
		hist, bin_edges = np.histogram(input_dict[file], bins=bins)
		newdict[file]=hist
	return newdict, bin_edges

def ribbon_plot(newdict, bin_edges,output_path):
	traces, xtickvals = [], []
	y_raw = bin_edges
	samples = sorted(list(newdict.keys()))
	sample_labels = [re.compile(r"\..*").sub("", m) for m in samples]
	for i in range(0, len(samples)):
		xtickvals.append((i+1)*2+0.5)
		z_raw = newdict[samples[i]]
		#print(samples[i],z_raw)
		#print(samples[i],sum(z_raw*range(len(z_raw)))/sum(z_raw))
		x,y,z = [],[],[]
		for j in range(1, len(z_raw)):
			z.append([z_raw[j], z_raw[j]])
			y.append([y_raw[j], y_raw[j]])
			x.append([(i+1)*2, (i+1)*2+1])
		traces.append(dict(z=z, x=x, y=y, showscale=False, type='surface'))
	layout = go.Layout(	title='Segment length distributions',
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

def scatter_plot(comm_args,weight_distr,length_distr,_lines=False,_maxx=False,_scatter=False):
	'''Outputs scatter plot image with colors depending on hte number of the input alignments.
	'''
	###   Defining Color scheme   ###
	ax = plt.subplot()
	if len(weight_distr) == 20:
		colors = matplotlib.cm.tab20(np.linspace(0, 1, len(weight_distr)))
	elif len(weight_distr) == 10:
		#In cases of just 10 alignments it assumes two sets of 5 each and uses the seismic/divergent gradient.
		#colors = matplotlib.cm.seismic(np.linspace(0, 1, len(weight_distr)))
		colors = matplotlib.cm.PRGn(np.linspace(0, 1, len(weight_distr)))
	elif len(weight_distr) > 20 and 'A_' in weight_distr.keys():
		#Creates a color mapping for groups of alignments defined with prepended A_, B_, C_ and so on.
		#Used in cases where there are too many alignments to properly discriminate with colors.
		#Uses the viridis gradient.
		colors=[]
		letters=[]
		for letter in sorted(weight_distr.keys()):
			letters.append(letter.split('_')[0])

		color_set = matplotlib.cm.viridis(np.linspace(0, 1, len(set(letters))))
		letter_color_map = dict(zip(sorted(set(letters)),color_set))

		for letter in letters:
			colors.append(letter_color_map[letter])
	else:
		colors = matplotlib.cm.tab20(np.linspace(0, 1, len(weight_distr)))

	###   Plotting   ###
	sorted_names = sorted(weight_distr.keys())
	for file, color in zip(sorted_names,colors):
		if _lines:
			degree_label = 180-round(np.rad2deg(slope(weight_distr[file][0],weight_distr[file][1])),2)
			plt.scatter(*zip(*weight_distr[file]), label=str(degree_label)+' '+re.sub(r'\.fas.*','',file), marker='.',color=color)
			plt.plot(*zip(*weight_distr[file]),color=color)
		if _scatter:
			abs_length = [n**2 for n in length_distr[file]]
			print(file, max(weight_distr[file],key=itemgetter(0)))
			plt.scatter(*zip(*weight_distr[file]), label=re.sub(r'\.fas.*','',file),marker='.',s=abs_length,color=color)
		if _maxx:
			plotlist=[]
			for x in zip(weight_distr[file].keys(),weight_distr[file].values()):
				plotlist.append((x[0],max(x[1],key=itemgetter(0))[0]))
			plt.scatter(*zip(*plotlist), label=re.sub(r'\.fas.*','',file), marker='.',color=color)
	#labelLines(plt.gca().get_lines(),align=False)
	#plt.plot((1,13),(4.5,0),color='black')
	if not comm_args.leg:
		lgnd = plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
		for n in range(len(weight_distr)):
			lgnd.legendHandles[n]._sizes = [30]
	plt.savefig(comm_args.output_path, dpi=600, bbox_inches='tight')
	return True

def csv_output(comm_args, file_to_data, weight_distr):
	'''
	Writes out data used for generating the plot in a csv file.
	'''
	#Assuming output path is always written with a .png at the end...
	with open(re.sub(r'.png','.csv',comm_args.output_path), mode ="w") as output_csv:
		csv_writer = csv.writer(output_csv)
		csv_writer.writerow(['File name', 'Segment length', 'Normalized segment length','Total segment weight', 'Alignment position'])
		for file, length_buckets in sorted(file_to_data.items()):
			for length, weights in sorted(length_buckets.items()):
				for single_weight in sorted(weights,key=itemgetter(0)): 
					csv_writer.writerow([file, length, single_weight[2],single_weight[0], single_weight[1]])
	return True


def main(commandline_args):
	comm_args = create_and_parse_argument_options(commandline_args)
	if comm_args.window:
		if not os.path.isfile(comm_args.alignment_path):
			raise ValueError("In case of specified window for sliding (-w argument), the alignment path  must be a single file!")
		SlidingWindow.main([comm_args.alignment_path,'-w '+str(comm_args.window)])
		sys.exit()

	group_dict={}
	aln_lengths={}
	for file in os.listdir(comm_args.alignment_path):
		if re.findall(r'(.*\/)(.*)(\.fasta|\.fas|\.fa)',comm_args.alignment_path+file):
			# print(file)
			out_dict={}
			alnindex_score,sliced_alns,number_of_aligned_positions=PhyMeas.main(['-a',comm_args.alignment_path+file, '-r', '-bl'])
			for x in sliced_alns:
				aln_lengths[file]=(sliced_alns[x].get_alignment_length(),number_of_aligned_positions)
				break
			for x in alnindex_score.keys():
				out_dict[x] = alnindex_score[x][0]
			group_dict[file] = out_dict
		else:
			raise ValueError("Directory must have only .fa, .fas or .fasta alignment files!")
	
	df = pd.DataFrame.from_dict(group_dict)
	length_distr, weight_distr, length_to_weight = make_length_distr(df,comm_args,group_dict,aln_lengths)

	if comm_args.ribbon_plot:
		lendict, len_bin_edges = make_hist (length_distr)
		ribbon_plot(lendict, len_bin_edges,comm_args.output_path)
		#weidict, wei_bin_edges = make_hist (weight_distr)
		#ribbon_plot(weidict, wei_bin_edges,comm_args.output_path)
	elif comm_args.scatter_plot:
		scatter_plot(comm_args,weight_distr,length_distr,_scatter=True)
	elif comm_args.multi_plot:
		scatter_plot(comm_args,length_to_weight,length_distr,_maxx=True)

		#for file in weight_distr.keys():
		#	maxweight[file]=[max(weight_distr[file],key=itemgetter(0)),max(weight_distr[file],key=itemgetter(1))]
		#	#print(file, slope(maxweight[file][0],maxweight[file][1]))
		#scatter_plot(comm_args,maxweight,length_distr,_lines=True)
	
	elif comm_args.multi_max:
		#Grouping max values
		max_weight_bylength={}
		csv_output(comm_args, length_to_weight, weight_distr)
		for file in length_to_weight.keys():
			plotlist=[]
			for x in zip(length_to_weight[file].keys(),length_to_weight[file].values()):
				plotlist.append((x[0],max(x[1],key=itemgetter(0))[0]))
				if x[0] in max_weight_bylength.keys():
					if max(x[1],key=itemgetter(0))[0] > max_weight_bylength[x[0]]:
						max_weight_bylength[x[0]]=max(x[1],key=itemgetter(0))[0]
					else:
						pass
				else:
					max_weight_bylength[x[0]]=max(x[1],key=itemgetter(0))[0]
		
		#Fitting line on max values
		#print(max_weight_bylength)
		#print(max_weight_bylength.keys(),max_weight_bylength.values())
		slope, intercept, r_value, p_value, std_err = stats.linregress(list(max_weight_bylength.keys()),list(max_weight_bylength.values()))
		print(slope, intercept, r_value, p_value, std_err)
		axes = plt.gca()
		plt.plot(list(max_weight_bylength.keys()),list(max_weight_bylength.values()), 'o', label='Maximal scores from all alignments')
		x_vals = np.array(axes.get_xlim())
		y_vals = intercept + slope * x_vals
		label_line = "Line fit with p value "+str(p_value)
		plt.plot(x_vals, y_vals, 'r', label=label_line)
		plt.legend()
		plt.savefig(comm_args.output_path, dpi=600, bbox_inches='tight')
	
	if comm_args.csv:
		csv_output(comm_args, length_to_weight, weight_distr)

	

if __name__ == "__main__":
	main(sys.argv[1:])
