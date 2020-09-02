#!/usr/bin/env python3
"""Calculate and visualize conservation between two groups of sequences from one alignment"""
import re, os, csv, sys, random, Bio.Align, argparse, random, math, matplotlib, ntpath
sys.path.append(os.path.dirname(os.path.abspath(__name__)))
matplotlib.use('Agg')
import numpy as np
from datetime import date
from Bio import AlignIO
from io import StringIO
from textwrap import wrap
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from Bio.SeqUtils import IUPACData
from bin.AlignmentGroup import AlignmentGroup
import bin.Sequence_Weight_from_Tree as Sequence_Weight_from_Tree
from bin.MatrixLoad import PAMLmatrix
from Bio.SubsMat import MatrixInfo

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o','--output_path', help='Output path', required=True)
    input_file = parser.add_mutually_exclusive_group(required=True)
    input_file.add_argument('-a','--alignment_paths', nargs='+', help='Path to alignment files. If given two files it will use mafft --merge to merge them in single alignment.', action=required_length(1,2))
    input_file.add_argument('-as','--alignment_string', help='Alignment string', type=str)
    parser.add_argument('-cg','--cut_gaps', help='Remove alignment positions with %% gaps greater than the specified value with gap_threshold.', action="store_true")
    parser.add_argument('-gt','--gap_threshold', help='Specify %% gaps per alignment position. (Default = the smaller between ((sequences of group1)/(all sequences) and (sequences of group2)/(all sequences)) minus 0.05)', type=float)
    parser.add_argument('-s','--structure_paths', nargs='+', help='Paths to structure files, for score calculation. Does not work with --nucleotide!')
    parser.add_argument('-sy','--structure_pymol', nargs='+', help='Paths to structure files, for plotting a pml.')
    parser.add_argument('-phy','--phylo_split', help='Split the alignment in two groups by constructing a tree instead of looking for _ separated strings.', action="store_true")
    parser.add_argument('-nc','--nucleotide', help='Use nucleotide matrix for score calculation', choices=['blastn', 'identity', 'trans'])
    parser.add_argument('-w','--weigh_sequences', help='Weigh sequences within each alignment group.', choices=['pairwise', 'voronoi'])
    output_type_group = parser.add_mutually_exclusive_group(required=True)
    output_type_group.add_argument('-p', '--plotit', help='Plots the calculated score as a bar graph for each alignment position.', action="store_true")
    output_type_group.add_argument('-pml', '--write_pml_script', help='Writes out a PyMOL coloring script for any structure files that have been defined. Choose between unix or windows style paths for the pymol script.', choices=['unix', 'windows'])
    output_type_group.add_argument('-r', '--return_within', help='To be used from within other python programs. Returns dictionary of alnpos->score.', action="store_true")
    output_type_group.add_argument('-csv', '--return_csv', help='Saves a csv with alignment position -> score.', action="store_true")
    output_type_group.add_argument('-jv', '--jalview_output', help='Saves an annotation file for Jalview.', action="store_true")
    entropy_group = parser.add_mutually_exclusive_group()
    entropy_group.add_argument('-mx','--substitution_matrix', help='Choose protein substitution matrix for score calculation.', choices=MatrixInfo.available_matrices)
    entropy_group.add_argument('-lg','--leegascuel', help='Use LG matrix for score calculation', action="store_true")
    entropy_group.add_argument('-e','--shannon_entropy', help='Use shannon entropy for conservation calculation.', action="store_true")
    entropy_group.add_argument('-rs','--reflected_shannon', help='Use shannon entropy for conservation calculation and reflect the result so that a fully random sequence will be scored as 0.', action="store_true")
    structure_option = parser.add_mutually_exclusive_group()
    structure_option.add_argument('-ss','--secondary_structure', help = 'Use substitution matrices derived from data dependent on the secondary structure assignment.', action="store_true")
    structure_option.add_argument('-be','--burried_exposed', help = 'Use substitution matrices derived from data dependent on the solvent accessability of a residue.', action="store_true")
    structure_option.add_argument('-ssbe','--both', help = 'Use substitution matrices derived from data dependent on both the secondary structure and the solvent accessability of a residue.', action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def required_length(nmin,nmax):
    '''Limiter for passed arguments.
    '''
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                    f=self.dest,nmin=nmin,nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength

def read_align(aln_path):
    '''Reads the fasta file and gets the sequences.
    '''
    alignments = AlignIO.read(open(aln_path), "fasta")
    for record in alignments:
        record.seq = record.seq.upper()
    return alignments

def deletefile(file_loc):
    '''Tries to delete provided file path.
    '''
    import subprocess
    try:
        subprocess.run(['rm', file_loc], check = True)
    except subprocess.CalledProcessError:
        raise IOError("When using mafft for merging two alignments working directory must be writable!")

def run_mafft(aln_paths):
    '''Tags separate alignments for TwinCons and merges them with mafft --merge.
    '''
    import warnings
    tempfiles = ['./tempsubMSAtable', './tempconcatfasta.fas', './tempmergedfasta.fas']
    for tempfile in tempfiles:
        if os.path.isfile(tempfile):
            warnings.warn(f"When using mafft for merging two alignments working directory must be free of file {tempfile}. Trying to delete the file.")
            deletefile(tempfile)
    list_with_alns = [read_align(aln_path) for aln_path in aln_paths]
    concatedfasta_handle = open("./tempconcatfasta.fas", "a")
    mergertable_ix = list()
    previous_len = 0
    for i, aln in enumerate(list_with_alns, 1):
        mergertable_ix.append((aln_paths[i-1],len(aln),previous_len))
        previous_len = len(aln)
        for seq in aln:
            seq.id = str(i)+"_"+seq.id
        AlignIO.write(aln, concatedfasta_handle, "fasta")
    concatedfasta_handle.close()
    mergertable = open("./tempsubMSAtable", "a")
    for aln_len in mergertable_ix:
        seq_list = [str(x) for x in list(range(aln_len[2]+1,aln_len[2]+aln_len[1]+1))]
        seq_nums = " ".join(seq_list)+" #"+aln_len[0]+"\n"
        mergertable.write(" "+seq_nums)
    mergertable.close()
    os.system("mafft --quiet --merge ./tempsubMSAtable ./tempconcatfasta.fas > ./tempmergedfasta.fas")
    merged_tagged_aln = read_align("./tempmergedfasta.fas")
    for tempfile in tempfiles:
        deletefile(tempfile)
    return merged_tagged_aln

def count_aligned_positions(aln_obj, gap_threshold):
    '''Counts how many positions are aligned (less than gap_threshold gaps)
    '''
    number_seqs = len(aln_obj)
    aligned_positions = 0
    extremely_gapped = dict()
    for i in range(0,aln_obj.get_alignment_length()):
        extremely_gapped[i+1] = 'True'
        if aln_obj[:,i].count('-')/number_seqs <= float(gap_threshold):
            aligned_positions+=1
            extremely_gapped[i+1] = 'False'
    if aligned_positions == 0:
        raise ValueError('Alignment:\n'+str(aln_obj)+'\nhas no positions with less than '+str(gap_threshold*100)+'% gaps!')
    return aligned_positions, extremely_gapped

def remove_extremely_gapped_regions(align, gap_perc, gap_mapping):
    '''Removes columns of alignment with more than gap_perc gaps.
    '''
    n = float(len(align[0]))
    i, x = 0, 0
    length=1
    while i < n:
        x = align[:, i].count('-')/len(align)                 #Get percentage of gaps in column
        if float(x) > abs(float(gap_perc)):
            if i == 0:
                align = align[:, 1:]
            elif i+1 == n:
                align = align[:, :i]
            else:
                align = align[:, :i] + align[:, i+1:]
            n -= 1                                            #  seq. 1 shorter
        else:                                                 #  nothing to delete, proceed
            i += 1
        length+=1
        if i in gap_mapping.keys():
            continue
        gap_mapping[i] = int(length)-1
    return gap_mapping, align, len(align[0])

def uniq_resi_list(aln_obj):
    '''
    Returns list of unique AA residues in the given MSA to be used for frequency iterator.
    Checks if the alignment has AA letters from the IUPAC extended_protein_letters.
    '''
    default_aa_sequence = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    hash_resi=dict()
    for alignment in aln_obj:
        for resi in alignment.seq:
            if not re.match(r'-|X',resi):
                hash_resi[resi]='null'
    if not all (x in IUPACData.extended_protein_letters for x in hash_resi.keys()):
        raise ValueError("Alignment:\n"+str(aln_obj)+"\nhas AA letters not found in the IUPAC extended list!")
    if len(Counter(hash_resi.keys())) > len(Counter(default_aa_sequence)):
        raise ValueError("Alignment has non-standard AAs:\n"+' ,'.join(hash_resi.keys()))
    return default_aa_sequence

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
        what = Bio.Align.MultipleSeqAlignment([])
        for entry in unsliced_aln_obj:
            if re.match(prot,entry.id.split("_")[0]):
                what.append(entry)
        sliced_dict[prot]=what
    return sliced_dict            #Iterate over the dict and create instances of AlignmentGroup

def determine_subs_matrix(comm_args):
    if comm_args.nucleotide:
        mx = nucl_matrix(comm_args.nucleotide)
    if comm_args.leegascuel or comm_args.structure_paths:
        mx = np.array(PAMLmatrix(str(os.path.dirname(__file__))+'/../matrices/LG.dat').lodd)
    if comm_args.substitution_matrix:
        mx = subs_matrix(comm_args.substitution_matrix)
    return mx, mx.min(), mx.max()

def gradientbars(bars, positivegradient, negativegradient, mx_min, mx_max):
        ax = bars[0].axes
        lim = ax.get_xlim()+ax.get_ylim()
        for bar in bars:
            bar.set_zorder(1)
            bar.set_facecolor("none")
            x,y = bar.get_xy()
            w, h = bar.get_width(), bar.get_height()
            if h >= 0:
                grad = np.atleast_2d(np.linspace(0,h/mx_max,256)).T
                ax.imshow(grad, extent=[x,x+w,y,y+h], cmap=plt.get_cmap(positivegradient), aspect="auto", norm=matplotlib.colors.NoNorm(vmin=0,vmax=1))
            else:            #different gradient for negative values
                grad = np.atleast_2d(np.linspace(0,h/mx_min,256)).T
                ax.imshow(grad, extent=[x,x+w,y,y+h], cmap=plt.get_cmap(negativegradient), aspect="auto", norm=matplotlib.colors.NoNorm(vmin=0,vmax=1))
        #ax.set_facecolor('Gray')
        ax.axis(lim)

def upsidedown_horizontal_gradient_bar(out_dict,group_names,comm_args, mx_min, mx_max):
    plt.rcParams['image.composite_image'] = False                    #So that bars are separate images
    plt.rcParams["figure.figsize"] = (8,8)
    fig, ax = plt.subplots()
    data=[]
    stdevdata=[]
    for x in sorted(out_dict.keys()):
        data.append(out_dict[x][0])
        if out_dict[x][1] == 'True':
            stdevdata.append(0.5)
        else:
            stdevdata.append(0)
    bar = ax.bar(range(1,len(data)+1),data, yerr=stdevdata,error_kw=dict(ecolor='gray', lw=0.25))

    #In case of no negative values BUG!
    if comm_args.reflected_shannon or comm_args.shannon_entropy:
        plt.yticks(np.arange(0,4.2, step=0.5))
        gradientbars(bar,'viridis','binary', mx_min, mx_max)
    else:
        plt.plot((0, len(data)+1), (1, 1), 'k-', linewidth=0.5)       #Horizontal line
        gradientbars(bar,'Greens','Purples', mx_min, mx_max)
    dpi_scaling = 3*len(out_dict)
    plt.savefig(comm_args.output_path+'.svg',format = 'svg',dpi=dpi_scaling)
    return True

def data_to_diverging_gradients(datapoint, maxdata, mindata, positivegradient, negativegradient):
    '''Maps a data point to a diverging colormap depending on whether its above or bellow 0.
    Returns a hex code.
    '''
    if datapoint == 'NA':
        return '#808080'
    if datapoint >= 0:
        grad = np.atleast_2d(np.linspace(0,datapoint/maxdata,256)).T
        rgb = plt.get_cmap(positivegradient)(grad[len(grad)-1])[0][:3]
    else:
        grad = np.atleast_2d(np.linspace(0,datapoint/mindata,256)).T
        rgb = plt.get_cmap(negativegradient)(grad[len(grad)-1])[0][:3]
    return matplotlib.colors.rgb2hex(rgb)

def gradients(data, positivegradient, negativegradient, mx_maxval, mx_minval):
    """Creates a dictionary with keys the alignment index and values a hex code for the color.
    """
    aln_index_hexcmap = {}
    aln_index=1
    for datapoint in data:
        aln_index_hexcmap[aln_index] = data_to_diverging_gradients(datapoint,
                                                                    mx_maxval,
                                                                    mx_minval,
                                                                    positivegradient,
                                                                    negativegradient)
        aln_index+=1
    return aln_index_hexcmap

def pymol_script_writer(out_dict, gapped_sliced_alns, comm_args, mx_minval, mx_maxval):
    """Creates the same gradients used for svg output and writes out a .pml file for PyMOL visualization.
    """
    from pathlib import Path, PureWindowsPath, PurePosixPath
    
    data = []
    for x in sorted(out_dict.keys()):
        data.append(out_dict[x][0])
    
    if comm_args.leegascuel or comm_args.nucleotide or comm_args.substitution_matrix:
        #alnindex_to_hexcolors = gradientbars(data,'Blues','Reds')
        alnindex_to_hexcolors = gradients(data,'Greens','Purples', mx_maxval, mx_minval)
    elif comm_args.reflected_shannon or comm_args.shannon_entropy:
        alnindex_to_hexcolors = gradients(data,'viridis','binary', mx_maxval, mx_minval)
    else:
        alnindex_to_hexcolors = gradients(data,'Greens','Purples', mx_maxval, mx_minval)

    group_names = list(gapped_sliced_alns.keys())
    #Open .pml file for structure coloring
    pml_output = open(comm_args.output_path+".pml","w")
    #pml_output = open("./"+'-'.join(sorted(group_names))+'.pml',"w")
    pml_output.write("set hash_max, 500\n\
        set cartoon_loop_radius,0.4\n\
        set cartoon_tube_radius,1\n\
        set cartoon_ladder_radius,0.6\n\
        set cartoon_oval_length,1.4\n\
        set cartoon_oval_width,0.6\n\
        set ray_opaque_background, off\n\
        bg_color black\n\
        set ray_trace_mode,1\n\
        set ray_shadows,0\n")

    #Bellow here needs fixing to properly do structures for plotting
    for alngroup_name in group_names:
        #Match groupnames with structure files
        current_path = [s for s in comm_args.structure_pymol if alngroup_name in s]
        
        if len(current_path) == 0:
            raise IOError("Cannot write PyMOL coloring script without at least single matching structure \
               and sequence!\nSequence:\t"+alngroup_name+"\nStructure:\t"+str(current_path))
        else:
            #We have to recalculate the structure to alignment mapping
            alngroup_name_object = AlignmentGroup(gapped_sliced_alns[alngroup_name],current_path[0])
            AlignmentGroup.add_struc_path(alngroup_name_object, current_path[0])
            struc_to_aln_index_mapping=AlignmentGroup.create_aln_struc_mapping_with_mafft(alngroup_name_object)
            #Open the structure file
            if comm_args.write_pml_script == 'unix':
                pml_path = PurePosixPath(current_path[0].replace(ntpath.dirname(comm_args.output_path), '')[1:])
            if comm_args.write_pml_script == 'windows':
                pml_path = PureWindowsPath(current_path[0].replace(ntpath.dirname(comm_args.output_path), '')[1:])
            pml_output.write(f"load {pml_path}, {alngroup_name}\n")
            #For each alignment position, color the appropriate residue with the hex transformed color from the gradient
            for aln_index in alnindex_to_hexcolors.keys():
                if aln_index in struc_to_aln_index_mapping:
                    hexcolors_appropriate_for_pml = alnindex_to_hexcolors[aln_index].replace('#','0x')
                    pml_output.write(f"color {hexcolors_appropriate_for_pml}, {alngroup_name} and resi {str(struc_to_aln_index_mapping[aln_index])}\n")
    pml_output.write(f"super {group_names[0]}, {group_names[1]}\n")
    return True

def jalview_output(output_dict, comm_args):
    '''Writes out a Jalview alignment annotation file (http://www.jalview.org/help/html/features/annotationsFormat.html)
    '''
    out_data=list()
    for x in sorted(output_dict.keys()):
        out_data.append(output_dict[x][0])
    
    jv_output = open(comm_args.output_path+".jlv","w")
    jv_output.write('JALVIEW_ANNOTATION\n')
    jv_output.write('# Created: '+str(date.today())+"\n")
    jv_output.write('# Contact: ppenev@gatech.edu\n')
    jv_output.write('BAR_GRAPH\tTWINCONS\t')
    for position in sorted(output_dict.keys(), key=abs):
        color_hex = data_to_diverging_gradients(output_dict[position][0], max(out_data), min(out_data), 'Greens', 'Purples')
        jv_output.write(str(output_dict[position][0])+'['+str(color_hex).replace('#','')+']|')
    return True

def shannon_entropy(alnObject, aa_list, commandline_args, alngroup_to_sequence_weight):
    '''
    Function to calcuate the reflected Shannon entropy per alignment column.
    Fully random sequence will have reflected entropy of 0, while fully conserved
    column will be around 4.322
    H=-sum_{i=1}^{M} P_i,log_2,P_i    (http://imed.med.ucm.es/Tools/svs_help.html)
    Hrefl = abs(log_2,M) +sum_{i=1}^{M} P_i,log_2,P_i
    For gapped regions a random selection is made from the present AAs.
    '''
    entropy_alignment = AlignmentGroup(alnObject)
    sequence_weights = dict()
    if len(alngroup_to_sequence_weight) != 0:
        sequence_weights = alngroup_to_sequence_weight['shannon']
    frequency_distribution = entropy_alignment.column_distribution_calculation(aa_list, 
                                                                               alnObject.get_alignment_length(),
                                                                               sequence_weights)
    
    entropyDict = dict()
    for column_ix in sorted(frequency_distribution.keys(), key=abs):
        entropy_list = list()
        for P_i in frequency_distribution[column_ix]:
            entropy_i = 0
            if P_i != 0:
                entropy_i = P_i*(math.log(P_i,2))
            entropy_list.append(entropy_i)
        if commandline_args.reflected_shannon:
            refl_shentr = abs(math.log(1/len(aa_list),2))+sum(entropy_list)
            entropyDict[column_ix]=refl_shentr
        elif commandline_args.shannon_entropy:
            sh_entropy = abs(sum(entropy_list))
            entropyDict[column_ix]=sh_entropy
    return entropyDict

def nucl_matrix(mx_def):
    '''Return a substitution matrix for nucleotides.
    '''
    nucl_sequence = ['A','U','C','G']
    if mx_def == 'identity':
        nuc_mx = np.array([[7,-5,-5,-5],[-5,7,-5,-5],[-5,-5,7,-5],[-5,-5,-5,7],])
    elif mx_def == 'blastn':
        nuc_mx = np.array([[5,-4,-4,-4],[-4,5,-4,-4],[-4,-4,5,-4],[-4,-4,-4,5],])
    elif mx_def == 'trans':
        nuc_mx = np.array([[6,-5,-5,-1],[-5,6,-1,-5],[-5,-1,6,-5],[-1,-5,-5,6],])

    testvr = np.repeat(1/len(nucl_sequence),len(nucl_sequence))
    baseline = float(testvr@np.array(nuc_mx)@testvr.T)
    revtestA=np.add(np.array(nuc_mx), abs(baseline))
    if int(testvr@revtestA@testvr.T) != 0:
        raise ValueError("Wasn't able to baseline the substitution matrix correctly!")
    return np.add(np.array(nuc_mx),abs(baseline))

def subs_matrix(matrix):
    '''Baseline and return a numpy form of the substitution matrix.
    '''
    aa_sequence = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    loddmx = []
    substitution_matrix = getattr(MatrixInfo,matrix)
    for aa in aa_sequence:
        linemx=[]
        for aa2 in aa_sequence:
            if (aa,aa2) in substitution_matrix:
                linemx.append(substitution_matrix[(aa,aa2)])
            else:
                linemx.append(substitution_matrix[(aa2,aa)])
        loddmx.append(linemx)
    testvr = np.repeat(1/len(aa_sequence),len(aa_sequence))
    baseline = float(testvr@np.array(loddmx)@testvr.T)
    revtestA=np.add(np.array(loddmx), abs(baseline))
    if int(testvr@revtestA@testvr.T) != 0:
        raise ValueError("Wasn't able to baseline the substitution matrix correctly!")
    return np.add(np.array(loddmx),abs(baseline))

def struc_anno_matrices (struc_anno):
    '''Returns a log odds matrix from a given name of a PAML type matrix'''
    return np.array(PAMLmatrix(str(os.path.dirname(__file__))+'/../matrices/'+struc_anno+'.dat').lodd)

def compute_score(aln_index_dict, groupnames, mx=None, struc_annotation=None):
    '''
    Computes transformation score between two groups, using substitution 
    matrices on the common structural elements between two groups.
    Returns a dictionary with key the alignment index and value the computed score.
    '''
    if struc_annotation and mx:
        raise IOError("Do not use structure defined matrices and sequence based matrices at the same time.")
    alnindex_score = defaultdict(dict)
    for aln_index in aln_index_dict:
        vr1 = np.array(aln_index_dict[aln_index][groupnames[0]])
        vr2 = np.array(aln_index_dict[aln_index][groupnames[1]])
        if struc_annotation:
            mx = np.array(PAMLmatrix(str(os.path.dirname(__file__))+'/../matrices/LG.dat').lodd)
            if aln_index in struc_annotation[groupnames[0]] and aln_index in struc_annotation[groupnames[1]]:
                common_chars = sorted(set(struc_annotation[groupnames[0]][aln_index]) & set (struc_annotation[groupnames[1]][aln_index]))
                if len(common_chars) > 0:
                    mx = struc_anno_matrices(''.join(common_chars))
        alnindex_score[aln_index] = vr1@mx@vr2.T
    return alnindex_score

def decision_maker(comm_args, alignIO_out, sliced_alns, aa_list, alngroup_to_sequence_weight, mx):
    '''Checks through the commandline options and does the appropriate frequency and score calculations.
    Returns a dictionary of alignment position -> computed score.
    '''

    if comm_args.shannon_entropy or comm_args.reflected_shannon:
        return shannon_entropy(alignIO_out, aa_list, comm_args, alngroup_to_sequence_weight)
    
    aln_index_dict = defaultdict(dict)
    struc_annotation = defaultdict(dict)
    for alngroup_name in sliced_alns:
        alngroup_name_object = AlignmentGroup(sliced_alns[alngroup_name])
        if comm_args.structure_paths:
            current_path = [s for s in comm_args.structure_paths if alngroup_name in ntpath.basename(s)]
            if len(current_path) == 0:
                raise IOError(f"When using structure-defined matrices the provided structure files must contain the name of the alignment group.\n\
                    The alignment group {alngroup_name} was not found in the structure file names {' '.join(comm_args.structure_paths)}")
            AlignmentGroup.add_struc_path(alngroup_name_object, current_path[0])
            struc_to_aln_index_mapping = AlignmentGroup.create_aln_struc_mapping_with_mafft(alngroup_name_object)
            if comm_args.secondary_structure:
                struc_annotation[alngroup_name] = AlignmentGroup.ss_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
            elif comm_args.burried_exposed:
                struc_annotation[alngroup_name] = AlignmentGroup.depth_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
            elif comm_args.both:
                struc_annotation[alngroup_name] = AlignmentGroup.both_map_creator(alngroup_name_object,struc_to_aln_index_mapping)
            else:
                raise IOError("When a structure is defined, one of the matrix options are required!")
        alnindex_col_distr = AlignmentGroup.column_distribution_calculation(alngroup_name_object,
                                                                            aa_list,
                                                                            alignIO_out.get_alignment_length(), 
                                                                            alngroup_to_sequence_weight[alngroup_name])
        for aln_index in alnindex_col_distr:
            aln_index_dict[aln_index][alngroup_name]=alnindex_col_distr[aln_index]
    if comm_args.structure_paths:
        return compute_score(aln_index_dict, list(struc_annotation.keys()), struc_annotation=struc_annotation)
    if comm_args.nucleotide:
        return compute_score(aln_index_dict, list(sliced_alns.keys()), mx=mx)
    if comm_args.leegascuel:
        return compute_score(aln_index_dict, list(sliced_alns.keys()), mx=mx)
    if comm_args.substitution_matrix:
        return compute_score(aln_index_dict, list(sliced_alns.keys()), mx=mx)

def main(commandline_arguments):
    '''Main entry point'''
    comm_args = create_and_parse_argument_options(commandline_arguments)
    if comm_args.cut_gaps and (comm_args.structure_pymol or comm_args.structure_paths):
        raise IOError("TwinCons can not take in this combination of arguments!\
    \nCombining gap removal (-cg) and structural mapping (-sy) or structure based matrices (-s) produces inconsistent alignment mapping!")
    if comm_args.alignment_string:
        alignIO_out_gapped = list(AlignIO.parse(StringIO(comm_args.alignment_string), 'fasta'))[0]
    elif len(comm_args.alignment_paths) == 1:
        alignIO_out_gapped=read_align(comm_args.alignment_paths[0])
    elif len(comm_args.alignment_paths) == 2:
        alignIO_out_gapped = run_mafft(comm_args.alignment_paths)
    randindex_norm = defaultdict(dict)
    gp_mapping = dict()

    subs_matrix, mx_minval, mx_maxval = determine_subs_matrix(comm_args)

    if comm_args.phylo_split:
        tree = Sequence_Weight_from_Tree.tree_construct(alignIO_out_gapped)
        deepestanc_to_child = Sequence_Weight_from_Tree.find_deepest_ancestors(tree)
        gapped_sliced_alns = Sequence_Weight_from_Tree.slice_by_anc(alignIO_out_gapped, deepestanc_to_child)
    else:
        deepestanc_to_child = {}
        gapped_sliced_alns = slice_by_name(alignIO_out_gapped)

    if len(gapped_sliced_alns.keys()) != 2:
        raise ValueError("For now does not support more than two groups! Offending groups are "+str(gapped_sliced_alns.keys()))

    if comm_args.gap_threshold is None:
        seq_records = list()
        for aln in gapped_sliced_alns:
            seq_records.append(gapped_sliced_alns[aln].__len__())
        comm_args.gap_threshold = round(min([seq_records[0]/(seq_records[0]+seq_records[1]),seq_records[1]/(seq_records[0]+seq_records[1])])-0.05,2)
    
    number_of_aligned_positions, extremely_gapped = count_aligned_positions(alignIO_out_gapped, comm_args.gap_threshold)
    if comm_args.cut_gaps:
        tempaln = alignIO_out_gapped[:,:]
        alignIO_out_gapped = Bio.Align.MultipleSeqAlignment([])
        gp_mapping, alignIO_out_gapped, alen = remove_extremely_gapped_regions(tempaln, float(comm_args.gap_threshold), gp_mapping)
    else:
        for i in range(1, alignIO_out_gapped.get_alignment_length()+1):
            gp_mapping[i] = i

    
    alngroup_to_sequence_weight = dict()
    for alngroup in gapped_sliced_alns:
        alngroup_to_sequence_weight[alngroup] = list()
        alngroup_to_sequence_weight['shannon'] = list()
        if comm_args.weigh_sequences:
            if comm_args.reflected_shannon or comm_args.shannon_entropy:
                alngroup_to_sequence_weight['shannon'] = Sequence_Weight_from_Tree.calculate_weight_vector(alignIO_out_gapped, algorithm=comm_args.weigh_sequences)
            else:
                alngroup_to_sequence_weight[alngroup] = Sequence_Weight_from_Tree.calculate_weight_vector(gapped_sliced_alns[alngroup], algorithm=comm_args.weigh_sequences)

    uniq_resis = uniq_resi_list(alignIO_out_gapped)
    if comm_args.nucleotide:
        for sequence in alignIO_out_gapped:
            if re.search('T', str(sequence.seq)):
                sequence.seq = sequence.seq.transcribe()
        uniq_resis = ['A','U','C','G']
    
    if comm_args.phylo_split:
        sliced_alns = Sequence_Weight_from_Tree.slice_by_anc(alignIO_out_gapped, deepestanc_to_child)
    else:
        sliced_alns = slice_by_name(alignIO_out_gapped)
    if len(sliced_alns.keys()) != 2:
        raise IOError("For now does not support more than two groups! Offending groups are "+str(sliced_alns.keys()))
    
    
    
    position_defined_scores = decision_maker(comm_args, 
                                            alignIO_out_gapped, 
                                            sliced_alns, 
                                            uniq_resis, 
                                            alngroup_to_sequence_weight, 
                                            subs_matrix)
    output_dict = dict()
    output_dict_pml = dict()
    for x in position_defined_scores.keys():                #If standard deviation is too big, set the result as 0
        if extremely_gapped[gp_mapping[x]] == 'True':
            output_dict[gp_mapping[x]] = (position_defined_scores[x], extremely_gapped[gp_mapping[x]])
            output_dict_pml[gp_mapping[x]] = ('NA', extremely_gapped[gp_mapping[x]])
            continue
        output_dict[gp_mapping[x]] = (position_defined_scores[x], extremely_gapped[gp_mapping[x]])
        output_dict_pml[gp_mapping[x]] = (position_defined_scores[x], extremely_gapped[gp_mapping[x]])
    
    if comm_args.plotit:                                    #for plotting
        upsidedown_horizontal_gradient_bar(output_dict, list(gapped_sliced_alns.keys()), comm_args, mx_minval, mx_maxval)
    elif comm_args.write_pml_script:
        pymol_script_writer(output_dict_pml, gapped_sliced_alns, comm_args, mx_minval, mx_maxval)
    elif comm_args.return_within:
        return output_dict, gapped_sliced_alns, number_of_aligned_positions, gp_mapping
    elif comm_args.return_csv:
        with open(comm_args.output_path+".csv", mode='w') as output_csv:
            csv_writer = csv.writer(output_csv, delimiter=',')
            csv_writer.writerow(["Alignment index", "Score", f"More than {comm_args.gap_threshold*100}% gaps"])
            for x in sorted(output_dict.keys(), key=abs):
                csv_writer.writerow([x, output_dict[int(x)][0], output_dict[int(x)][1]])
    elif comm_args.jalview_output:
        jalview_output(output_dict, comm_args)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))