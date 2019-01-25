import re, random
from Bio.Seq import MutableSeq
from Bio.PDB import DSSP
from Bio.PDB import PDBParser
from Bio.PDB import ResidueDepth
from collections import Counter
"""
Contains class for alignment groups
"""

class AlginmentGroup:
	"""
	Class for a single group within an alignment.
	Must pass the alignment object and optionally a structure object or secondary structure string
	"""

	def __init__(self,aln_obj, struc_file=None, sstruc_str=None):
		self.aln_obj = aln_obj
		self.struc_file = struc_file if struc_file is not None else None
		self.sstruc_str = sstruc_str if sstruc_str is not None else None

	def create_aln_struc_mapping(self):
		'''
		Make a dictionary complementing the locations in 
		the alignment fasta file and the structure pdb file
		
		Check and think what to do when no mapping is possible 
		since we are passing the wrong struc (nonexistent in this group)
		'''
		alignDict = {}
		struc_name = re.findall(r'(.*\/)(.*)(_.*.pdb)',self.struc_file)[0][1]
		for alignment in self.aln_obj: 							# Iterating through each sequence
			alignDict[alignment.id] = alignment.seq 			# Making a dictionary of the sequences
		dictseq={} 												# gives the complement between the location in alns (key) and the location in pdb (value)
		i=0
		a=1
		for alignment in self.aln_obj:
			if re.search(struc_name, alignment.id) is not None: # else error as wrong structure was passed
				anchor_seq=alignment 							# anchor_seq will hold the sequence data from the name of the pdb file
				for x in anchor_seq:
					if (x=="-"):
						dictseq[a] = []
						dictseq[a].append(0)
						a+=1
					else:
						i+=1
						dictseq[a] = []
						dictseq[a].append(i)
						a+=1
		dictList={} 							# gives the complement between the location in pdb (key) and the location in alns (value) after removing gaps
		for k,v in dictseq.items():
			if v[0]==0:
				next
			else:
				dictList[k]=v[0]
		return dictList
	
	def _freq_iterator(self, column, aa_list):
		'''
		Calculates frequency of each AA in the column.
		'''
		col_aalist=[]
		for aa in aa_list:
			#print(aa, column.count(aa)/len(column.replace("-", "")))
			col_aalist.append(column.count(aa)/len(column))
		#print()
		return col_aalist

	def randomize_gaps(self, aa_list):
		'''
		Substitutes gaps in the alignment with a random choice from the present AAs.
		Should be an option to either use absolute random or random sampled from the distribution of the sequence.
		'''
		for aln in self.aln_obj:
			i = 0
			newaln=MutableSeq(str(aln.seq))
			aln_no_dash = str(aln.seq).replace('-', '')
			distribution = Counter(aln_no_dash)
			choice_distr = []
			for aa in aa_list:
				choice_distr.append(distribution[aa]/len(aln_no_dash))
			for resi in aln.seq:
				if resi == '-':
					#newaln[i] = choice(aa_list, p=choice_distr)
					newaln[i] = random.choice(aa_list)
				i+=1
			aln.seq=newaln
		return True
	
	def column_distribution_calculation(self, aa_list, alen):
		column_distr={}
		i=0
		while i < alen:
			col_aalist=self._freq_iterator(self.aln_obj[:, i], sorted(aa_list))
			i+=1
			column_distr[i] = col_aalist
		return column_distr

	def ss_map_creator(self):
		parser = PDBParser()
		structure = parser.get_structure('current_structure',self.struc_file)
		model = structure[0]
		dssp = DSSP(model, self.struc_file)
		a_key = list(dssp.keys())[2]
		print (dssp[a_key])

		rd = ResidueDepth(model)

	def _return_alignment_obj(self):
		"""
		Returns current alignment object of this group
		"""
		return self.aln_obj