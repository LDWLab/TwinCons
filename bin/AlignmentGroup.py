import re
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

	def create_struc_aln_mapping(self):
		'''
		Make a dictionary complementing the locations in 
		the alignment fasta file and the structure pdb file
		'''
		alignDict = {}
		struc_name = re.findall(r'(.*\/)(.*)(_.*.pdb)',self.struc_file)[0][1]
		for alignment in self.aln_obj: 							# Iterating through each sequence
			alignDict[alignment.id] = alignment.seq 			# Making a dictionary of the sequences
		dictseq={} 												# gives the complement between the location in alns (key) and the location in pdb (value)
		i=0
		a=1
		for alignment in self.aln_obj:
			if re.search(struc_name, alignment.id) is not None:
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
				dictList[v[0]]=k
		return dictList