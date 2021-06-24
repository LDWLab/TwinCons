import re, ntpath
from Bio import SeqIO
from Bio.PDB import DSSP
from Bio.PDB import PDBParser
#from Bio.PDB import ResidueDepth
'''Contains class for alignment groups'''

class AlignmentGroup:
    '''
    Class for a single group within an alignment.
    Must pass the alignment object and optionally a 
    structure object or secondary structure string
    '''
    DSSP_code_mycode = {'H':'H','B':'S','E':'S','G':'H','I':'H','T':'O','S':'O','-':'O'}
    def __init__(self,aln_obj, struc_path=None, sstruc_str=None):
        self.aln_obj = aln_obj
        
        self.struc_path = struc_path if struc_path is not None else None
        self.sstruc_str = sstruc_str if sstruc_str is not None else None

    def add_struc_path(self, struc_path):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.PDB import PDBParser, MMCIFParser
        from Bio.SeqUtils import seq1

        self.struc_path = struc_path
        if ntpath.splitext(self.struc_path)[1] == ".pdb":
            parser = PDBParser()
        elif ntpath.splitext(self.struc_path)[1] == ".cif":
            parser = MMCIFParser()
        else:
            raise IOError("Unrecognized structure file type! Please use .pdb or .cif files!")
        
        structure = parser.get_structure("none", self.struc_path)
        chains = list()
        for chain in structure.get_chains():
            chains.append(chain)
        if len(chains) != 1:
            raise IOError(f"When using structure files, they need to have a single chain!")
        sequence = str()
        seq_ix_mapping = dict()
        untrue_seq_ix = 1
        residues = list(chains[0].get_residues())
        for resi in residues:
            resi_id = resi.get_id()
            if not re.match(r' ', resi_id[2]):
                continue
            if re.match(r'^H_', resi_id[0]):
                continue
            sequence += resi.get_resname().replace(' ','')
            seq_ix_mapping[untrue_seq_ix] = int(resi.get_id()[1])
            untrue_seq_ix += 1

        if len(seq1(residues[seq_ix_mapping[1]].get_resname().replace(' ',''))) != 0:
            sequence = seq1(sequence)
        self.seq_ix_mapping = seq_ix_mapping
        self.struc_seq = SeqRecord(Seq(sequence))

    def create_aln_struc_mapping_with_mafft(self):
        from subprocess import Popen, PIPE
        from Bio import AlignIO
        from os import remove, path
        from warnings import warn

        aln_group_path = f"{path.dirname(path.realpath(__file__))}/TWCtempAln.txt"
        pdb_seq_path = f"{path.dirname(path.realpath(__file__))}/TWCtempStrucSeq.txt"
        mappingFileName = pdb_seq_path + ".map"
        tempfiles = [aln_group_path, pdb_seq_path, mappingFileName]
        for tempf in tempfiles:
            if path.isfile(tempf):
                warn(f"When using mafft to make structural mapping the working directory must be free of file {tempf}. Trying to delete the file.")
                remove(tempf)
                if path.isfile(tempf):
                    raise IOError(f"Couldn't delete the file {tempf} please remove it manually!")

        aln_group_fh = open(aln_group_path, "w")
        AlignIO.write(self.aln_obj, aln_group_fh, "fasta")
        aln_group_fh.close()

        pdb_seq_fh = open(pdb_seq_path, "w")
        SeqIO.write(self.struc_seq, pdb_seq_fh, "fasta")
        pdb_seq_fh.close()

        pipe = Popen(f"mafft --quiet --addfull {pdb_seq_path} --mapout {aln_group_path}; cat {mappingFileName}", stdout=PIPE, shell=True)
        output = pipe.communicate()[0]
        mapping_file = output.decode("ascii").split('\n#')[1]
        groupName = output.decode('ascii').split('>')[1].split('_')[0]
        firstLine = True
        mapping, bad_map_positions, fail_map = dict(), 0, False
        for line in mapping_file.split('\n'):
            if firstLine:
                firstLine = False
                continue
            row = line.split(', ')
            if len(row) < 3:
                continue
            if row[2] == '-':
                bad_map_positions += 1
                continue
            if row[1] == '-':
                fail_map = True
            mapping[int(row[2])] = self.seq_ix_mapping[int(row[1])]
        for tempf in tempfiles:
            remove(tempf)
        if fail_map:
            raise ValueError(f"Mapping between structure file {self.struc_path} and group {groupName} did not work properly!")
        if bad_map_positions > 0:
            warn(f"Mapping between structure file {self.struc_path} and group {groupName} is poor!\n Continue with caution!")
        self.mapping = mapping
        return mapping

    def _freq_iterator(self, column, aa_list, weight_aa_distr):
        '''Calculates gap adjusted frequency of each AA in the column.'''
        #Still doesn't handle ambiguous letters well
        if type(aa_list) == list:
            aa_list = ''.join(aa_list)
        if len(aa_list) >= 20:
            abs_length = 20
            adjsuted_column_list = ['-' if resi=='X' else resi for resi in column]
            all_residues = aa_list.replace('X', '')
        else:
            abs_length = 4
            adjsuted_column_list = ['-' if resi=='N' else resi for resi in column]
            aa_list.replace('N', '')

        M   =  len(adjsuted_column_list)
        
        #Gap adjustment
        num_gaps = adjsuted_column_list.count('-')
        if '-' in weight_aa_distr.keys():
            num_gaps = weight_aa_distr['-']*M
        gap_freq = num_gaps/abs_length
        frequency_list = list()
        
        # Number of residues in column
        for base in aa_list:
            n_i = adjsuted_column_list.count(base) # Number of residues of type i
            if base in weight_aa_distr.keys():     # In case of weighted
                n_i = weight_aa_distr[base]*M
            n_i += gap_freq
            P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
            frequency_list.append(P_i)
        return frequency_list

    def column_distribution_calculation(self, aa_list, alignment_length, seq_weights):
        '''Calculates AA distribution for the current alignment column'''
        column_distr = dict()
        col_ix = 0
        while col_ix < alignment_length:
            col_aalist = list()
            weighted_distr = dict()
            if len(seq_weights) > 0:
                col_aalist = list()
                row_ix = 0
                for col_aa in self.aln_obj[:, col_ix]:
                    if col_aa not in weighted_distr.keys():
                        weighted_distr[col_aa] = float()
                    weighted_distr[col_aa] += seq_weights[row_ix]
                    row_ix += 1
            col_aalist = self._freq_iterator(self.aln_obj[:, col_ix], aa_list, weighted_distr)
            col_ix += 1
            column_distr[col_ix] = col_aalist
        return column_distr

    def structure_loader(self,struc_to_aln_index_mapping):
        inv_map = {v: k for k, v in struc_to_aln_index_mapping.items()}
        parser = PDBParser()
        structure = parser.get_structure('current_structure',self.struc_path)
        model = structure[0]
        return inv_map, model

    def ss_map_creator(self,struc_to_aln_index_mapping):
        '''
        Connects the alignment mapping index and the secondary structural
        assignments from DSSP.
        '''
        ss_aln_index_map={}
        inv_map, model = self.structure_loader(struc_to_aln_index_mapping)
        dssp = DSSP(model, self.struc_path)
        for a_key in list(dssp.keys()):
            ss_aln_index_map[inv_map[a_key[1][1]]] = self.DSSP_code_mycode[dssp[a_key][2]]
        return ss_aln_index_map

    def depth_map_creator(self, struc_to_aln_index_mapping):
        '''Connects the alignment mapping index and the residue depth'''

        res_depth_aln_index_map={}
        inv_map, model = self.structure_loader(struc_to_aln_index_mapping)
        dssp = DSSP(model, self.struc_path)
        #rd = ResidueDepth(model)
        for a_key in list(dssp.keys()):
            if dssp[a_key][3] > 0.2:
                res_depth_aln_index_map[inv_map[a_key[1][1]]]='E'
            else:
                res_depth_aln_index_map[inv_map[a_key[1][1]]]='B'
        return res_depth_aln_index_map

    def both_map_creator(self, struc_to_aln_index_mapping):
        '''Connects the alignment mapping index and the residue depth'''
        sda={}
        inv_map, model = self.structure_loader(struc_to_aln_index_mapping)
        try:
            dssp = DSSP(model, self.struc_path)
        except OSError as e:
            raise OSError("DSSP failed with the following error:\n"+e)
        for a_key in list(dssp.keys()):
            if a_key[1][1] in inv_map.keys():
                if dssp[a_key][3] > 0.2:
                    sda[inv_map[a_key[1][1]]]='E'+self.DSSP_code_mycode[dssp[a_key][2]]
                else:
                    sda[inv_map[a_key[1][1]]]='B'+self.DSSP_code_mycode[dssp[a_key][2]]
        return sda

    def _return_alignment_obj(self):
        '''Returns current alignment object of this group'''
        return self.aln_obj