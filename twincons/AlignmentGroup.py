import os, re, ntpath, shutil
import numpy as np
from Bio import SeqIO
from Bio.PDB import DSSP
from Bio.PDB import PDBParser
from twincons.twcSupportFunctions import alignment_array
'''Contains class for alignment groups'''

def locate_dssp_data():
    '''
    DSSP 4 reads its mmCIF dictionaries through libcifpp, whose conda build does not find them
    unless LIBCIFPP_DATA_DIR is set. When it is unset, point it at <prefix>/share/libcifpp next to
    the mkdssp executable, if the dictionaries are there.
    '''
    if 'LIBCIFPP_DATA_DIR' in os.environ:
        return
    executable = shutil.which('mkdssp') or shutil.which('dssp')
    if executable is None:
        return
    prefix = os.path.dirname(os.path.dirname(os.path.realpath(executable)))
    data_dir = os.path.join(prefix, 'share', 'libcifpp')
    if os.path.isfile(os.path.join(data_dir, 'mmcif_pdbx.dic')):
        os.environ['LIBCIFPP_DATA_DIR'] = data_dir

class AlignmentGroup:
    '''
    Class for a single group within an alignment.
    Must pass the alignment object and optionally a structure object
    and a sequence distribution (used for gap adjustment). When no
    sequence distribution is passed a uniform distribution is assumed.
    '''
    # DSSP 4 added P (polyproline II helix), grouped here with turns and coil.
    DSSP_code_mycode = {'H':'H','B':'S','E':'S','G':'H','I':'H','T':'O','S':'O','P':'O','-':'O'}
    def __init__(self, aln_obj, seq_distribution=None, struc_path=None):
        self.aln_obj = aln_obj
        self.uniq_resi_list = self._determineUniqResis(aln_obj)
        if seq_distribution is not None:
            if isinstance(seq_distribution, np.ndarray):
                self.seq_distribution = {self.uniq_resi_list[i] : seq_distribution[i] for i in range(len(seq_distribution))}
            elif isinstance(seq_distribution, dict):
                self.seq_distribution = seq_distribution
            else:
                raise IOError("Incorrect type of seq_distribution passed. Must be np.array or dict.")
        else:
            tempStorage = str()
            for entry in aln_obj:
                tempStorage += str(entry.seq).replace('-','').replace('\n','')
            self.seq_distribution = {i : tempStorage.count(i)/len(tempStorage) for i in set(tempStorage)}
        self.struc_path = struc_path

    def validateType(self, string, alphabet='protein'):
        '''Check that a string only contains values from an alphabet'''
        alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 
                     'rna': re.compile('^[acgun]*$', re.I),
                 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}

        if alphabets[alphabet].search(string) is not None:
             return True
        else:
             return False

    def _determineUniqResis(self, aln_obj):
        '''Used to determine the unique residues of the alignment.'''
        tempSeq = ''.join([str(x.seq).replace('-','') for x in aln_obj])
        if self.validateType(tempSeq):
            return ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        elif self.validateType(tempSeq, alphabet='rna'):
            return ['A','U','C','G']
        elif self.validateType(tempSeq, alphabet='dna'):
            return ['A','T','C','G']
        else:
            raise IOError("Wasn't able to determine the type of alignment. Do you have weird characters in the alignment?")

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
            raise IOError("When using structure files, they need to have a single chain!")
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
            if re.match(r'W', resi_id[0]):
                continue
            sequence += resi.get_resname().replace(' ','')
            seq_ix_mapping[untrue_seq_ix] = int(resi.get_id()[1])
            untrue_seq_ix += 1

        if len(seq1(residues[seq_ix_mapping[1]].get_resname().replace(' ',''))) != 0:
            sequence = seq1(sequence)
        self.seq_ix_mapping = seq_ix_mapping
        self.struc_seq = SeqRecord(Seq(sequence))

    def create_aln_struc_mapping_with_mafft(self):
        import os
        import subprocess
        import tempfile
        from Bio import AlignIO
        from warnings import warn
        from twincons.twcSupportFunctions import find_executable

        mafft = find_executable('mafft', 'to map structure residues onto the alignment')
        with tempfile.TemporaryDirectory(prefix='twincons_') as temp_dir:
            aln_group_path = os.path.join(temp_dir, 'aln_group.fas')
            pdb_seq_path = os.path.join(temp_dir, 'struc_seq.fas')
            with open(aln_group_path, "w") as aln_group_fh:
                AlignIO.write(self.aln_obj, aln_group_fh, "fasta")
            with open(pdb_seq_path, "w") as pdb_seq_fh:
                SeqIO.write(self.struc_seq, pdb_seq_fh, "fasta")

            result = subprocess.run([mafft, '--quiet', '--addfull', pdb_seq_path, '--mapout', aln_group_path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            map_path = pdb_seq_path + ".map"
            if result.returncode != 0 or not os.path.isfile(map_path):
                raise OSError(f"mafft --addfull failed with exit code {result.returncode}:\n{result.stderr}")
            with open(map_path) as map_handle:
                map_text = map_handle.read()
        mapping_file = map_text.split('\n#')[1]
        groupName = result.stdout.split('>')[1].split('_')[0]
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
        if fail_map:
            raise ValueError(f"Mapping between structure file {self.struc_path} and group {groupName} did not work properly!")
        if bad_map_positions > 0:
            warn(f"Mapping between structure file {self.struc_path} and group {groupName} is poor!\n Continue with caution!")
        self.mapping = mapping
        return mapping

    def column_distribution_calculation(self, aa_list, alignment_length, seq_weights):
        '''
        Returns {column index (1-based): gap adjusted frequency of each residue in aa_list}.
        Ambiguous residues (X for proteins, N for nucleotides) count as gaps. Gaps are
        redistributed according to seq_distribution, or uniformly for residues missing from it.
        With seq_weights, residue counts are replaced by the summed weights of the sequences
        carrying that residue, scaled by the number of sequences.
        '''
        aa_string = ''.join(aa_list)
        if len(aa_string) >= 20:
            abs_length, ambiguous = 20, b'X'
        else:
            abs_length, ambiguous = 4, b'N'
        seqs = alignment_array(self.aln_obj)[:, :alignment_length]
        M = seqs.shape[0]
        adjusted = np.where(seqs == ambiguous, b'-', seqs)
        num_gaps = (adjusted == b'-').sum(axis=0)
        weighted = len(seq_weights) > 0
        if weighted:
            weights = np.asarray(seq_weights, dtype=float)[:, None]
            def summed_weights(char):
                # Sequential sum over sequences keeps the floating point result of per-row accumulation.
                present = seqs == char
                return present.any(axis=0), (present * weights).cumsum(axis=0)[-1]
            gap_present, gap_weight = summed_weights(b'-')
            num_gaps = np.where(gap_present, gap_weight*M, num_gaps)
        frequencies = np.empty((seqs.shape[1], len(aa_string)))
        for i, base in enumerate(aa_string):
            n_i = (adjusted == base.encode('ascii')).sum(axis=0)
            if weighted:
                base_present, base_weight = summed_weights(base.encode('ascii'))
                n_i = np.where(base_present, base_weight*M, n_i)
            if base in self.seq_distribution:
                n_i = n_i + self.seq_distribution[base]*num_gaps
            else:
                n_i = n_i + num_gaps/abs_length
            frequencies[:, i] = n_i/float(M)
        return {col_ix: column.tolist() for col_ix, column in enumerate(frequencies, 1)}

    def structure_loader(self,struc_to_aln_index_mapping):
        locate_dssp_data()
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
            raise OSError(f"DSSP failed with the following error:\n{e}") from e
        for a_key in list(dssp.keys()):
            if a_key[1][1] in inv_map.keys():
                if dssp[a_key][3] > 0.2:
                    sda[inv_map[a_key[1][1]]]='E'+self.DSSP_code_mycode[dssp[a_key][2]]
                else:
                    sda[inv_map[a_key[1][1]]]='B'+self.DSSP_code_mycode[dssp[a_key][2]]
        return sda

    def getAAfrequenciesList (self):
        return [self.seq_distribution.get(aa, 0.0) for aa in self.uniq_resi_list]