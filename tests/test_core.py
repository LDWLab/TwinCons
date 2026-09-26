#!/usr/bin/env python3
import os
import tempfile
import unittest

import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from twincons import MatrixInfo, TwinCons
from twincons.AlignmentGroup import AlignmentGroup
from twincons.MatrixLoad import matrix_path
from twincons.twcSupportFunctions import slice_by_name

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ALIGNMENT_PATH = os.path.join(TEST_DIR, 'input_test_data', 'alns', 'uL02ab_txid_tagged.fas')


def make_alignment(records):
    return MultipleSeqAlignment([SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in records])


class TestSliceByName(unittest.TestCase):
    def test_group_name_that_prefixes_another_group(self):
        aln = make_alignment([('A_1', 'ACDE'), ('A_2', 'ACDF'), ('AB_1', 'GCDE'), ('AB_2', 'GCDF')])
        groups = slice_by_name(aln)
        self.assertEqual([rec.id for rec in groups['A']], ['A_1', 'A_2'])
        self.assertEqual([rec.id for rec in groups['AB']], ['AB_1', 'AB_2'])

    def test_group_name_with_regex_characters(self):
        aln = make_alignment([('G.1_a', 'ACDE'), ('GX1_b', 'ACDF')])
        groups = slice_by_name(aln)
        self.assertEqual([rec.id for rec in groups['G.1']], ['G.1_a'])

    def test_groups_keep_alignment_order(self):
        aln = make_alignment([('Zeta_1', 'ACDE'), ('Alpha_1', 'ACDF'), ('Zeta_2', 'GCDE')])
        self.assertEqual(list(slice_by_name(aln)), ['Zeta', 'Alpha'])


class TestSubstitutionMatrices(unittest.TestCase):
    def run_csv(self, *options):
        with tempfile.TemporaryDirectory() as output_dir:
            output_path = os.path.join(output_dir, 'out')
            TwinCons.main(['-a', ALIGNMENT_PATH, '-csv', '-o', output_path] + list(options))
            return os.path.getsize(output_path + '.csv')

    def test_pam_matrix_with_background_frequencies(self):
        self.assertGreater(self.run_csv('-mx', 'pam250'), 0)

    def test_matrix_without_background_frequencies_uses_uniform_baseline(self):
        self.assertGreater(self.run_csv('-mx', 'gonnet', '-bn', 'uniform'), 0)

    def test_matrix_without_background_frequencies_explains_fix(self):
        with self.assertRaisesRegex(IOError, '-bn uniform'):
            self.run_csv('-mx', 'gonnet')

    def test_argument_parsing_does_not_mutate_matrix_list(self):
        before = list(MatrixInfo.available_matrices)
        TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-csv'])
        TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-csv'])
        self.assertEqual(MatrixInfo.available_matrices, before)

    def test_packaged_matrices_exist(self):
        for parts in (['LG.dat'], ['BLOSUM', 'blosum62.out'], ['structureDerived', 'BEHOS.dat'], ['jp', 'blosum62.dat']):
            self.assertTrue(os.path.isfile(matrix_path(*parts)), parts)


class TestAlignmentGroup(unittest.TestCase):
    def test_seq_distribution_accepts_dict(self):
        aln = make_alignment([('A_1', 'ACDE'), ('A_2', 'ACDF')])
        distribution = {'A': 0.5, 'C': 0.5}
        self.assertEqual(AlignmentGroup(aln, seq_distribution=distribution).seq_distribution, distribution)

    def test_seq_distribution_accepts_array(self):
        aln = make_alignment([('A_1', 'ACDE'), ('A_2', 'ACDF')])
        group = AlignmentGroup(aln, seq_distribution=np.repeat(0.05, 20))
        self.assertEqual(group.seq_distribution['A'], 0.05)


if __name__ == '__main__':
    unittest.main()
