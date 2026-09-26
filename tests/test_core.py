#!/usr/bin/env python3
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from twincons import MatrixInfo, SequenceWeightFromTree, TwinCons
from twincons.AlignmentGroup import AlignmentGroup
from twincons.CompositionalAdjustment import CompositionalAdjustmentError, adjust_matrix
from twincons.MatrixLoad import matrix_path
from twincons.twcSupportFunctions import read_align, slice_by_name

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


@unittest.skipUnless(shutil.which('mafft'), 'mafft is not installed')
class TestMergedAlignment(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.group_paths = list()
        for name, group in slice_by_name(read_align(ALIGNMENT_PATH)).items():
            path = os.path.join(self.temp_dir.name, f'{name}.fas')
            AlignIO.write(group, path, 'fasta')
            self.group_paths.append(path)

    def path(self, name):
        return os.path.join(self.temp_dir.name, name)

    def test_merge_only(self):
        self.assertIsNone(TwinCons.main(['-a'] + self.group_paths + ['-ma', self.path('merged.fas')]))
        merged = read_align(self.path('merged.fas'))
        self.assertEqual(len(merged), 122)
        self.assertEqual(list(slice_by_name(merged)), ['1', '2'])

    def test_merged_alignment_can_be_rescored(self):
        TwinCons.main(['-a'] + self.group_paths + ['-ma', self.path('merged.fas'), '-lg', '-csv', '-o', self.path('merged')])
        TwinCons.main(['-a', self.path('merged.fas'), '-lg', '-csv', '-o', self.path('rescored')])
        with open(self.path('merged.csv')) as merged, open(self.path('rescored.csv')) as rescored:
            self.assertEqual(merged.read(), rescored.read())

    def test_no_files_left_without_merged_alignment_option(self):
        previous = os.getcwd()
        os.chdir(self.temp_dir.name)
        self.addCleanup(os.chdir, previous)
        TwinCons.main(['-a'] + self.group_paths + ['-lg', '-csv', '-o', self.path('scores')])
        self.assertEqual(sorted(os.listdir(self.temp_dir.name)), sorted([os.path.basename(p) for p in self.group_paths] + ['scores.csv']))


class TestArgumentValidation(unittest.TestCase):
    def test_merged_alignment_needs_two_files(self):
        with self.assertRaises(SystemExit), mock.patch('sys.stderr'):
            TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-ma', 'merged.fas'])

    def test_voronoi_samples(self):
        args = TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-lg', '-csv'])
        self.assertEqual(args.voronoi_samples, SequenceWeightFromTree.DEFAULT_VORONOI_SAMPLES)
        args = TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-lg', '-csv', '-w', 'voronoi', '-vs', '500'])
        self.assertEqual(args.voronoi_samples, 500)
        with self.assertRaises(SystemExit), mock.patch('sys.stderr'):
            TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-lg', '-csv', '-vs', '0'])

    def test_output_option_required(self):
        with self.assertRaises(SystemExit), mock.patch('sys.stderr'):
            TwinCons.create_and_parse_argument_options(['-a', ALIGNMENT_PATH, '-lg'])


class TestExternalPrograms(unittest.TestCase):
    def test_missing_mafft_is_reported(self):
        with mock.patch('shutil.which', return_value=None):
            with self.assertRaisesRegex(OSError, 'mafft was not found on PATH'):
                TwinCons.run_mafft([ALIGNMENT_PATH, ALIGNMENT_PATH])



class TestCompositionalAdjustment(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.groups = list(slice_by_name(read_align(ALIGNMENT_PATH)).values())
        cls.freqs = [AlignmentGroup(group).getAAfrequenciesList() for group in cls.groups]
        cls.length = cls.groups[0].get_alignment_length()
        cls.blosum62 = np.loadtxt(matrix_path('jp', 'blosum62.dat'))

    def test_matches_original_newton_direct_solve(self):
        # Output of the formerly bundled newton_direct_solve binary for the same inputs (9 decimals).
        with open(os.path.join(TEST_DIR, 'input_test_data', 'compositional_adjustment', 'newton_direct_solve_blosum62_uL02ab.txt')) as fh:
            expected = np.array([[float(value) for value in line.split()] for line in fh if line.strip()])
        result = adjust_matrix(self.blosum62, self.freqs[0], self.freqs[1], self.length, self.length)
        np.testing.assert_allclose(result, expected, atol=1e-8, rtol=0)

    def test_group_missing_an_amino_acid(self):
        aln = make_alignment([('A_1', 'ACDEFGHIK'), ('A_2', 'ACDEFGHIL'), ('B_1', 'MNPQRSTVY'), ('B_2', 'MNPQRSTVW')])
        frequencies = AlignmentGroup(slice_by_name(aln)['A']).getAAfrequenciesList()
        self.assertEqual(frequencies[AlignmentGroup(aln).uniq_resi_list.index('W')], 0.0)
        self.assertEqual(adjust_matrix(self.blosum62, frequencies, self.freqs[1], 9, 9).shape, (20, 20))

    def test_non_convergence_is_reported(self):
        uniform = [0.05] * 20
        with self.assertRaises(CompositionalAdjustmentError):
            adjust_matrix(np.loadtxt(matrix_path('jp', 'blosum35.dat')), uniform, uniform, 500, 500)

    def test_command_line(self):
        with tempfile.TemporaryDirectory() as output_dir:
            TwinCons.main(['-a', ALIGNMENT_PATH, '-mx', 'blosum62', '-ca', '-csv', '-o', os.path.join(output_dir, 'out')])
            self.assertGreater(os.path.getsize(os.path.join(output_dir, 'out.csv')), 0)


class TestVectorizedAlgorithms(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aln = read_align(ALIGNMENT_PATH)[:25]

    def test_distance_matrix_matches_biopython(self):
        for model in ('identity', 'blosum62'):
            calculator = DistanceCalculator(model)
            expected = calculator.get_distance(self.aln)
            result = SequenceWeightFromTree.distance_matrix(self.aln, calculator)
            self.assertEqual(result.names, expected.names)
            self.assertEqual(result.matrix, expected.matrix)

    def test_leaf_distance_sums_match_tree_distances(self):
        tree = SequenceWeightFromTree.tree_construct(self.aln)
        names = [record.id for record in self.aln]
        expected = [sum(tree.distance(a, b) for b in names) for a in names]
        np.testing.assert_allclose(SequenceWeightFromTree.leaf_distance_sums(tree, names), expected, rtol=1e-12)

    def test_voronoi_weights_are_reproducible_distribution(self):
        np.random.seed(0)
        first = SequenceWeightFromTree.calculate_weight_vector(self.aln, algorithm='voronoi', repeat=200)
        np.random.seed(0)
        second = SequenceWeightFromTree.calculate_weight_vector(self.aln, algorithm='voronoi', repeat=200)
        self.assertEqual(first, second)
        self.assertAlmostEqual(sum(first), 1.0)

    def test_remove_extremely_gapped_regions(self):
        aln = make_alignment([('A_1', 'A--CD-'), ('A_2', 'A-EC--'), ('B_1', 'A--CDE'), ('B_2', 'AG-C--')])
        mapping, trimmed, length = TwinCons.remove_extremely_gapped_regions(aln, 0.5, {})
        self.assertEqual([str(record.seq) for record in trimmed], ['ACD', 'AC-', 'ACD', 'AC-'])
        self.assertEqual([record.id for record in trimmed], ['A_1', 'A_2', 'B_1', 'B_2'])
        self.assertEqual(length, 3)
        self.assertEqual(mapping, {1: 1, 2: 4, 3: 5})


class TestAlignmentGroup(unittest.TestCase):
    def test_seq_distribution_accepts_dict(self):
        aln = make_alignment([('A_1', 'ACDE'), ('A_2', 'ACDF')])
        distribution = {'A': 0.5, 'C': 0.5}
        self.assertEqual(AlignmentGroup(aln, seq_distribution=distribution).seq_distribution, distribution)

    def test_all_dssp4_secondary_structure_codes_are_mapped(self):
        self.assertLessEqual(set('HBEGIPTS-'), set(AlignmentGroup.DSSP_code_mycode))

    def test_seq_distribution_accepts_array(self):
        aln = make_alignment([('A_1', 'ACDE'), ('A_2', 'ACDF')])
        group = AlignmentGroup(aln, seq_distribution=np.repeat(0.05, 20))
        self.assertEqual(group.seq_distribution['A'], 0.05)


if __name__ == '__main__':
    unittest.main()
