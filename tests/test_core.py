#!/usr/bin/env python3
import io
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from twincons import MatrixInfo, SequenceWeightFromTree, TwinCons
from twincons.AlignmentGroup import AlignmentGroup, locate_dssp_data
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


class TestOtherAlignments(unittest.TestCase):
    RNA_PATH = os.path.join(TEST_DIR, 'input_test_data', 'alns', 'AB_LSU_rRNA.fa')
    PROTEIN_PATH = os.path.join(TEST_DIR, 'input_test_data', 'alns', 'bS01-RNAP7Ca.fa')

    def scores(self, *args):
        output_dict = TwinCons.main(list(args) + ['-r'])[0]
        return [output_dict[position][0] for position in sorted(output_dict)]

    def test_nucleotide_matrix(self):
        scores = self.scores('-a', self.RNA_PATH, '-nc', '-mx', 'blastn')
        self.assertEqual(len(scores), read_align(self.RNA_PATH).get_alignment_length())
        self.assertTrue(np.isfinite(scores).all())

    def test_nucleotide_entropy_with_gap_removal(self):
        scores = self.scores('-a', self.RNA_PATH, '-nc', '-rs', '-cg')
        self.assertLess(len(scores), read_align(self.RNA_PATH).get_alignment_length())
        self.assertTrue(np.isfinite(scores).all())

    def test_groups_from_phylogenetic_tree(self):
        by_name = self.scores('-a', self.PROTEIN_PATH, '-lg')
        by_tree = self.scores('-a', self.PROTEIN_PATH, '-lg', '-phy')
        self.assertEqual(len(by_name), len(by_tree))
        self.assertTrue(np.isfinite(by_tree).all())


class TestClustalWWeights(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.group = slice_by_name(read_align(ALIGNMENT_PATH))['uL02b']

    def weights(self, aln):
        return np.array(SequenceWeightFromTree.calculate_weight_vector(aln, algorithm='clustalw'))

    def with_copies(self, index, copies):
        record = self.group[index]
        return MultipleSeqAlignment(list(self.group) + [SeqRecord(record.seq, id=f'{record.id}_copy{k}') for k in range(copies)])

    def test_hand_computed_tree(self):
        # A and B share the internal branch of length 2, so each gets 1 from it.
        tree = Phylo.read(io.StringIO('((A:1,B:1):2,C:3);'), 'newick')
        self.assertEqual(SequenceWeightFromTree.branch_sharing_weights(tree, ['A', 'B', 'C']), [2.0, 2.0, 3.0])

    def test_gapless_identity_distance(self):
        aln = make_alignment([('A_full1', 'ACDEFGHIKLMN'), ('A_full2', 'ACDEFGHIKLMQ'), ('A_frag1', 'ACDEF-------'),
                              ('A_frag2', 'WYVTS-------'), ('A_frag3', '-------IKLMN')])
        distances = SequenceWeightFromTree.gapless_identity_distance_matrix(aln)
        self.assertAlmostEqual(distances['A_full1', 'A_full2'], 1/12)
        self.assertEqual(distances['A_frag1', 'A_full1'], 0)
        self.assertEqual(distances['A_frag1', 'A_frag2'], 1)
        self.assertEqual(distances['A_frag1', 'A_frag3'], 1)

    def test_shared_gaps_do_not_make_fragments_similar(self):
        aln = make_alignment([('A_full1', 'ACDEFGHIKLMN'), ('A_full2', 'ACDEFGHIKLMQ'), ('A_full3', 'ACDEYGHIKRMN'),
                              ('A_frag1', 'ACDEF-------'), ('A_frag2', 'WYVTS-------')])
        weights = self.weights(aln)
        # frag1 repeats full1's residues, frag2's residues are unique; their shared gaps are irrelevant.
        self.assertLessEqual(weights[3], weights[0])
        self.assertTrue((weights[4] > weights[:4]).all())

    def test_matches_path_sum_definition(self):
        tree = SequenceWeightFromTree.tree_from_distances(SequenceWeightFromTree.gapless_identity_distance_matrix(self.group))
        leaves_below = {clade: len(clade.get_terminals()) for clade in tree.find_clades()}
        expected = [sum((clade.branch_length or 0) / leaves_below[clade] for clade in tree.get_path(record.id))
                    for record in self.group]
        result = SequenceWeightFromTree.branch_sharing_weights(tree, [record.id for record in self.group])
        np.testing.assert_allclose(result, expected, rtol=1e-12)

    def test_weights_are_positive_and_normalized(self):
        weights = self.weights(self.group)
        self.assertEqual(len(weights), len(self.group))
        self.assertTrue((weights > 0).all())
        self.assertAlmostEqual(weights.sum(), 1.0)

    def test_identical_sequences_share_weight_equally(self):
        weights = self.weights(self.with_copies(30, 3))
        copies = weights[[30, len(self.group), len(self.group) + 1, len(self.group) + 2]]
        np.testing.assert_allclose(copies, copies[0], rtol=1e-12)

    def test_duplicates_do_not_accumulate_weight(self):
        original = self.weights(self.group)[30]
        weights = self.weights(self.with_copies(30, 3))
        together = weights[[30, len(self.group), len(self.group) + 1, len(self.group) + 2]].sum()
        # Pairwise weights give four copies about 3.7x the original weight on this group.
        self.assertLess(together, 1.5 * original)

    def test_divergent_sequence_outweighs_redundant_clade(self):
        aln = make_alignment([('A_1', 'ACDEFGHIKL'), ('A_2', 'ACDEFGHIKM'), ('A_3', 'ACDEFGHIKN'),
                              ('A_4', 'ACDEFGHIKP'), ('A_5', 'WYVTSRQPNM')])
        weights = self.weights(aln)
        self.assertTrue((weights[4] > weights[:4]).all())

    def test_identical_sequences_only(self):
        aln = make_alignment([('A_1', 'ACDE'), ('A_2', 'ACDE'), ('A_3', 'ACDE')])
        for algorithm in ('pairwise', 'clustalw'):
            self.assertEqual(SequenceWeightFromTree.calculate_weight_vector(aln, algorithm=algorithm), [1/3] * 3)

    def test_single_sequence(self):
        aln = make_alignment([('A_1', 'ACDE')])
        self.assertEqual(SequenceWeightFromTree.calculate_weight_vector(aln, algorithm='clustalw'), [1.0])

    def test_scoring_with_clustalw_weights(self):
        def scores(*options):
            output_dict = TwinCons.main(['-a', ALIGNMENT_PATH, '-r'] + list(options))[0]
            return np.array([output_dict[position][0] for position in sorted(output_dict)])
        for matrix in (['-lg'], ['-rs']):
            weighted = scores(*matrix, '-w', 'clustalw')
            self.assertTrue(np.isfinite(weighted).all())
            self.assertFalse(np.allclose(weighted, scores(*matrix)))


class TestDsspDataLocation(unittest.TestCase):
    def fake_prefix(self, with_dictionary=True):
        prefix = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, prefix)
        os.makedirs(os.path.join(prefix, 'bin'))
        executable = os.path.join(prefix, 'bin', 'mkdssp')
        with open(executable, 'w') as fh:
            fh.write('')
        os.chmod(executable, 0o755)
        os.makedirs(os.path.join(prefix, 'share', 'libcifpp'))
        if with_dictionary:
            with open(os.path.join(prefix, 'share', 'libcifpp', 'mmcif_pdbx.dic'), 'w') as fh:
                fh.write('')
        return prefix

    def environment(self, prefix, **extra):
        return mock.patch.dict(os.environ, {'PATH': os.path.join(prefix, 'bin'), **extra}, clear=True)

    def test_points_dssp_at_its_dictionaries(self):
        prefix = self.fake_prefix()
        with self.environment(prefix):
            locate_dssp_data()
            self.assertEqual(os.path.realpath(os.environ['LIBCIFPP_DATA_DIR']),
                             os.path.realpath(os.path.join(prefix, 'share', 'libcifpp')))

    def test_keeps_an_existing_setting(self):
        prefix = self.fake_prefix()
        with self.environment(prefix, LIBCIFPP_DATA_DIR='/somewhere/else'):
            locate_dssp_data()
            self.assertEqual(os.environ['LIBCIFPP_DATA_DIR'], '/somewhere/else')

    def test_leaves_unset_without_dictionaries(self):
        prefix = self.fake_prefix(with_dictionary=False)
        with self.environment(prefix):
            locate_dssp_data()
            self.assertNotIn('LIBCIFPP_DATA_DIR', os.environ)


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
