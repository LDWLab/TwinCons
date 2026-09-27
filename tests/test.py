#!/usr/bin/env python3
import itertools
import os
import tempfile
import unittest

from twincons import TwinCons

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ALIGNMENT_PATH = os.path.join(TEST_DIR, 'input_test_data', 'alns', 'uL02ab_txid_tagged.fas')
OUTPUT_EXTENSION = {'-csv': '.csv', '-jv': '.jlv', '-p': '.svg'}


class TestTwinCons(unittest.TestCase):
    output_type_args = ['-csv', '-jv', '-p']
    entropy_type = ['-lg', '-e', '-rs', ['-mx', 'blosum62']]
    structure_mx = ['-ss', '-be', '-ssbe']
    nucleotide_mx = ['blastn', 'identity', 'trans']
    weigh_algorithms = ['', ['-w', 'pairwise']]

    def test_TWC_pseq_params(self):
        with tempfile.TemporaryDirectory() as output_dir:
            combinations = itertools.product(self.entropy_type, self.output_type_args, self.weigh_algorithms)
            for entropy_arg, output_arg, weigh_arg in combinations:
                argset = list()
                for arg in (entropy_arg, output_arg, weigh_arg):
                    if isinstance(arg, list):
                        argset.extend(arg)
                    elif arg:
                        argset.append(arg)
                out_file_name = '_'.join(argset).replace('-', '')
                output_path = os.path.join(output_dir, out_file_name)
                args_for_twc = ['-a', ALIGNMENT_PATH, '-o', output_path] + argset
                with self.subTest(args=argset):
                    TwinCons.main(args_for_twc)
                    output_file = output_path + OUTPUT_EXTENSION[output_arg]
                    self.assertGreater(os.path.getsize(output_file), 0)


if __name__ == '__main__':
    unittest.main()
