#!/usr/bin/env python3
import unittest, re, pathlib, sys, os
from Bio.SubsMat import MatrixInfo

sys.path.append(os.path.dirname(os.path.abspath(__name__)))
from bin import TwinCons

class TestTwinCons(unittest.TestCase):
    output_type_args = [['-pml', 'unix'], '-csv', '-jv', '-p']
    entropy_type = ['-lg', '-e', '-rs', ['-mx', 'blosum62']]
    structure_mx = ['-ss', '-be', '-ssbe']
    nucleotide_mx = ['blastn', 'identity', 'trans']
    weigh_algorithms = ['', ['-w', 'pairwise'], ['-w', 'voronoi']]

    def test_TWC_pseq_params(self):
        args_for_twc = ['-a', './tests/input_test_data/alns/uL02ab_txid_tagged.fas']
        args_file_name = []
        output_files = list()
        for w_alg in self.weigh_algorithms:
            w_arg = ['-w', w_alg]
            if w_alg != '':
                args_file_name.extend(w_arg)
            for out_arg in self.output_type_args:
                args_file_name.extend(out_arg)
                for entr_arg in self.entropy_type:
                    if re.match('_', entr_arg):
                        args_file_name.extend(entr_arg.split('_'))
                    else:
                        args_file_name.append(entr_arg)
                    out_file_name = '_'.join(args_file_name).replace('-', '')
                    args_file_name.extend(['-o', f'./tests/output_test_data/{out_file_name}'])
                    output_files.append(f'./tests/output_test_data/{out_file_name}')
                    args_for_twc.extend(args_file_name)
                    TwinCons.main(args_for_twc)
                    args_for_twc = ['-a', './tests/input_test_data/alns/uL02ab_txid_tagged.fas']


if __name__ == '__main__':
    unittest.main()