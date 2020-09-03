#!/usr/bin/env python3
import unittest, re, pathlib, sys, os, itertools, filecmp
from Bio.SubsMat import MatrixInfo

sys.path.append(os.path.dirname(os.path.abspath(__name__)))
from bin import TwinCons

class TestTwinCons(unittest.TestCase):
    output_type_args = ['-csv', '-jv', '-p']
    entropy_type = ['-lg', '-e', '-rs', ['-mx', 'blosum62']]
    structure_mx = ['-ss', '-be', '-ssbe']
    nucleotide_mx = ['blastn', 'identity', 'trans']
    weigh_algorithms = ['', ['-w', 'pairwise']]

    def test_TWC_pseq_params(self):
        args_for_twc = ['-a', './tests/input_test_data/alns/uL02ab_txid_tagged.fas']
        output_files = list()
        argument_combinations = list(itertools.product(self.entropy_type, self.output_type_args, self.weigh_algorithms))
        arguments = [[tup for tup in x if tup] for x in argument_combinations if x]
        flat_args = list()
        for argset in arguments:
            tempargset = list()
            for arg in argset:
                if type(arg) == list:
                    tempargset.extend(arg)
                else:
                    tempargset.append(arg)
            flat_args.append(tempargset)
        
        for argset in flat_args:
            out_file_name = '_'.join(argset).replace('-', '')
            argset.extend(['-o', f'./tests/output_test_data/{out_file_name}'])
            output_files.append(f'./tests/output_test_data/{out_file_name}')
            args_for_twc.extend(argset)
            TwinCons.main(args_for_twc)
            args_for_twc = ['-a', './tests/input_test_data/alns/uL02ab_txid_tagged.fas']
    filecmp.cmp('file1.txt', 'file1.txt')

if __name__ == '__main__':
    unittest.main()