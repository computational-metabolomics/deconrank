from __future__ import (
    print_function
)
import unittest
import os
from deconrank import Deconrank
from deconrank.scoring import load_score_table
import example_data_s

def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)

class TestDeconrank(unittest.TestCase):

    def test_grouping(self):
        '''Test the the grouping stage of the deconrank workflow'''
        dr = Deconrank(to_test_data('peaklist_positive.csv'))
        dr.group()
        self.assertEqual(dr.adduct_groups, example_data_s.adduct_groups)

        self.assertEqual(dr.isotope_groups, example_data_s.isotope_groups)

        self.assertEqual(dr.header, example_data_s.header)

        self.assertEqual(dr.grouped_features, example_data_s.grouped_features)

        self.assertEqual(dr.grouped_features_nm, example_data_s.grouped_features_nm)


    def test_scoring(self):
        '''Test the the scoring stage of the deconrank workflow'''
        dr = Deconrank(to_test_data('peaklist_positive.csv'))
        self.assertEqual(dr.scores, example_data_s.score_rules)
        self.assertEqual(dr.polarity, 'POS')

        dr.group()
        dr.score()

        l = dr.scored_adduct_list

        # check each row has the same number of elements
        self.assertTrue(all(len(i) == len(l[0]) for i in l))

        # Check that the second tier is scored as 11
        st = [i for i in l if i[4] == 'NA']
        self.assertTrue(all(i[2] == 11 for i in st))
        cnames = ('peakID', 'mz', 'intensity', 'isotopes', 'adduct', 'groupid', 'bestAdductScore', 'clustn',
                  'medianPurity','intp', 'addp', 'clustnp', 'pp', 'intensityw', 'adductw', 'clustnw', 'purityw',
                  'totalS', 'excluded', 'excludedFinal')


        self.assertEqual(len(dr.d_table), 140)
        self.assertEqual(dr.d_table.dtype.names, cnames)

        self.assertEqual(list(dr.d_table['bestAdductScore']), example_data_s.bas)
        self.assertEqual(list(dr.d_table['excluded']), [0] * 140)

        l = [round(i, 3) for i in dr.d_table['totalS']]

        self.assertEqual(list(dr.d_table['bestAdductScore']), example_data_s.bas)


    def test_filter(self):
        '''Test the the filter stage of the deconrank workflow'''
        dr = Deconrank(to_test_data('peaklist_positive.csv'))

        dr.group()
        dr.score()
        dr.filter()

        self.assertEqual(list(dr.d_table['excluded']), example_data_s.excluded)

    def test_dims_targets(self):
        '''Test the the filter stage of the deconrank workflow'''
        dr = Deconrank(to_test_data('A01_Polar_Daph_WAX1_Phenyl_LCMS_Neg_DIMS_annotated.csv'),
                       out_dir=to_test_data('results'))
        dr.group()
        dr.score()
        dr.filter()
        dr.dims_targets()

        self.assertEqual(list(dr.d_table['excludedFinal']), example_data_s.excluded_final_dims)

    def test_cli_dims(self):
        '''Test command line for DIMS'''
        import subprocess
        import csv

        in_pth = os.path.abspath(to_test_data('A01_Polar_Daph_WAX1_Phenyl_LCMS_Neg_DIMS_annotated.csv'))
        out_dir = os.path.abspath(to_test_data('results_cli_dims'))

        spth = os.path.join(os.path.dirname((os.path.dirname(os.path.realpath(__file__)))), 'deconrank', 'deconrank.py')

        cmd = 'python {spth} -i {in_pth} -o {out_pth} '.format(in_pth=in_pth, out_pth=out_dir, spth=spth)

        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output, err = popen.communicate()
        print(output)
        print(err)

        with open(os.path.join(out_dir, 'A01_Polar_Daph_WAX1_Phenyl_LCMS_Neg_DIMS_annotated_scores.csv'), 'rb') as csvfile:
            r = csv.reader(csvfile, delimiter=',')
            next(r, None)
            l = [tuple(i) for i in r]

        d_table = load_score_table(l)

        self.assertEqual(list(d_table['excludedFinal']), example_data_s.excluded_final_dims)












if __name__ == '__main__':
    unittest.main()