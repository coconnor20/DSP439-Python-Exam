import unittest
from exam_four import *

class TestKmers(unittest.TestCase):
    
    def test_possiblekmers(self):
        test_str = 'ATTTGGATT'
        pre_comp_kmers = [4,8,7,6,5,4,3,2,1]
        post_comp_kmers = []
        for i in range(1, len(test_str)+1):
            post_comp_kmers.append(get_possible_kmers(i, test_str))

        self.assertEqual(pre_comp_kmers, post_comp_kmers)

    def test_observedkmers(self):
        test_str = 'ATTTGGATT'
        pre_comp_kmers = [3,5,6,6,5,4,3,2,1]
        post_comp_kmers = []
        for i in range(1, len(test_str)+1):
            post_comp_kmers.append(len(get_observed_kmers(i, test_str)))

        self.assertEqual(pre_comp_kmers, post_comp_kmers)

    def test_linguisticcomplexity(self):
        test_str = 'ATTTGGATT'
        pre_comp_lc = 0.875
        post_comp_kmers_o = []
        post_comp_kmers_p = []
        for i in range(1, len(test_str)+1):
            post_comp_kmers_p.append(get_possible_kmers(i, test_str))
            post_comp_kmers_o.append(len(get_observed_kmers(i, test_str)))
        post_calc_lc = calc_linguistic_comp(sum(post_comp_kmers_o), sum(post_comp_kmers_p))
        self.assertEqual(pre_comp_lc, post_calc_lc)

if __name__ == '__main__':
    unittest.main()