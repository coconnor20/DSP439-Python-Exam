import unittest
import pandas as pd
import sys

def parse_input(path):
    '''
    Input Parsing Function

    Args: 
    Path(str) -> Filepath to a textfile containg genetic strings

    Returns:
    List of Genetic strings
    '''
    strs = []
    with open(path, "r") as input:
        chars = ["A", "G", "T", "C"]
        lines = input.readlines()
        stripped_strs = [i.strip('\n') for i in lines]
        stripped_strs = [i.strip(';') for i in stripped_strs]
        for i in stripped_strs:
            if i[len(i)-1] not in chars:
                i = i[-1]
    return stripped_strs

def get_observed_kmers(k, s):
    '''
    Get observerd Kmers Function

    Args: 

    K(int) -> size of string

    S(str) -> Full String

    Returns:

    A list of all the observed K-Mers


    '''
    seen = {}
    for i in range(0, len(s)):
        # k-mers are just a slice of the list from i to i+k
        s_str = s[i:i+k]
        if s_str in seen:
            continue
        else:
            #get rid of non valid k-mers at the end of the list
            if len(s_str) == k:
                seen.update({s_str: True})
            else:
                continue
    return list(seen.keys())


def get_possible_kmers(k, s):
    '''
    Get Possible Kmers Function

    Args: 

    K(int) -> size of string
    
    S(str) -> Full String

    Returns:
    
    the minimum of the two possible k-mer calculations for a genomic string
    either cardinality of string - k + 1
    OR
    4^k
    '''
    return min((len(s)-k)+1, 4**k)

def calc_linguistic_comp(observed, possible):
    '''
    Calculates Linguistic Complexity

    Args:

    observed(int) -> sum of observerd kmers 
    
    possible(int) -> sum of possible kmers

    Returns:

    float represntaion of the linguistic complexity
    '''
    return observed/possible

def make_dataframe(possible, observed, k):
    df = pd.DataFrame(zip(k, observed, possible), columns=["K", "Observed", "Possible"])
    return df

def main():
    path = sys.argv[1]
    input_strs = parse_input(path)
    possible_kmers = []
    observed_kmers = []
    k_list = []
    dfs = []
    print(input_strs)
    for i in input_strs:
        for j in range(1, len(i)+1):
            possible_kmers.append(get_possible_kmers(j, i))
            observed_kmers.append(len(get_observed_kmers(j, i)))
            k_list.append(j)
        dfs.append(make_dataframe(possible_kmers.copy(), observed_kmers.copy(), k_list.copy()))
        print(calc_linguistic_comp(sum(observed_kmers), sum(possible_kmers)))
        possible_kmers.clear()
        observed_kmers.clear()
        k_list.clear()
    
    for i in range(0, len(dfs)):
        fname = "kmers_" + str(i) + ".csv"
        f = open(fname, "w")
        f.write(dfs[i].to_csv(index=False))
    

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


if __name__ == "__main__":
    main()