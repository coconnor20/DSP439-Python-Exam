import re
import pandas as pd
import sys

def parse_input(path):
    '''
    Input Parsing Function

    Args: 
    Path -> Filepath to a textfile containg genetic strings

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

    K -> size of string

    S -> Full String

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

    K -> size of string
    
    S -> Full String

    Returns:
    
    the minimum of the two possible k-mer calculations for a genomic string
    either cardinality of string - k + 1
    OR
    4^k
    '''
    return min((len(s)-k)+1, 4**k)

def calc_linguistic_comp(observed, possible):
    r'''
    Calculates Linguistic Complexity

    Args:
    observed(int) -> sum of observerd kmers \t
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
    

if __name__ == "__main__":
    main()