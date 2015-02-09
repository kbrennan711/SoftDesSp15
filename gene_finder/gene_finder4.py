# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Kelly Brennan

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
from load import load_seq


def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('a') #Test when input is unexpected
    ' '
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == "A":
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "G":
        return "C"
    elif nucleotide == "C":
        return "G"
    else:
        return " " 


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        
        These tests are sufficient for testing this function because there no boundary conditions

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    i = len(dna)-1
    rev_comp = ""
    while i >= 0:
        rev_comp += get_complement(dna[i])
        i -= 1 
    return rev_comp


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string

    >>> rest_of_ORF("ATGTGAA")
    'ATG'

    Test a longer strand without stop codon that is not a multiple of 3
    >>> rest_of_ORF("ATGAGAAGT") 
    'ATGAGAAGT'

    >>> rest_of_ORF("ATGAGAGAGAT")
    'ATGAGAGAGAT'
    """
    index = 0
    for index in range(0, len(dna), 3): #(start, end, increment)
        codon = dna[index:index+3]
        if codon in ['TAG', 'TAA', 'TGA']:
            return dna[0:index]
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    #Returns empty list when start codon is not present
    >>> find_all_ORFs_oneframe('AATGTAGATAGATGTGCCC')
    []

    #Returns correct ORFs when ATG is not first in DNA strand
    >>> find_all_ORFs_oneframe('AGTATGTAGATAGAAATGTGCCC')
    ['ATG', 'ATGTGCCC']
    """
    res = []
    index = 0
    while index < len(dna):
        codon = dna[index:index+3]
        if codon == 'ATG':
            res.append(rest_of_ORF(dna[index:]))
            index += len(res[-1]) #test when the DNA strand contains more than two ORFs
        else:
            index += 3
    return res


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    Test is sufficient because it tests all three reading frames
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']

    #Returns empty list when ATG is the not present in strand
    >>> find_all_ORFs("GAGAGAGAGAG")
    []

    #Test without stop codon. Returns DNA to the end of the strand
    >>> find_all_ORFs("ATGCATGAATG")
    ['ATGCATGAATG', 'ATGAATG', 'ATG']
    """
    ORFS = []
    if 'ATG' in dna:
        for i in range(3):
            ORFS += find_all_ORFs_oneframe(dna[i:])
    return ORFS


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCAT")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    Returns empty list if input is nothing
    >>> find_all_ORFs_both_strands("")
    []

    Test the three frames on the complement strand and one frame on reverse complement strand
    >>> find_all_ORFs_both_strands('ATGCATGAATGTAG')
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG', 'ATGCAT']

    Test the three frames of the reverse complement strand and one frame on complement strand
    >>> find_all_ORFs_both_strands('CTACATTCATGCAT')
    ['ATGCAT', 'ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    res = []
    res += find_all_ORFs(dna)
    res += find_all_ORFs(get_reverse_complement(dna))
    return res


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'

    Has to determine the longest DNA sequence out of four items in output of find_all_ORFs_both_strands
    >>> longest_ORF("TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCAGCATTGGTGTACATCTT")
    'ATGTACACCAATGCTGGCAGCCATACCTGGTTTAAGTCATTTATTGTTCAGAAACACAGTCA'

    Returns empty string when no ORFs are present in the DNA sequence
    >> longest_ORF("GAGAGAGA")

    """
    if len(find_all_ORFs_both_strands(dna)) != 0:
        return max(find_all_ORFs_both_strands(dna), key=len)
    return ""

def longest_ORF_noncoding(dna, num_trials): #Question: What is this really supposed to do?
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 

    >>> longest_ORF_noncoding('ATGCGAATGTAGCATCA', 200)
    17

    Longest ORF is 0 in non-coding region
    >>> longest_ORF_noncoding("GAGAGA", 1)
    0
    """
    shuffled_ORFs = []
    for i in range(num_trials):
        x = shuffle_string(dna)
        shuffled_ORFs.append(longest_ORF(x))
    return len(max(shuffled_ORFs, key = len))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'

        Returns nothing when there is no start codon.
    """
    start = dna[0:3]
    res = []
    if start == "ATG":
        res.append('M')
        for index in range(3, len(dna)-2, 3): # -2 to deal with nonmultiples of 3
            codon = dna[index:index+3]
            amino_acid = aa_table[codon]
            res.append(amino_acid)
    return ''.join(res)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        This function takes as input a sequence of DNA.  
        Use your longest_ORF_noncoding on the input DNA sequence to compute a conservative 
        threshold for distinguishing between genes and non-genes by running longest_ORF_noncoding 
        for 1500 trials.  Next, find all open reading frames on both strands, and then return a 
        list containing the amino acid sequence encoded by any open reading frames that are longer
        than the  threshold computed using longest_ORF_noncoding. 
    """
    threshold = longest_ORF_noncoding(dna, 3000) 
    print threshold
    all_ORFs = find_all_ORFs_both_strands(dna)
    filtered_ORFs = []
    primary_AA = []
    for i in all_ORFs:
        if len(i) >= threshold:
            filtered_ORFs.append(i)
    for i in filtered_ORFs:
        primary_AA.append(coding_strand_to_AA(i))
    return primary_AA


dna = load_seq("./data/X73525.fa")
print gene_finder(dna)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #doctest.run_docstring_examples(find_all_ORFs_both_strands, globals())
#>>>>>>> daf89177ec30bf2d3eb6f14c538e0a94cf857a77