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

### YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('a') #Test when input is unexpected
    'X'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    if nucleotide == "G":
        return "C"
    if nucleotide == "C":
        return "G"
    else:
        return "X" 

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    This test is sufficient for testing this function. 
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
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    #Test a longer strand without stop codon that is not a multiple of 3
    >>> rest_of_ORG("ATGAGAGAGAT")
    'ATGAGAGAGAT'
    """
    index = 0
    i = 0
    while i < len(dna): #Can divide by 3?
        codon = dna[index:index+3]
        index += 3
        if codon in ['TAG', 'TAA', 'TGA']:
            #print 'Found a stop'
            return dna[0:index-3]
        else: 
            i += 1
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
    >>> find_all_ORFs_oneframe('AATGTAGATAGATGTGCCC')
    []
    >>> find_all_ORFs_oneframe('AGTATGTAGATAGAAATGTGCCC')
    'ATGTAGATAGAA', 'ATGTGCCC']
    The second test tests when the start codon is not present.
    Third test shows function works when ATG is not first codon

    """
    #Determine how many 'ATG' are in a multiple of three
    index = 0
    i = 0
    while i < len(dna):
        codon = dna[index:index+3]
        index += 3
        if codon == 'ATG':
            print 'Found a start'
            return dna[0:dna.find()]
        else: 
            i += 1
    return dna

    # codons = len(dna)/3 #How many codons are in the strand
    # ATG_index = []#Initialize lists
    # ORFS = []
    # for i in range(codons): 
    #     c_ind = i #codon index
    #     begin = c_ind*3
    #     end = (c_ind+1)*3
    #     if dna[begin:end] == 'ATG':
    #         #print begin, end
    #         ATG_index.append(begin)
    #     else:
    #         pass
    # for i in range(len(ATG_index)):
    #     if i == (len(ATG_index)-1): 
    #         ORFS.append(dna[ATG_index[i]:])
    #     else: 
    #         ORFS.append(dna[ATG_index[i]:ATG_index[i+1]])
    # return ORFS

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("GAGAGAGAGAG") #String that contains no ATG
    []
    """

    if 'ATG' in dna:
        print dna.count('ATG')
        index =[]
        for i in range(dna.count('ATG')):
            frame = dna[i:]
            index.append(frame.find('ATG')+i)
            print frame
        print index
    else:
        print 'no'

#find_all_ORFs("ATGCATGAATGTAG")

# def find_all_ORFs_both_strands(dna):
#     """ Finds all non-nested open reading frames in the given DNA sequence on both
#         strands.
        
#         dna: a DNA sequence
#         returns: a list of non-nested ORFs
#     >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
#     ['ATGCGAATG', 'ATGCTACATTCGCAT']
#     """
#     # TODO: implement this
#     pass


# def longest_ORF(dna):
#     """ Finds the longest ORF on both strands of the specified DNA and returns it
#         as a string
#     >>> longest_ORF("ATGCGAATGTAGCATCAAA")
#     'ATGCTACATTCGCAT'
#     """
#     # TODO: implement this
#     pass


# def longest_ORF_noncoding(dna, num_trials):
#     """ Computes the maximum length of the longest ORF over num_trials shuffles
#         of the specfied DNA sequence
        
#         dna: a DNA sequence
#         num_trials: the number of random shuffles
#         returns: the maximum length longest ORF """
#     # TODO: implement this
#     pass

# def coding_strand_to_AA(dna):
#     """ Computes the Protein encoded by a sequence of DNA.  This function
#         does not check for start and stop codons (it assumes that the input
#         DNA sequence represents an protein coding region).
        
#         dna: a DNA sequence represented as a string
#         returns: a string containing the sequence of amino acids encoded by the
#                  the input DNA fragment

#         >>> coding_strand_to_AA("ATGCGA")
#         'MR'
#         >>> coding_strand_to_AA("ATGCCCGCTTT")
#         'MPA'
#     """
#     # TODO: implement this
#     pass

# def gene_finder(dna, threshold):
#     """ Returns the amino acid sequences coded by all genes that have an ORF
#         larger than the specified threshold.
        
#         dna: a DNA sequence
#         threshold: the minimum length of the ORF for it to be considered a valid
#                    gene.
#         returns: a list of all amino acid sequences whose ORFs meet the minimum
#                  length specified.
#     """
#     # TODO: implement this
#     pass

# if __name__ == "__main__": #Program only runs if it is imported into the command line
#     import doctest
#     doctest.testmod()
