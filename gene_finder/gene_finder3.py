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
    ' '
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
        return " " 
#swap = {"A":"T", "T":"A", "C":"G", "G":"c"}
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

    #Test a longer strand without stop codon that is not a multiple of 3
    >>>rest_of_ORF("ATGAGAAGT") 
    'ATGAGAAGT'

    >>> rest_of_ORF("ATGAGAGAGAT")
    'ATGAGAGAGAT'
    """
    index = 0
    i = 0
    while i < len(dna)/3: #Can divide by 3?
        codon = dna[index:index+3]
        index += 3
        if codon in ['TAG', 'TAA', 'TGA']:
            return dna[0:index-3]
        else: 
            i += 1
    return dna
#print rest_of_ORF("ATGAGAAGT") 

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
    #Determine how many 'ATG's are in a multiple of three
    index = 0
    i = 0
    start = []
    res = []
    while i < len(dna):
        codon = dna[index:index+3]
        index += 3
        if codon == 'ATG':
            res.append(rest_of_ORF(dna[index-3:]))
            index += len(res[-1]) #test when the DNA strand contains more than two ORFs
        else: 
            i += 1
    return res

# x = shuffle_string("ATGCATGAATGTAGATAGATGTGCCC")
# print x
# print find_all_ORFs_oneframe(x)
#print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
#print find_all_ORFs_oneframe("ATGCGAATGTAGCATCAAA")

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    #Test is sufficient because it tests all three reading frames
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
        f = 0
        while f < 3:
            ORFS += find_all_ORFs_oneframe(dna[f:])
            f += 1
        return ORFS
    else:
        return ORFS
# x = shuffle_string("ATGCATGAATGTAGATAGATGTGCCC")
# print x
# print find_all_ORFs(x)
#print "Find ORFS:", 
#print find_all_ORFs("ATGCGAATGTAGCATCAAA")

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

    Test all six frames of the DNA
    >>> find_all_ORFs_both_strands('TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCGATGGTTATCTT')
    ['ATGCATGAATGGTA', 'ATGAATGGTA', 'ATGGTA', 'ATGCAT']
    """
    # TODO: implement this
    res = []
    if not len(find_all_ORFs(dna)) == 0: 
        res += find_all_ORFs(dna)
    if not len(find_all_ORFs(get_reverse_complement(dna))) == 0:
        #print 'Entered loop'
        #print get_reverse_complement(dna)
        res += find_all_ORFs(get_reverse_complement(dna))
    return res
#print get_reverse_complement("TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCAGCATTGGTGTACATCTT")
#print find_all_ORFs_both_strands("TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCAGCATTGGTGTACATCTT")
    #Long sequence: "GGGATCGATGCCCCTTAAAGAGTTTACATATTGCTGGAGGCGTTAACCCCGGAT")

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        #TEST - what happens when the length of both strands are equal?
        #Cannot be given an empty string

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'

    Has to determine the longest DNA sequence out of four items in the list.
    >>> longest_ORF("TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCAGCATTGGTGTACATCTT")
    ATGTACACCAATGCTGGCAGCCATACCTGGTTTAAGTCATTTATTGTTCAGAAACACAGTCA

    Returns empty string when no ORFs are present in the DNA sequence
    >> longest_ORF("GAGAGAGA")

    """
    if len(find_all_ORFs_both_strands(dna)) == 0:
        return ""
    else: 
        return max(find_all_ORFs_both_strands(dna), key = len)

#print longest_ORF("GAGAGAGA")
#print longest_ORF("TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCAGCATTGGTGTACATCTT")

def longest_ORF_noncoding(dna, num_trials): #Question: What is this really supposed to do?
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 

    >>> longest_ORF_noncoding('ATGCGAATGTAGCATCA', 10)
    17

    Longest ORF is 0 in non-coding region
    >>> longest_ORF_noncoding("GAGAGA", 1)
    0
    """
    # TODO: implement this
    shuffled_ORFs = []
    i = 0
    while i < num_trials: 
        x = shuffle_string(dna)
        #print i
        #print 'Shuffled DNA', x
        #print 'Longest ORF', longest_ORF(x)
        shuffled_ORFs.append(longest_ORF(x))
        i += 1
    #print "Shuffled DNA list",shuffled_ORFs
    #print max(shuffled_ORFs, key = len)
    return len(max(shuffled_ORFs, key = len))
    #print longest_ORF(shuffled_ORFs)
    #return len(longest_ORF(shuffled_ORFs))

#print longest_ORF_noncoding('ATGCGAATGTAGCATCA', 5)
#print longest_ORF_noncoding("GAGAGAGA", 1)
#print longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 1)
#print 'DNA input length', len('ATGCGAATGTAGCATCA')

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MetArg'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MetProAla'

        Returns nothing when there is no start codon.
    """
    index = 0
    i = 0
    start = dna[index:index+3]
    res = []
    if start == "ATG":
        res = ['Met']
        while i < len(dna)/3 - 1:
            index += 3
            codon = dna[index:index+3]
            if codon in ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']:
                res.append('Arg')
            if codon in ['GGT', 'GGC', 'GGA', 'GGG']:
                res.append('Gly')
            if codon in ['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG']:
                res.append('Ser')
            if codon in ['TGT', 'TGC']:
                res.append('Cys')
            if codon in ['TGG']:
                res.append('Typ') #What do I do here?
            if codon in ['TAT', 'TAC']:
                res.append('Tyr')
            if codon in ['CAT', 'CAC']:
                res.append('His')
            if codon in ['CAA', 'CAG']:
                res.append('Gin')
            if codon in ['AAT', 'AAC']:
                res.append('Asn')
            if codon in ['AAA', 'AAG']:
                res.append('Lys')
            if codon in ['GAT', 'GAC']:
                res.append('Asp')
            if codon in ['GAA', 'GAG']:
                res.append('Glu')
            if codon in ['GCT', 'GCC', 'GCA', 'GCG']:
                res.append('Ala')
            if codon in ['ACT', 'ACC', 'ACA', 'ACG']:
                res.append('Thr')
            if codon in ['CCT', 'CCC', 'CCA', 'CCG']:
                res.append('Pro')
            if codon in ['TTT', 'TTC']:
                res.append('Phe')
            if codon in ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']:
                res.append('Leu')
            if codon in ['ATT', 'ATC', 'ATA']:
                res.append('Ile')
            if codon in ['GTT', 'GTC', 'GTA', 'GTG']:
                res.append('Val')
            else:
                pass
            i += 1
    else:
        pass
    return ''.join(res) #joins the items in the list

# print coding_strand_to_AA('ATGCGA')
# print len(coding_strand_to_AA("ATGCGA"))/3
#print coding_strand_to_AA("AGAGAGAGAGAGA")
#print coding_strand_to_AA("ATGCCCGCTTT")
#print coding_strand_to_AA("ATGCGA")

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
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500) #Should be set to 1500
    print threshold
    All_ORFs = find_all_ORFs_both_strands(dna)
    Filtered_ORFs = []
    primary_AA = []
    for i in All_ORFs:
        if len(i) >= threshold:
            Filtered_ORFs.append(i)
        else:
            pass
    #print Filtered_ORFs
    for i in Filtered_ORFs:
        primary_AA.append(coding_strand_to_AA(i))
    return primary_AA

DNA_sequence = raw_input("What DNA sequence?\n")
gene_finder("DNA_sequence")

# if __name__ == "__main__":
#     import doctest
#     doctest.run_docstring_examples(coding_strand_to_AA, globals())
#>>>>>>> daf89177ec30bf2d3eb6f14c538e0a94cf857a77