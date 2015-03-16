class DNASequence(object):
    """ Represents a sequence of DNA """
    def __init__(self, nucleotides):
        """ constructs a DNASequence with the specified nucleotides.
             nucleotides: the nucleotides represented as a string of
                          capital letters consisting of A's, C's, G's, and T's """
        self.nucleotides = nucleotides
 
    def __str__(self):
        """ Returns a string containing the nucleotides in the DNASequence
        >>> seq = DNASequence("TTTGCC")
        >>> print seq
        TTTGCC
        """
        return self.nucleotides

    def get_reverse_complement(self):
        """ Returns the reverse complement DNA sequence represented
            as an object of type DNASequence

            >>> seq = DNASequence("ATGC")
            >>> rev = seq.get_reverse_complement()
            >>> print rev
            GCAT
            >>> print type(rev)
            <class '__main__.DNASequence'>
            """
        rev_comp = []
        for nucleotide in self.nucleotides:
            if nucleotide == "A":
                rev_comp.append('T')
            elif nucleotide == "T":
                rev_comp.append('A')
            elif nucleotide == "G":
                rev_comp.append('C')
            elif nucleotide == "C":
                rev_comp.append('G')
            else:
                rev_comp.append(' ')
        return ''.join(rev_comp)[::-1]

    def get_proportion_ACGT(self):
        """ Computes the proportion of nucleotides in the DNA sequence
            that are 'A', 'C', 'G', and 'T'
            returns: a dictionary where each key is a nucleotide and the
                corresponding value is the proportion of nucleotides in the
            DNA sequence that are that nucleotide.
            (NOTE: this doctest will not necessarily always pass due to key
                    re-ordering don't worry about matching the order)
        # >>> seq = DNASequence("AAGAGCGCTA")
        # >>> d = seq.get_proportion_ACGT()
        # >>> print (d['A'], d['C'], d['G'], d['T'])
        # (0.4, 0.2, 0.3, 0.1)
        """
        prop_dict = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}
        for nucleotide in self.nucleotides:
            if nucleotide in prop_dict:
                prop_dict[nucleotide] += 1.0/len(nucleotides)       
        return prop_dict

# ing empty list
#     for c in s:
#         if c not in d:
#             d[c] = 1
#         else:
#             d[c] += 1
#     return d

seq = DNASequence("AAGAGCGCTA")
d = seq.get_proportion_ACGT()
print (d['A'], d['C'], d['G'], d['T'])

# seq = DNASequence("TTTGCC")
# print seq
# rev = seq.get_reverse_complement()
# print rev
# print type(rev)
#<class '__main__.DNASequence'>

# if __name__ == '__main__':
#     import doctest
#     doctest.testmod()
