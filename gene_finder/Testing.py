"""
Testing
"""

def find_reverse_complement(dna):
	"""Return the reverse complement DNA sequence"""
	swap = {"A":"T", "T":"A", "C":"G", "G":"c"}
	return swap.join(dna)


find_reverse_complement("ATG")