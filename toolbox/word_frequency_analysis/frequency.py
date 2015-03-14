""" 
Kelly Brennan
Software Design
Professors Paul Ruvolo and Ben Hill 
Spring 2015

Analyzes the word frequencies in a book downloaded from
Project Gutenberg 
"""

import string

def get_word_list(file_name):
	""" Reads the specified project Gutenberg book.  Header comments,
		punctuation, and whitespace are stripped away.  The function
		returns a list of the words used in the book as a list.
		All words are converted to lower case. 

	>>> get_word_list('doc_test.txt')
	['cat', 'dog', 'moo']
	"""
	full_text = open(file_name)
	lines = full_text.readlines()
	line_begin = 0
	line_end = len(lines)-1
	while lines[line_begin].find('START OF THIS PROJECT GUTENBERG EBOOK') == -1:
		line_begin += 1
	while lines[line_end].find('END OF THIS PROJECT GUTENBERG EBOOK') == -1:
		line_end -= 1
	text = lines[line_begin+1:line_end-1] #Plus and minus one to get rid of header comment
	words_list = []
	for line in text:
		if line == '\n':
			pass
		else: 
			for word in line.split():
				word = word.strip(string.punctuation + string.whitespace)
				word = word.lower()
				words_list.append(word)
	return words_list

def get_top_n_words(word_list, n):
	""" Takes a list of words as input and returns a list of the n most frequently
		occurring words ordered from most to least frequently occurring.

		word_list: a list of words (assumed to all be in lower case with no
					punctuation
		n: the number of words to return
		returns: a list of n most frequently occurring words ordered from most
				 frequently to least frequentlyoccurring

	>>> get_top_n_words(['a', 'a', 'b'], 2)
	[(2, 'a'), (1, 'b')]
	"""
	d = dict() 
	most_common = []
	for word in word_list:
		if word not in d:
			d[word] = 1
		else:
			d[word] += 1
	for key, value in d.items():
		most_common.append((value, key))
	most_common.sort(reverse = True)
	return most_common[0:n]

# print get_top_n_words(get_word_list('persuasion.txt'), 10)

if __name__ == "__main__":
    import doctest
    doctest.testmod()