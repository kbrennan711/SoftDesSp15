"""
Mini-Project 3: Text Mining and Analysis
Kelly Brennan

Sonnet generator

"""

import sys
import pdb
import nltk
from pattern.web import *
from random import choice

sonnets_full_text = open('ShakespeareanSonnets.txt').read()
begin_index = sonnets_full_text.find('III')
end_index = sonnets_full_text.find('IV')
sonnets_text = sonnets_full_text[begin_index+1:end_index] 

def create_line_list(text):
	"""Creates a list of the lines in the text

	input: text
	output: list of lines in Shakespearean sonnets

	>>> create_line_list(sonnets_full_text[begin_index:2208])
	['III', 'Look in thy glass and tell the face thou viewest', 'Now is the time that face should form another;', 'Whose fresh repair if now thou not renewest,']
	"""
	all_lines = text.splitlines()
	line_list = []
	for item in all_lines:
		if item != '':
			line_list.append(item.strip(' '))
	return line_list


def last_words(text):
	""" Creates a list of the last words in each line

	input: list of lines in text
	output: list of the last word in each line_list with all letters lowercase and without punctuation

	>>> last_words(['Hello, how are you?', 'I am learning', 'about Text Mining and Analysis!', '', '', ''])
	['you', 'learning', 'analysis']
	"""
	last = []
	for line in text[:-3]: #Three spaces at the end of the sonnets
		if line in ['I', 'II', 'III', '  ']:
			pass
		else: 
			word = line.split()[-1]
			if word[-1] in ['.', '?', '!', ';', ':', ',', "'"]:
				last.append(word[:-1].lower())
			else:
				last.append(word.lower())
	return last

def beginning_of_line(line):
	""" Cuts the last word of the line

	input: line
	output: line without the last word

	>>> beginning_of_line('From fairest creatures we desire increase,')
	'From fairest creatures we desire'
	"""
	res = line.split()[:-1]
	return ' '.join(res)


def rhyme(word_input, level):
	"""
	Input: word_input = word; level = measure of how good the rhymes should be 
	(higher the level, the more the closer the phoneme of the words must be for function to return True)
	output: word and phoneme

	This does not include a doctest because the output is a VERY large list of words that rhyme with the input word.
	*Adapted from: http://stackoverflow.com/questions/25714531/find-rhyme-using-nltk-in-python 
	"""
	entries = nltk.corpus.cmudict.entries()
	syllables = [(word, syl) for word, syl in entries if word == word_input]
	rhymes = []
	for (word, syllable) in syllables:
		rhymes += [word for word, pron in entries if pron[-level:] == syllable[-level:]]
	return set(rhymes)	


def doTheyRhyme(word1, word2):
	"""Check if the words rhyme. This function ensures that the two word input are truly rhymes. It makes sure that the words are not the same

	Input: two words
	Output: Determination of whether or not the words rhyme
	*Adapted from: http://stackoverflow.com/questions/25714531/find-rhyme-using-nltk-in-python

	>>> doTheyRhyme('rat', 'bat')
	True
	>>> doTheyRhyme('color', 'cat')
	False
	>>> doTheyRhyme('necessary', 'unnecessary')
	False
	"""
	if word1 == word2:
		return False
	elif word1.find(word2) == len(word1) - len(word2):
		return False
	elif word2.find(word1) == len(word2) - len(word1):
		return False
	return word1 in rhyme(word2, 2)


def create_rhyming_dictionary(word_list):
	"""
	>>> create_rhyming_dictionary(['bat', 'cat', 'glue', 'rat', 'blue'])
	{'bat': ['bat', 'cat', 'rat'], 'glue': ['glue', 'blue']}
	>>> create_rhyming_dictionary(['hello', 'goodbye'])
	{'hello': ['hello'], 'goodbye': ['goodbye']}
	"""
	level = 1
	rhyming_dict = {word_list[0]:[word_list[0]]}
	for word in word_list[1:]:
		rhymed = False
		for key in rhyming_dict.keys():
			if doTheyRhyme(word, key):
				rhyming_dict.get(key).append(word)
				rhymed = True
		if rhymed == False:
			rhyming_dict[word] = [word]
	return rhyming_dict


rhyming_dictionary = create_rhyming_dictionary(last_words(create_line_list(sonnets_text)))


def rhyming_word(dictionary):
	""" Shakespearean sonnet rhyming pattern: abab cdcd efef gg
	input: rhyming_dictionary created
	ouput: a random key that contains a rhyming family

	doctest limited to extreme case because the key generated is random.
	>>> rhyming_word({'hello':['hello'], 'goodbye':['goodbye']})
	'Not Possible'
	"""
	useful_keys = []
	for key in dictionary:
		# print dictionary.get(key)
		if len(dictionary.get(key)) > 1:
			useful_keys.append(key)
	if useful_keys == []:
		return 'Not Possible'
	return choice(useful_keys)


def generate_sonnet(text_list, dictionary):
	""" Generates last two lines of Shakespearean sonnet (gg) randomly

	input: list of lines without last word, rhyming dictionary
	output: tuple of sonnet lines where each item is a line

	doctest limited to extreme case because of random generation of lines and rhyming words.
	>>> generate_sonnet(create_line_list(sonnets_text), {'hello':['hello'], 'goodbye':['goodbye']})
	'Sorry, sonnet generation not possible.'
	"""
	if rhyming_word(dictionary) == 'Not Possible':
		return 'Sorry, sonnet generation not possible.'
	a_key = rhyming_word(dictionary)
	a1 = choice(dictionary.get(a_key))
	a2 = choice(dictionary.get(a_key))
	while a1 == a2:
		a2 = choice(dictionary.get(a_key)) #a1 and a2 will always be different words
	l1 = beginning_of_line(choice(text_list[1:-1]))
	l2 = beginning_of_line(choice(text_list[1:-1]))
	while l1 == l2:
		l2 = beginning_of_line(choice(text_list[1:-1]))
	return [l1 + ' ' + a1, l2 + ' ' + a2]


print generate_sonnet(create_line_list(sonnets_text), rhyming_dictionary)


if __name__ == "__main__":
    import doctest
    doctest.testmod()