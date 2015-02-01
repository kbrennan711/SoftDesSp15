"""
First Transversal loop 
"""

string = "Hello"

index = len(string)-1

while index >= 0:
	letter = string[index]
	print letter
	index = index - 1

#First for loop
prefixes = 'JKLMNOPQ'
suffix = 'ack'

for letter in prefixes:
	if letter == 'O' or letter == 'Q':
		print letter + 'u' + suffix
	else:
		print letter + suffix