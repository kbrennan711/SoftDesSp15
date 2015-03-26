import random
import string

"""Testing levenshtein_distance function"""
word1 = 'HowDoYouDo'
word2 = 'Kellyz'

known = {}

# def levenshtein_distance(s1, s2):
#     print (s1,s2)
# 	if len(s1) == 0:
# 		return len(s2)
# 	elif len(s2) == 0:
# 		return len(s1)
# 	elif (s1, s2) in known:
# 		return known[(s1,s2)]
# 	known[(s1, s2)] = min([int(s1[0] != s2[0]) + levenshtein_distance(s1[1:], s2[1:]), 1 + levenshtein_distance(s1[1:], s2), 1 + levenshtein_distance(s1, s2[1:])])
# 	return known[(s1,s2)]


# print levenshtein_distance(word1, word2)
# print known


""" Testing TwoPointCrossover Function"""
def TwoPointCrossover(parent1, parent2):
	x = random.randint(1,min(len(parent1), len(parent2))/3)
	new1 = parent2[0:x] + parent1[x:-x] + parent2[-x:]
	new2 = parent1[0:x] + parent2[x:-x] + parent1[-x:]
	return (new1, new2)

TwoPointCrossover(word1, word2)

"""Testing Mutation function"""

VALID_CHARS = string.ascii_uppercase + " "
message = 'ASNAKXTGBKXQUKDRANLM'

def mutate_text(message, prob_ins=1.0, prob_del=0.0, prob_sub=0.0):
    """
    Given a Message and independent probabilities for each mutation type,
    return a length 1 tuple containing the mutated Message.

    Possible mutations are:
        Insertion:      Insert a random (legal) character somewhere into
                        the Message
        Deletion:       Delete one of the characters from the Message
        Substitution:   Replace one character of the Message with a random
                        (legal) character
    """
    if random.random() < prob_ins:
        print len(message)
        i = random.randint(0, len(message)-1)
        message.insert(i, random.choice(list(VALID_CHARS)))
        print len(message)
        return message
    elif random.random() < prob_del:
        i = random.randint(0, len(message)-1)
        print i
        print message.pop(i)
        return message
    elif random.random() < prob_sub:
        i = random.randint(0, len(message)-1)
        message.pop(i)
        message.insert(i, random.choice(list(VALID_CHARS))) 
        return 

    # TODO: Also implement deletion and substitution mutations
    # HINT: Message objects inherit from list, so they also inherit
    #       useful list methods
    # HINT: You probably want to use the VALID_CHARS global variable

    return (message, )   # Length 1 tuple, required by DEAP


print mutate_text(list(message))