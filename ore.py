from clt import CLT
from sage.all import *
import random

# Order revealing encryption

X = {}
Y = {}

X[0] = Matrix([[0,1,0,0,0],
			   [0,0,0,0,0],
			   [0,0,0,0,0],
			   [0,0,0,1,0],
			   [0,0,0,0,1]])

X[1] = Matrix([[0,0,1,0,0],
			   [0,0,0,0,0],
			   [0,0,0,0,0],
			   [0,0,0,1,0],
			   [0,0,0,0,1]])

Y[0] = Matrix([[0,0,0,0,0],
			   [1,0,0,0,0],
			   [0,0,0,0,1],
			   [0,0,0,1,0],
			   [0,0,0,0,1]])

Y[1] = Matrix([[0,0,0,0,0],
			   [0,0,0,1,0],
			   [1,0,0,0,0],
			   [0,0,0,1,0],
			   [0,0,0,0,1]])

e_1 = Matrix([[1,0,0,0,0]])


def random_number(b):
	''' returns a random b-bit number and its representation as a list of bits '''
	x_ = random.randint(0,2**4 - 1)
	x = map(int, list(bin(x_)[2:]))

	# padding with zeros
	while len(x) < b:
		x.insert(0,0)

	return x_, x

if __name__=="__main__":
	# run tests

	no_tests = 50
	passes = 0
	b = 16

	for i in range(no_tests):
		x_, x = random_number(b)
		y_, y = random_number(b)

		result = e_1

		# perform computation
		for j in range(len(x)):
			result *= X[x[j]]*Y[y[j]]

		fail = False

		if sum(result[0]) != 1:
			fail = True
		elif result[0][3] == 1:
			if x_ < y_:
				passes += 1
			else:
				fail = True
		elif result[0][0] == 1:
			if x_ == y_:
				passes += 1
			else:
				fail = True
		elif result[0][4] == 1:
			if x_ > y_:
				passes += 1
			else:
				fail = True

		if fail:
			print "Fail"
			print "Result:", result
			print "x:",x," y:",y
			print

	print "Passes: " + str(passes) + "/" + str(no_tests)
