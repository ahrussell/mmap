from clt import CLT
from sage.all import *
import random

# Order revealing encryption

# util functions
def e(i):
	e_i = zero_vector(5)
	e_i[i-1] = 1

	return e_i

def bits(n,b):
	x = map(int, list(bin(n)[2:]))

	# padding with zeros
	while len(x) < b:
		x.insert(0,0)

	return x

def random_number(b):
	''' returns a random b-bit number and its representation as a list of bits '''
	x_ = random.randint(0,2**4 - 1)

	return x_, bits(x_,b)

def random_zmod_star(q):
	c = Zmod(q).random_element()

	while c == 0:
		c = Zmod(q).random_element()

	return c

I_5 = identity_matrix(5)


class ORE:

	def __init__(self, n):
		self.n = n

		# automaton matrices

		self.X = {}
		self.Y = {}

		self.X[0] = Matrix([[0,1,0,0,0],
					        [0,0,0,0,0],
					        [0,0,0,0,0],
					        [0,0,0,1,0],
					        [0,0,0,0,1]])

		self.X[1] = Matrix([[0,0,1,0,0],
					        [0,0,0,0,0],
					        [0,0,0,0,0],
					        [0,0,0,1,0],
					        [0,0,0,0,1]])

		self.Y[0] = Matrix([[0,0,0,0,0],
			      		    [1,0,0,0,0],
		     			    [0,0,0,0,1],
					        [0,0,0,1,0],
				     	    [0,0,0,0,1]])

		self.Y[1] = Matrix([[0,0,0,0,0],
			      		    [0,0,0,1,0],
					        [1,0,0,0,0],
					        [0,0,0,1,0],
					        [0,0,0,0,1]])

		# fix prime q
		self.q = random_prime(200)
		q = self.q

		GLq = GL(5, Zmod(q))
		Zq = Zmod(q)

		R = [e(1)]
		RInv = ["RInv_0"]

		for i in range(0,2*n-1):
			Ri = GLq.random_element()

			R.append(Ri)
			RInv.append(Ri.inverse())

		# R_2n = e_5
		R.append(e(5))
		RInv.append(e(5))

		self.secret_key = R, RInv # and index sets

	def encrypt(self, x):
		R = self.secret_key[0]
		RInv = self.secret_key[1]
		n = self.n

		x = bits(x,n)
		x.insert(0,0)

		X = []
		Y = []

		iden = identity_matrix(5)
		for i in range(1,n+1):
			a = random_zmod_star(self.q)* I_5
			b = random_zmod_star(self.q)* I_5

			X.append(a * R[2*i-2]*self.X[x[i]]*RInv[2*i-1])
			Y.append(b * R[2*i-1]*self.Y[x[i]]*RInv[2*i])

		# for all i, encode entries of X_i under index set S_i and entries of Y_i under index set T_i

		return (X,Y)

	def comp(self, c_x,c_y):
		X,_ = c_x
		_,Y = c_y

		z = I_5

		for i in range(n):
			z *=  X[i] * Y[i]

		return z 

if __name__=="__main__":

	no_tests = 500
	passes = 0
	n = 64

	ore = ORE(n)

	for i in range(no_tests):
		x = random.randint(0,2**n-1)
		y = random.randint(0,2**n-1)

		c_x = ore.encrypt(x)
		c_y = ore.encrypt(y)

		result = ore.comp(c_x,c_y)

		fail = False

		if result == 0:
			if x <= y:
				passes += 1
			else:
				fail = True
		else:
			if x > y:
				passes += 1
			else:
				fail = True
		
		if fail:
			print "Fail"
			print "Result:", result
			print "x:",x," y:",y
			print

	print "Passes: " + str(passes) + "/" + str(no_tests)
