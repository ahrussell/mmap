import clt
import time

current_time = lambda:time.time()

def linearcomb(omega, x0):
		n = len(omega)
		matrix = block_matrix([[identity_matrix(n), matrix(omega).transpose()],[matrix(ZZ, 1, n, 0), x0]])
		print matrix
		while True:
			temp = list(matrix.column(n))
			index = next((i for i, x in enumerate(temp) if x), None)
			minimum = temp[index]
			count = 0
			for i in range (self.n+1):
				if (temp[i] != 0):
					count+=1
				 	if(temp[i] < minimum):
						minimum = temp[i]
						index = i
			if count == 1:
				break
			for i in range (n+1):
				if (i != index):
					self.matrix.add_multiple_of_row(i, index, -(matrix[i, n]//minimum))
			temp = list(matrix.column(n))
		matrix = matrix.delete_rows([index, n])
		matrix = matrix.delete_columns([n])
		matrix = matrix.LLL()
		print "Reduced Matrix"
		print matrix
		return 0

def main():
	lam = 20
	k = 5
	print "CLT (lambda="+str(lam)+", k="+str(k)+")"
	params = clt.CLT.set_params(lam, k)
	mmap = clt.CLT(params)
	num_encodings = 10
	omega = vector(ZZ, [0 for i in range(num_encodings)])
	for i in range(self.n):
		self.omega[i] = (encodings[i])
	mmap.run(k, of_zero = True)
   	c = current_time()
	print "Computing matrix reduction:"
	Attack.__linearcomb__(atk)
	print "Time:", current_time()-c

if __name__=="__main__":
	main()

# attack psuedocode

# generate a bunch of level-k encodings of 0 c_i
# multiply all c_i by p_zt yielding w_i for each c_i
# compute lattice L orthogonal to w = (w_1, ..., w_l)
    # run LLL/BKZ on basis B of L
# get all vectors {u} in B that are short enough
# compute lattice L_r orthogonal to L with basis B_r
    # run LLL/BKZ on basis B_r

# for r in reduced_B_r
    # compute s such that <s, r> = 0 (mod x_0)
    # if <s, w> = 0 (mod x_0)
        # continue to next loop iteration
    # compute p = gcd(<s, w>, x_0)
    # p should be a factor of x_0

# repeat algorithm until x_0 is fully factored
