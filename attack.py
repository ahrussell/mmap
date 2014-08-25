from sage.all import *
import clt
import time

current_time = lambda:time.time()

def linearcomb(omega, x0):
		n = len(omega)
		mat = block_matrix([[identity_matrix(n), matrix(omega).transpose()],[matrix(ZZ, 1, n, 0), x0]])
		while True:
			temp = list(mat.column(n))
			index = next((i for i, x in enumerate(temp) if x), None)
			minimum = temp[index]
			count = 0
			for i in range (n+1):
				if (temp[i] != 0):
					count+=1
				 	if(temp[i] < minimum):
						minimum = temp[i]
						index = i
			if count == 1:
				break
			for i in range (n+1):
				if (i != index):
					mat.add_multiple_of_row(i, index, -(mat[i, n]//minimum))
			temp = list(mat.column(n))
		mat = mat.delete_rows([index, n])
		mat = mat.delete_columns([n])
		mat = mat.LLL()
		print "Reduced Matrix"
		print mat
		return mat

def main():
	lam = 2
	k = 5
	print "CLT (lambda="+str(lam)+", k="+str(k)+")"
	params = clt.CLT.set_params(lam, k)
	mmap = clt.CLT(params)
	num_encodings = 10
	omega = vector(ZZ, [0 for i in range(num_encodings)])
	for i in range(num_encodings):
		omega[i] = Zmod(mmap.x0)(mmap.run(k, of_zero = True) * mmap.p_zt)
   	c = current_time()
	print "Computing matrix reduction:"
	u = linearcomb(omega, mmap.x0)
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
