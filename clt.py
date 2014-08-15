from sage.all import *

import time

print "CLT multilinear map implementation using SAGE" 

current_time = lambda:time.time()

#Constants
hBits = 80  # length of h_i's
alpha = 80  # length of g_i's
N = 30      # number of primes
pBits = 400 # length of primes p_i
rho = 41    # length of r_i's
bound = 160 # bound to determine zero test 

class MMP():
    
    def __init__(self, kappa):

        self.x0 = ZZ(1)
        self.p = [0 for i in range(N)]
        
        print "generate p_i's and x0: "
       
        for i in range(N):
            self.p[i] = next_prime(ZZ.random_element(2**pBits))
            
        self.x0 = prod(self.p[i] for i in range(N))

        print "generate crtCoeff_i's: "

        self.crtCoeff = [0 for i in range(N)]
        for i in range(N):
            Q = self.x0/self.p[i]

            self.crtCoeff[i] = ZZ(Zmod(self.p[i])(Q)**(-1))
            self.crtCoeff[i] = self.crtCoeff[i]*Q

        print "generate the g_i's: "
        self.g = [0 for i in range(N)]
        for i in range(N):
            self.g[i] = next_prime(ZZ.random_element(2**alpha))

        print "generate z and zinv: "
        while True:
            self.z = ZZ.random_element(self.x0)  
            try:
                self.zinv = ZZ(Zmod(self.x0)(self.z)**(-1))
                break
            except ZeroDivisionError:
                error = 1

        print "generate y: "
        self.y = self.encrypt(1,rho,1)


        print "generate zero tester p_zt: "
        zkappa = 1
        self.p_zt = 0
        for i in range(kappa):
            zkappa = Zmod(self.x0)(zkappa*self.z)
        for i in range(N):
            t_res = ZZ(Zmod(self.p[i])(self.g[i])**(-1))
            t_res = ZZ(Zmod(self.p[i])(t_res*zkappa))*ZZ.random_element(2**hBits)*(self.x0/self.p[i])
            self.p_zt = self.p_zt + t_res
        self.p_zt = Zmod(self.x0)(self.p_zt)

    def encrypt(self,m,nSize,level):
        res = 0
        for i in range(N):
            res = res + (m + self.g[i]*ZZ.random_element(2**nSize))*self.crtCoeff[i]
        res = Zmod(self.x0)(res)
        for j in range(level):
            res = Zmod(self.x0)(res*self.zinv)
        return res

    def sample(self,k):
        m = ZZ.random_element(2**alpha)
        c = self.encrypt(m,rho,k)
        return Zmod(self.x0)(c)
    
    def is_zero(self,c):
        w = Zmod(self.x0)(c*self.p_zt)
        if self.numBits(w) < self.numBits(self.x0) - bound:
            return 0
        else:
            return 1

    def zero_test(self,val,deg):

       # for i in range(
       return 0
       
    def numBits(self,x):
        return len(ZZ(x).digits(2)) 

if __name__=="__main__":

        kappa = 5

        begin = current_time()

        print "setup"
        c = current_time()
        mmap = MMP(kappa)
        print "time:", current_time() - c

        print "generate level 1 encodings"
        c = current_time()
        encodings = [mmap.sample(1) for i in range(kappa)]
        print "time:", current_time() - c

        print "multiply encodings to get level kappa encoding"
        c = current_time()
        result = 1
        for e in encodings:
            result *= e
        print "time:", current_time() - c

        print "zero test"
        c = current_time()
        mmap.is_zero(result)
        print "time:", current_time() - c

        print "total time:", current_time() - begin

