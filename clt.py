from sage.all import *

import time

print "CLT multilinear map implementation using SAGE" 

current_time = lambda:time.time()

# #Constants
# hBits = 80  # length of h_i's
# alpha = 80  # length of g_i's
# N = 540      # number of primes
# pBits = 1838 # length of primes p_i
# rho = 41    # length of r_i's
# bound = 160 # bound to determine zero test 

class MMP():
    
    def __init__(self, lam, k):

        # set parameters
        self.alpha = lam
        self.beta = lam
        self.rho = lam
        self.mu = lam + self.alpha + self.rho

        rho_f = k * (self.mu + self.rho + self.alpha + 2) + self.rho + 1

        self.eta = rho_f + self.alpha + 2*self.beta + lam + 8
        self.bound = self.eta - self.beta - rho_f - lam - 3

        self.n = self.eta * ZZ(round(log(float(lam),2.0)))
        print "N=", self.n
        print "eta=", self.eta
        print "bound", self.bound

        self.x0 = ZZ(1)
        
        print "generate primes"
        primes = [random_prime(2**self.eta, proof=False) for i in range(self.n)]
        
        self.x0 = prod(primes)

        print "generate crtCoeff_i's: "

        self.coeff = [ZZ((self.x0/p_i) * ZZ(Zmod(p_i)(self.x0/p_i)**(-1))) for p_i in primes]

        print "generate the g_i's: "
        self.g = [random_prime(2**self.alpha, proof=False) for i in range(self.n)]

        print "generate z and zinv: "
        while True:
            z = ZZ.random_element(self.x0)  
            try:
                self.zinv = ZZ(Zmod(self.x0)(z)**(-1))
                break
            except ZeroDivisionError:
                ''' Error occurred '''

        print "generate zero tester p_zt: "
        zk = Zmod(self.x0)(1)
        self.p_zt = 0
        for i in range(k):
            zk *= Zmod(self.x0)(z)
        for i in range(self.n):
            self.p_zt += Zmod(self.x0)(ZZ(Zmod(primes[i])(self.g[i])**(-1) * Zmod(primes[i])(zk)) * ZZ.random_element(2**self.beta) * (self.x0/primes[i]))

        self.p_zt = Zmod(self.x0)(self.p_zt)

    def sample(self):
        m = [ZZ.random_element(2**self.alpha) for i in range(self.n)]

        c = Zmod(self.x0)(0)
        for i in range(self.n):
            r_i = ZZ.random_element(2**self.rho)
            c += Zmod(self.x0)((m[i] + self.g[i] * r_i) * self.coeff[i])

        return Zmod(self.x0)(c * self.zinv)
    
    def is_zero(self,c):
        num_bits = lambda x: len(ZZ(x).digits(2))

        w = Zmod(self.x0)(c*self.p_zt)
        if num_bits(w) < num_bits(self.x0) - self.bound:
            return 0
        else:
            return 1

if __name__=="__main__":

        lam = 52
        k = 5

        begin = current_time()

        print "setup"
        c = current_time()
        mmap = MMP(lam,k)
        print "time:", current_time() - c

        print "generate level 1 encodings"
        c = current_time()
        encodings = [mmap.sample() for i in range(k)]
        print "time:", current_time() - c

        print "multiply"
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

