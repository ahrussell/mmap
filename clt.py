from sage.all import *

import math

from util import *

class MMP():
    
    @staticmethod
    def set_params(lam, k):
        alpha = lam # bitsize of g_i
        beta = lam # bitsize of h_i
        rho = lam # bitsize of r_i

        rho_f = k * (rho + alpha + 2) + rho + 1 # max bitsize of r_i at level-k
        eta = rho_f + alpha + 2*beta + lam + 8 # bitsize of primes p_i
        bound = eta - beta - rho_f - lam - 3 # bitsize of message to extract with p_zt

        n = 6*lam # number of primes

        return (alpha, beta, rho, eta, bound, n)

    def __init__(self, params):

        # set parameters
        (self.alpha, beta, self.rho, eta, self.bound, self.n) = params

        self.x0 = ZZ(1)
        
        print "generate primes"
        primes = [random_prime(2**eta, proof=False) for i in range(self.n)]
        
        self.x0 = prod(primes)

        print "generate crt coeff: "

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
            self.p_zt += Zmod(self.x0)(ZZ(Zmod(primes[i])(self.g[i])**(-1) * Zmod(primes[i])(zk)) * ZZ.random_element(2**beta) * (self.x0/primes[i]))

        self.p_zt = Zmod(self.x0)(self.p_zt)

    def encode(self,m,level):
        c = Zmod(self.x0)(0)

        for i in range(self.n):
            r_i = ZZ.random_element(2**self.rho)
            c += Zmod(self.x0)((m[i] + self.g[i] * r_i) * self.coeff[i])

        return Zmod(self.x0)(c * self.zinv**level)

    def sample(self):
        m = [ZZ.random_element(self.g[i]) for i in range(self.n)]

        return self.encode(m, 1)
    
    def is_zero(self,c):
        w = Zmod(self.x0)(c*self.p_zt)
        return w < (self.x0 >> self.bound)

if __name__=="__main__":

        lam = 52
        k = 5
        params = MMP.set_params(lam, k)

        begin = current_time()

        print "setup"
        c = current_time()
        mmap = MMP(params)
        print "time:", current_time() - c

        print "generate level 1 encodings"
        c = current_time()
        encodings = [mmap.sample() for i in range(k)]
        # zero = mmap.encrypt([0 for i in range(mmap.n)], 1)
        print "time:", current_time() - c

        print "multiply"
        c = current_time()
        result = 1
        for e in encodings:
            result *= e
        # result *= zero
        print "time:", current_time() - c

        print "zero test"
        c = current_time()
        is_0 = mmap.is_zero(result)
        print "time:", current_time() - c

        print "total time:", current_time() - begin

        print is_0
