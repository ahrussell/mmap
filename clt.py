from sage.all import *

import math
import random as rand

from mmp import MMP
from util import *

class CLT(MMP):
    
    @staticmethod
    def set_params(lam, k):
        alpha = lam # bitsize of g_i
        beta = lam # bitsize of h_i
        rho = lam # bitsize of r_i

        rho_f = k * (rho + alpha + 2) + rho + 1 # max bitsize of r_i at level-k
        eta = rho_f + alpha + 2*beta + lam + 8 # bitsize of primes p_i
        bound = eta - beta - rho_f - lam - 3 # bitsize of message to extract with p_zt

        n = 6*lam # number of primes
        k_temp = k
        return (alpha, beta, rho, eta, bound, n, k_temp)

    #@profile(LOG, "setup")
    def __init__(self, params):

        # set parameters
        (self.alpha, beta, self.rho, self.eta, self.bound, self.n, self.k) = params

        eta = self.eta
        self.x0 = ZZ(1)
        
        primes = [random_prime(2**eta, proof=False) for i in range(self.n)]
        
        self.x0 = prod(primes)

        # generate CRT coefficients
        self.coeff = [ZZ((self.x0/p_i) * ZZ(Zmod(p_i)(self.x0/p_i)**(-1))) for p_i in primes]

        # generate secret g_i
        self.g = [random_prime(2**self.alpha, proof=False) for i in range(self.n)]

        # generate z and z^(-1)
        while True:
            z = ZZ.random_element(self.x0)  
            try:
                self.zinv = ZZ(Zmod(self.x0)(z)**(-1))
                break
            except ZeroDivisionError:
                ''' Error occurred '''

        # generate p_zt
        zk = Zmod(self.x0)(1)
        self.p_zt = 0
        for i in range(self.k):
            zk *= Zmod(self.x0)(z)
        for i in range(self.n):
            self.p_zt += Zmod(self.x0)(ZZ(Zmod(primes[i])(self.g[i])**(-1) * Zmod(primes[i])(zk)) * ZZ.random_element(2**beta) * (self.x0/primes[i]))

        self.p_zt = Zmod(self.x0)(self.p_zt)

    def encode(self,m,level):
        ''' Encodes a vector m in ZZ^n '''

        c = Zmod(self.x0)(0)

        for i in range(self.n):
            r_i = ZZ.random_element(2**self.rho)
            c += Zmod(self.x0)((m[i] + self.g[i] * r_i) * self.coeff[i])

        return Zmod(self.x0)(c * self.zinv**level)

    def sample(self):
        m = [ZZ.random_element(self.g[i]) for i in range(self.n)]

        return self.encode(m, 1)

    def zero(self):
        return self.encode([0 for i in range(self.n)], 1)
    
    def is_zero(self,c):
        w = abs(mod_near(c*self.p_zt, self.x0))
        return w < (self.x0 >> self.bound)

if __name__=="__main__":
    lam = 20
    k = 5
    print "CLT (lambda="+str(lam)+", k="+str(k)+")"
    params = CLT.set_params(lam, k)
    mmap = CLT(params)

    mmap.run(k, of_zero = rand.choice([True, False]))
