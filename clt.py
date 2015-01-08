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
        return (alpha, beta, rho, rho_f, eta, bound, n, k_temp)

    def setup(self, params):
        # set parameters
        (self.alpha, self.beta, self.rho, self.rho_f, self.eta, self.bound, self.n, self.k) = params

        self.x0 = ZZ(1)
        
        self.primes = [random_prime(2**self.eta, lbound = 2**(self.eta - 1), proof=False) for i in range(self.n)]
        primes = self.primes
        
        self.x0 = prod(primes)

        # generate CRT coefficients
        self.coeff = [ZZ((self.x0/p_i) * ZZ(Zmod(p_i)(self.x0/p_i)**(-1))) for p_i in primes]

        # generate secret g_i
        self.g = [random_prime(2**self.alpha, proof=False) for i in range(self.n)]

    #@profile(LOG, "setup")
    def __init__(self, params):
        self.setup(params)

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
            self.p_zt += Zmod(self.x0)(ZZ(Zmod(self.primes[i])(self.g[i])**(-1) * Zmod(self.primes[i])(zk)) * ZZ.random_element(2**self.beta) * (self.x0/self.primes[i]))

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
        ''' level-1 encoding of 0 '''
        return self.encode([0 for i in range(self.n)], 1)
    
    def is_zero(self,c):
        w = abs(mod_near(c*self.p_zt, self.x0))
        return w < (self.x0 >> self.bound)

class CLTA(CLT):

    def __init__(self, params):
        self.setup(params) # common setup steps

        # generate zs and zs^(-1)
        zs = []
        zinvs = []

        for i in range(self.k):
            while True:
                z_i = ZZ.random_element(self.x0)  
                try:
                    zinv_i = ZZ(Zmod(self.x0)(z_i)**(-1))
                    break
                except ZeroDivisionError:
                    ''' Error occurred '''

            zs.append(z_i)
            zinvs.append(zinv_i)

        self.zinvs = zinvs

        # generate p_zt
        zk = Zmod(self.x0)(1)
        self.p_zt = 0
        for z in zs:
            zk *= Zmod(self.x0)(z)
        for i in range(self.n):
            self.p_zt += Zmod(self.x0)(ZZ(Zmod(self.primes[i])(self.g[i])**(-1) * Zmod(self.primes[i])(zk)) * ZZ.random_element(2**self.beta) * (self.x0/self.primes[i]))

        self.p_zt = Zmod(self.x0)(self.p_zt)

    def generate(self):
        m = [ZZ.random_element(self.g[i]) for i in range(self.n)]

        ss = []

        for i in range(self.k):
            ss.append(self.encode(m,[i]))

        return ss

    def zero(self):
        ''' level-k encoding of 0 '''

        m = [ZZ.random_element(self.g[i]) for i in range(self.n)]

        ss = []

        for i in range(self.k-1):
            ss.append(self.encode(m,[i]))

        ss.append(self.encode([0 for i in range(self.n)], [k-1]))

        return prod(ss)

    def encode(self, m, S):
        ''' encodes a vector m (in ZZ^n) to index set S '''

        c = Zmod(self.x0)(0)

        for i in range(self.n):
            r_i = ZZ.random_element(2**self.rho)
            c += Zmod(self.x0)((m[i] + self.g[i] * r_i) * self.coeff[i])

        zinv = prod([self.zinvs[i] for i in S])

        return Zmod(self.x0)(c * zinv)

    def run(self, k, of_zero = False):
        ''' Generates k level-1 encodings, multiplies them'''
        ''' Will return a level-k encoding of 0 if of_zero=True '''

        if of_zero:
            return self.zero()
        else:
            return self.multiply(self.generate())


if __name__=="__main__":
    lam = 20
    k = 5
    print "CLT (lambda="+str(lam)+", k="+str(k)+")"
    params = CLTA.set_params(lam, k)
    mmap = CLTA(params)

    print mmap.test_mmap(8,50)