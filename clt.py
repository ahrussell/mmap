import math
import random as rand

from sage.rings.all import ZZ, Zmod
from sage.misc.misc_c import prod
from sage.crypto.util import random_prime
from sage.modules.free_module_element import vector, zero_vector

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

        n = (lam**2)*k # number of primes
        return (alpha, beta, rho, rho_f, eta, bound, n, k)

    @profile(LOG, "setup")
    def __init__(self, params, asym=False):
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

        # generate zs and zs^(-1)
        z = []
        zinv = []

        # generate z and z^(-1)
        if not asym:
            while True:
                z = ZZ.random_element(self.x0)  
                try:
                    zinv = ZZ(Zmod(self.x0)(z)**(-1))
                    break
                except ZeroDivisionError:
                    ''' Error occurred, retry sampling '''

            z, self.zinv = zip(*[(z,zinv) for i in range(self.k)])

        else: # asymmetric version
            for i in range(self.k):
                while True:
                    z_i = ZZ.random_element(self.x0)  
                    try:
                        zinv_i = ZZ(Zmod(self.x0)(z_i)**(-1))
                        break
                    except ZeroDivisionError:
                        ''' Error occurred, retry sampling '''

                z.append(z_i)
                zinv.append(zinv_i)

            self.zinv = zinv

        # generate p_zt
        zk = Zmod(self.x0)(1)
        self.p_zt = 0
        for z_i in z:
            zk *= Zmod(self.x0)(z_i)
        for i in range(self.n):
            self.p_zt += Zmod(self.x0)(ZZ(Zmod(self.primes[i])(self.g[i])**(-1) * Zmod(self.primes[i])(zk)) * ZZ.random_element(2**self.beta) * (self.x0/self.primes[i]))

        self.p_zt = Zmod(self.x0)(self.p_zt)

    def encode(self, m, S):
        ''' encodes a vector m (in ZZ^n) to index set S '''

        c = Zmod(self.x0)(0)

        for i in range(self.n):
            r_i = ZZ.random_element(2**self.rho)
            c += Zmod(self.x0)((m[i] + self.g[i] * r_i) * self.coeff[i])

        zinv = prod([self.zinv[i] for i in S])

        return Zmod(self.x0)(c * zinv)

    def sample(self, S):
        ''' sample an element at index set S '''
        m = [ZZ.random_element(self.g[i]) for i in range(self.n)]

        return self.encode(m, S)

    def zero(self,S):
        ''' encoding of 0 at index S '''
        return self.encode(zero_vector(self.n), S)

    def is_zero(self,c):
        w = abs(mod_near(c*self.p_zt, self.x0))
        return w < (self.x0 >> self.bound)

