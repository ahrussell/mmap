from sage.all_cmdline import *
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler as DGSL

from mmp import MMP
from util import *
import norms

import random as rand

class GGH(MMP):
    @staticmethod
    def set_params(lam, k):
        n = pow(2, ceil(log(lam**2 * k)/log(2))) # dim of poly ring, closest power of 2 to k(lam^2)
        q = next_prime(ZZ(2)**(8*k*lam) * n**k, proof=False) # prime modulus

        sigma = int(sqrt(lam * n))
        sigma_prime = lam * int(n**(1.5))

        return (n, q, sigma, sigma_prime, k)

    @profile(LOG, "setup")
    def __init__(self, params, asym=False):

        (self.n, self.q, sigma, self.sigma_prime, self.k) = params

        S = PolynomialRing(ZZ, 'x')
        self.R = S.quotient_ring(S.ideal(x**self.n + 1))

        Sq = PolynomialRing(Zmod(self.q), 'x')
        self.Rq = Sq.quotient_ring(Sq.ideal(x**self.n + 1))

        # draw z_is uniformly from Rq and compute its inverse in Rq
        if asym:
            self.z = [self.Rq.random_element() for i in range(self.k)]
            self.zinv = [z_i**(-1) for z_i in self.z]
        else: # or do symmetric version
            z = self.Rq.random_element()
            self.z = [z for i in range(self.k)]

            zinv = z**(-1)
            self.zinv = [zinv for z_i in self.z]

        Sk = PolynomialRing(QQ, 'x')
        K = Sk.quotient_ring(Sk.ideal(x**self.n + 1))

        self.D_sigma = lambda: self.Rq(list(DGSL(ZZ**self.n, sigma)()))
        self.D_sigmap = lambda center: self.Rq(list(DGSL(ZZ**self.n, self.sigma_prime, c=center)()))

        # draw g (in Rq) repeatedly from a Gaussian distribution of Z^n (with param sigma)
        # until g^(-1) in QQ[x]/<x^n + 1> is small (< n^2)
        # while True:
        #     l = self.D_sigma
        #     ginv_K = K(l)**(-1)
        #     ginv_size = vector(ginv_K).norm()

        #     if ginv_size < self.n**2:
        #         self.g = self.Rq(l)
        #         self.ginv = g**(-1)
        #         break

        # don't check if g^(-1) in K is small because inverting g in K is expensive
        # and it's probably small anyway
        self.g = self.D_sigma()
        self.ginv = self.g**(-1)

        # compute zero-testing parameter p_zt
        # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        self.h = self.Rq(list(DGSL(ZZ**self.n, round(sqrt(self.q)))()))

        # create p_zt
        self.p_zt = self.ginv * self.h * prod(self.z)

    def encode(self, m, S):
        ''' encodes a vector m (in Zmod(q)^n) to index set S '''

        m = vector(Zmod(self.q),m)
        c = self.Rq(list(self.D_sigmap(self.g) + m))

        zinv = prod([self.zinv[i] for i in S])

        return self.Rq(c * zinv)

    def sample(self,i):
        # draw an element of Rq from a Gaussian distribution of Z^n (with param sigmaprime)
        # multiply by z_i^(-1) (or z^(-1) in symmetric version)

        return self.D_sigmap(zero_vector(self.n)) * self.zinv[i]

    def zero(self,i):
        ''' Level-1 encoding of 0 of index i'''
        return self.g * self.zinv[i]

    def is_zero(self, c):
        w = self.Rq(c) * self.p_zt

        return (norms.linf(w,self.q) < ZZ(RR(self.q)**(.75)))
