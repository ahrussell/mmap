import random as rand

from sage.all_cmdline import *
from sage.all import ZZ
from sage.rings.all import RR, QQ, PolynomialRing, Zmod
from sage.rings.arith import next_prime
from sage.functions.all import log, ceil, sqrt
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector, zero_vector
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler as DGSL

from mmp import MMP
from util import *
import norms

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

        S, x = PolynomialRing(ZZ, 'x').objgen()
        self.R = S.quotient_ring(S.ideal(x**self.n + 1))

        Sq = PolynomialRing(Zmod(self.q), 'x')
        self.Rq = Sq.quotient_ring(Sq.ideal(x**self.n + 1))

        # draw z_is uniformly from Rq and compute its inverse in Rq
        if asym:
            z = [self.Rq.random_element() for i in range(self.k)]
            self.zinv = [z_i**(-1) for z_i in z]
        else: # or do symmetric version
            z = self.Rq.random_element()
            zinv = z**(-1)
            z, self.zinv = zip(*[(z,zinv) for i in range(self.k)])

        Sk = PolynomialRing(QQ, 'x')
        K = Sk.quotient_ring(Sk.ideal(x**self.n + 1))

        # set up some discrete Gaussians
        DGSL_sigma = DGSL(ZZ**self.n, sigma)
        self.D_sigma = lambda: self.Rq(list(DGSL_sigma()))

        self.D_sigmap = lambda center: self.Rq(list(DGSL(ZZ**self.n, self.sigma_prime, c=center)()))
        
        # sigmap centered at zero
        DGSL_sigmap_zero = DGSL(ZZ**self.n, self.sigma_prime)
        self.D_sigmap_zero = lambda: self.Rq(list(DGSL_sigmap_zero()))

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

        # set up a discrete Gaussian centered at 0
        g = vector(ZZ, mod_near_poly(self.g,self.q))
        DGSL_sigmap_g = DGSL(ZZ**self.n, self.sigma_prime, c=g)
        self.D_sigmap_g = lambda: self.Rq(list(DGSL_sigmap_g()))

        # compute zero-testing parameter p_zt
        # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        self.h = self.Rq(list(DGSL(ZZ**self.n, round(sqrt(self.q)))()))

        # create p_zt
        self.p_zt = self.ginv * self.h * prod(z)

    def encode(self, m, S):
        ''' encodes a vector m (in Zmod(q)^n) to index set S '''

        m = vector(Zmod(self.q),m)
        # print 1
        # B = (ZZ**self.n).matrix()
        # print 0
        # w = B.solve_left(vector(self.g))
        # print 2
        # if w in ZZ**self.n:
        #     print 3
        # else:
        #     print 4

        s = vector(Zmod(self.q),self.D_sigmap_g())
        c = self.Rq(list(s + m))

        zinv = prod([self.zinv[i] for i in S])

        return self.Rq(c * zinv)

    def sample(self,S):
        # draw an element of Rq from a Gaussian distribution of Z^n (with param sigmaprime)
        # then encode at index set S

        print "sample"
        return self.D_sigmap_zero() * prod([self.zinv[i] for i in S])

    def is_zero(self, c):
        w = self.Rq(c) * self.p_zt

        return (norms.linf(w,self.q) < ZZ(RR(self.q)**(.75)))
