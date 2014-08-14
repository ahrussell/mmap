from sage.all_cmdline import *

import time

current_time = lambda: time.time()

class MMP():
    def __init__(self, lam, k):
        self.n = lam**2 * k # dimension of polynomial ring
        self.m = 2 # number of primes
        sigma = sqrt(lam*self.n)
        sigma_prime = lam * self.n * sqrt(self.n)

        primes = [131, 193] # list of primes
        self.x0 = prod(primes)
        # compute CRT coefficients
        coeff = [(self.x0 / p_i) * ZZ(Zmod(p_i)(self.x0 / p_i)**(-1)) for p_i in primes]

        S = PolynomialRing(ZZ, 'x')
        self.R = S.quotient_ring(S.ideal(x**self.n + 1))

        Sx0 = PolynomialRing(Zmod(self.x0), 'x')
        self.Rx0 = Sx0.quotient_ring(Sx0.ideal(x**self.n + 1))

        # draw z uniformly from Rq and compute its inverse in Rq
        self.z = self.Rx0.random_element()
        self.zinv = self.z**(-1)

        # c = current_time()
        # Sk = PolynomialRing(QQ, 'x')
        # K = Sk.quotient_ring(Sk.ideal([x**self.n + 1]))
        # # draw g (in Rq) repeatedly from a Gaussian distribution of Z^n (with param sigma)
        # # until g^(-1) in QQ[x]/<x^n + 1> is small (< n^2)
        # self.g = self.Rq.random_element()
        # self.ginv = self.g**(-1)

        # c = current_time()
        # # compute zero-testing parameter p_zt
        # # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        # h = self.Rq.random_element()

        # # create p_zt
        # self.p_zt = self.ginv * h

        # for i in range(k):
        #     self.p_zt *= self.zinv

    def sample(self):
        # draw an element from a Gaussian distribution of Z^n (with param sigmaprime)
        # convert this element into an element of R (or Rq)
        # multiply by z^(-1)

        return 0
        # return self.Rq.random_element() * self.zinv 

    def is_zero(self, c):
        # w = c * self.p_zt

        # return (max(w) < ZZ(self.q**(3/4)))
        return False

if __name__=="__main__":
    lam = 10 
    k = 5

    mmap = MMP(lam, k)
    encodings = [mmap.sample() for i in range(k)]

    result = 1

    t = current_time()
    for c in encodings:
        result *= c

    mmap.is_zero(result)
# 
