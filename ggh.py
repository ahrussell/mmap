from sage.all_cmdline import *

import time

current_time = lambda: time.time()

class MMP():
    def __init__(self, lam, k):
        print "set up rings"
        c = current_time()
        self.n = lam**2 * k
        self.q = Primes().next(ZZ(2**(8*k*lam) * self.n**k))
        sigma = sqrt(lam*self.n)
        sigma_prime = lam * self.n * sqrt(self.n)

        S = PolynomialRing(ZZ, 'x')
        self.R = S.quotient_ring(S.ideal(x**self.n + 1))

        Sq = PolynomialRing(Zmod(self.q), 'x')
        self.Rq = Sq.quotient_ring(Sq.ideal(x**self.n + 1))
        print "time: ", current_time() - c

        print "generate z"
        c = current_time()
        # draw z uniformly from Rq and compute its inverse in Rq
        self.z = self.Rq.random_element()
        self.zinv = self.z**(-1)
        print "time: ", current_time() - c

        print "generate g"
        c = current_time()
        Sk = PolynomialRing(QQ, 'x')
        K = Sk.quotient_ring(Sk.ideal(x**self.n + 1))
        # draw g (in Rq) repeatedly from a Gaussian distribution of Z^n (with param sigma)
        # until g^(-1) in QQ[x]/<x^n + 1> is small (< n^2)
        self.g = self.Rq.random_element()
        self.ginv = self.g**(-1)
        print "time: ", current_time() - c

        print "generate p_zt"
        c = current_time()
        # compute zero-testing parameter p_zt
        # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        h = self.Rq.random_element()

        # create p_zt
        self.p_zt = self.ginv * h

        for i in range(k):
            self.p_zt *= self.zinv
        print "time: ", current_time() - c

    def sample(self):
        # draw an element from a Gaussian distribution of Z^n (with param sigmaprime)
        # convert this element into an element of R (or Rq)
        # multiply by z^(-1)

        return self.Rq.random_element() * self.zinv 

    def is_zero(self, c):
        w = c * self.p_zt

        return (max(w) < ZZ(self.q**(3/4)))

if __name__=="__main__":
    lam = 50 
    k = 10

    mmap = MMP(lam, k)
    print "generate encodings"
    c = current_time()
    encodings = [mmap.sample() for i in range(k)]
    print "time: ", current_time() - c

    result = 1

    print "multiply encodings"
    t = current_time()
    for c in encodings:
        result *= c
    print "time: ", current_time() - t

    print "zero test"
    c = current_time()
    mmap.is_zero(result)
    print "time: ", current_time() - c
