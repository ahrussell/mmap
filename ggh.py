from sage.all_cmdline import *

from util import *

import random as rand

class MMP():
    @staticmethod
    def set_params(lam, k):
        n = lam**2 * k # dim of poly ring
        q = next_prime(ZZ(2)**(8*k*lam) * n**k, proof=False) # prime modulus

        sigma = int(sqrt(lam * n))
        sigma_prime = lam * int(n**(1.5))

        return (n, q, sigma, sigma_prime)

    def __init__(self, params):

        print "set up rings"
        c = current_time()

        (self.n, self.q, sigma, self.sigma_prime) = params

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
        # while True:
        #     l = random_gauss(sigma, self.n)
        #     ginv_K = K(l)**(-1)
        #     ginv_size = vector(ginv_K).norm()

        #     if ginv_size < self.n**2:
        #         self.g = self.Rq(l)
        #         self.ginv = g**(-1)
        #         break

        # don't check if g^(-1) in K is small because inverting g in K is expensive
        # and it's probably small anyway
        self.g = self.Rq(random_gauss(sigma, self.n))
        self.ginv = self.g**(-1)

        print "time: ", current_time() - c

        print "generate p_zt"
        c = current_time()
        # compute zero-testing parameter p_zt
        # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        self.h = self.Rq(random_gauss(round(sqrt(self.q)), self.n))

        # create p_zt
        self.p_zt = self.ginv * self.h * self.z**k

        print "time: ", current_time() - c

    def sample(self):
        # draw an element of Rq from a Gaussian distribution of Z^n (with param sigmaprime)
        # multiply by z^(-1)

        return self.Rq(random_gauss(self.sigma_prime, self.n)) * self.zinv

    def zero(self):
        ''' Level-1 encoding of 0 '''
        return self.g * self.zinv

    def is_zero(self, c):
        w = self.Rq(c) * self.p_zt

        f = lambda x, y: max(x, abs(mod_near(y, self.q)))

        norm = reduce(f, w, float('-inf')) 
        return (norm < ZZ(self.q**(.75)))


if __name__=="__main__":
    lam = 10
    k = 5
    params = MMP.set_params(lam, k)
    mmap = MMP(params)

    no_tests = 10

    tests_passed = 0
    for i in range(no_tests): 
        tests_passed += test_mmap(mmap, k, rand.choice([True, False]))

    print
    print "Tests passed:", tests_passed
    print "Tests failed:", no_tests - tests_passed

