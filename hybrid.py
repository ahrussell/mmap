from sage.all_cmdline import *

import time

current_time = lambda: time.time()

class MMP():
    @staticmethod
    def set_params(lam, k):
        return (int(sqrt(lam*lam*k)), int(sqrt(lam)))


    def __init__(self, params):
        # n = dim of poly ring, m = number of primes
        (self.n, self.m) = params

        primes = [random_prime(2**(lam*k*8), proof=False) for i in range(self.m)] # list of primes
        self.x0 = prod(primes)
        # compute CRT coefficients
        self.coeff = [ZZ((self.x0 / p_i) * ZZ(Zmod(p_i)(self.x0 / p_i)**(-1))) for p_i in primes]

        S = PolynomialRing(ZZ, 'x')
        self.R = S.quotient_ring(S.ideal(x**self.n + 1))

        Sx0 = PolynomialRing(Zmod(self.x0), 'x')
        self.Rx0 = Sx0.quotient_ring(Sx0.ideal(x**self.n + 1))

        Sprimes = [PolynomialRing(Zmod(p_i), 'x') for p_i in primes]
        self.Rprimes = [Spi.quotient_ring(Spi.ideal(x**self.n + 1)) for Spi in Sprimes]

        # draw z and its inverse, defined via CRT
        a = [Rpi.random_element() for Rpi in self.Rprimes]
        self.z = sum([self.Rx0(a[i]) * self.coeff[i] for i in range(self.m)])
        self.zinv = sum([self.Rx0(a[i]**(-1)) * self.coeff[i] for i in range(self.m)])

        # draw g (in Rq) repeatedly from a Gaussian distribution of Z^n (with param sigma)
        # until g^(-1) in QQ[x]/<x^n + 1> is small (< n^2)
        b = [Rpi.random_element() for Rpi in self.Rprimes]
        self.g = sum([self.Rx0(b[i]) * self.coeff[i] for i in range(self.m)])
        self.ginv = sum([self.Rx0(b[i]**(-1)) * self.coeff[i] for i in range(self.m)])

        # compute zero-testing parameter p_zt
        # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        h = self.Rx0.random_element()

        # create p_zt
        self.p_zt = self.Rx0(self.ginv * h)

        for i in range(k):
            self.p_zt *= self.z

    def sample(self):
        # draw an element from a Gaussian distribution of Z^n (with param sigmaprime)
        # convert this element into an element of R (or Rq)
        # multiply by z^(-1)

        a = [Rpi.random_element() for Rpi in self.Rprimes]
        c = sum([self.Rx0(a[i]) * self.coeff[i] for i in range(self.m)])

        return c * self.zinv

    def is_zero(self, c):
        w = c * self.p_zt

        return (max(w) < ZZ(self.x0**(3/4)))

if __name__=="__main__":
    lam = 60
    k = 5

    params = MMP.set_params(lam,k)

    print "setup"; c = current_time()
    mmap = MMP(params)
    encodings = [mmap.sample() for i in range(k)]
    print "time:", current_time() - c

    result = 1

    print "multiply"; t = current_time()
    for c in encodings:
        result *= c
    print "time:", current_time() - t

    print mmap.is_zero(result)
