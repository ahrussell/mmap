from sage.all_cmdline import *

from util import *

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
        while True:
            g = self.Rq(random_gauss(sigma, self.n))
            ginv_size = vector(K(list(g))**(-1)).norm()

            if ginv_size < self.n**2:
                self.ginv = g**(-1)
                break

        print "time: ", current_time() - c

        print "generate p_zt"
        c = current_time()
        # compute zero-testing parameter p_zt
        # randomly draw h (in Rq) from a discrete Gaussian with param q^(1/2)
        h = self.Rq(random_gauss(int(sqrt(self.q)), self.n))

        # create p_zt
        self.p_zt = self.ginv * h

        for i in range(k):
            self.p_zt *= self.z
        print "time: ", current_time() - c

    def sample(self):
        # draw an element of Rq from a Gaussian distribution of Z^n (with param sigmaprime)
        # multiply by z^(-1)

        return self.Rq(random_gauss(self.sigma_prime, self.n)) * self.zinv 

    def is_zero(self, c):
        w = c * self.p_zt

        return (max(w) < ZZ(self.q**(3/4)))

if __name__=="__main__":
    lam = 30 
    k = 5
    params = MMP.set_params(lam, k)

    mmap = MMP(params)
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

