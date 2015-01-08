import argparse
import random as rand

from ggh import GGH
from clt import CLT
from util import *

if __name__=="__main__":
    lam = 10
    k = 5
    default_implementation = "ggh"

    parser = argparse.ArgumentParser(description='''Run a multilinear map.
        You can specify the implementation and parameters used, as well as run some tests.''')

    parser.add_argument('-m', dest='implementation', choices=["clt", "ggh"], default=default_implementation, 
        help='The implementation to use. The default is ' + default_implementation + '.')
    parser.add_argument('-lam', dest='lam', default=lam, type=int, 
        help='Security parameter, default is ' + str(lam))
    parser.add_argument('-k', dest='k', default=k, type=int, 
        help='Multilinearity parameter, default is ' + str(k))
    parser.add_argument('--asym', dest='asym', action="store_true",
        help='Asymmetric version')
    parser.set_defaults(asymm=False)

    parser.add_argument('--tests', dest='tests', metavar='N', default=0, type=int, 
        help='Tests the maps and zero testing parameter with N tests')

    opts = parser.parse_args()

    lam = opts.lam
    k = opts.k
    asym = opts.asym
    maps = {"clt": CLT, "ggh": GGH}

    sym = "Asymmetric" if asym else "Symmetric"
    print opts.implementation.upper() + " (lambda="+str(lam)+", k="+str(k)+") - "+sym+"\n"

    M = maps[opts.implementation]
    params = M.set_params(lam, k)
    mmap = M(params, asym=asym)

    if opts.tests:
        passes = mmap.test_mmap(k, opts.tests)
        print "\nTests passed: " + str(passes) + "; tests failed: " + str(opts.tests - passes)
    else:
        mmap.run(k, of_zero = rand.choice([True, False]))



