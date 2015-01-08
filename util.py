import time
import sys

from sage.rings.all import ZZ

LOG = sys.stdout

def mod_near_poly(x,q):
    ''' mod_near the coefficients of x (mod q) '''

    return map(lambda y: mod_near(y,q), x)

def mod_near(a,b):
    ''' returns a mod b where -b/2 <= a < b/2 '''
    a = ZZ(a)
    b = ZZ(b)

    return a - b*((2*a+b) // (2*b))

def profile(log_file, message):
    def timed(func):
        def wrap(*args, **kwargs):
            if type(log_file) != 'str':
                start = time.time()
                result = func(*args, **kwargs)
                log_file.write(message+": " + str(time.time() - start) + "\n")
            else:
                with open(log_file, 'a') as f:
                    f.write(message + "\n")
                    start = time.time()
                    result = func(*args, **kwargs)
                    f.write("time: " + str(time.time() - start) + "\n")

            return result
        return wrap
    return timed


