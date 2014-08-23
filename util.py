from sage.all import *

import time
import sys

LOG = sys.stdout

def random_gauss(stddev, dim):
    return list(random_vector(Zmod(3*stddev), dim))

def current_time():
    return time.time()

def mod_near(a,b):
    ''' returns a mod b where -b/2 <= a < b/2 '''

    x = ZZ(Zmod(b)(a))

    if x < round(b / 2):
        return x
    else:
        return x - ZZ(b)

def profile(log_file, message):
    def timed(func):
        def wrap(*args, **kwargs):
            if type(log_file) != 'str':
                log_file.write(message + "\n")
                start = time.time()
                result = func(*args, **kwargs)
                log_file.write("time: " + str(time.time() - start) + "\n\n")
            else:
                with open(log_file, 'a') as f:
                    f.write(message + "\n")
                    start = time.time()
                    result = func(*args, **kwargs)
                    f.write("time: " + str(time.time() - start) + "\n")

            return result
        return wrap
    return timed


