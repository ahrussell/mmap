from sage.all import *

import time

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
            with open(log_file, 'a') as f:
                f.write(message + "\n")
                start = time.time()
                result = func(*args, **kwargs)
                f.write("time: " + str(time.time() - start) + "\n")

            return result
        return wrap
    return timed

def test_mmap(mmap, k, with_zero = False):
    begin = current_time()

    print "generate encodings"
    c = current_time()
    encodings = [mmap.sample() for i in range(k - with_zero)]
    print "time: ", current_time() - c

    result = 1

    print "multiply encodings"
    t = current_time()
    for c in encodings:
        result *= c

    if with_zero:
        result *= mmap.zero()

    print "time: ", current_time() - t

    print "zero test"
    c = current_time()
    is_0 = mmap.is_zero(result)
    print "time: ", current_time() - c

    print "total time: ", current_time() - begin

    return is_0 == with_zero # should always be True

