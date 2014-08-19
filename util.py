from sage.all import *

import time

def random_gauss(stddev, dim):
    return list(random_vector(Zmod(3*stddev), dim))

def current_time():
    return time.time()

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

