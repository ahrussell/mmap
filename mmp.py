from util import *

import random as rand

class MMP():
    ''' Abstract MMP class '''

    def __init__(self):
        raise Exception("abstract class")

    #@profile(LOG, "generate")
    def generate(self, k):
        ''' Generate k level-1 encodings '''
        return [self.sample() for i in range(k)]

    #@profile(LOG, "multiply")
    def multiply(self, encodings, *args):
        ''' Multiplies encodings and args '''

        result = 1
        for c in encodings:
            result *= c

        return result

    def run(self, k, of_zero = False):
        ''' Generates k level-1 encodings, multiplies them'''
        ''' Will return a level-k encoding of 0 if of_zero=True '''

        encodings = self.generate(k - of_zero)

        if of_zero:
            encodings.append(self.zero())

        return self.multiply(encodings)

    def test_mmap(self, k, no_tests):

        passes = 0

        for i in range(no_tests):
            with_zero = rand.choice([True, False])
            result = self.run(k, of_zero = with_zero)

            passes += self.is_zero(result) == with_zero # add 1 if it passes

        return passes
