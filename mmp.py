from util import *

import random as rand

class MMP():
    ''' Abstract MMP class '''

    def __init__(self):
        raise Exception("abstract class")

    @profile(LOG, "generate")
    def generate(self, l):
        ''' Generate l level-1 encodings '''
        return [self.sample() for i in range(l)]

    @profile(LOG, "multiply")
    def multiply(self, encodings, *args):
        ''' Multiplies encodings and args '''

        result = 1
        for c in encodings:
            result *= c

        return result

    def run(self, l, of_zero = False):
        ''' Generates l level-1 encodings, multiplies them'''
        ''' Will return a level-l encoding of 0 if of_zero=True '''

        encodings = self.generate(l - of_zero)

        if of_zero:
            encodings.append(self.zero())

        return self.multiply(encodings)

    def test(self, no_tests):

        passes = 0

        for i in range(no_tests):
            zero = rand.choice([True, False])
            result = self.run(self.k, of_zero = zero)

            passes += self.is_zero(result) == zero # add 1 if it passes

        return passes
