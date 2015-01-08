from util import mod_near
from sage.all import ZZ,RR

def linf(v,q):
	''' l-inf norm (mod q) '''
	return reduce(lambda x,y: max(x, abs(mod_near(y,q))), v, float('-inf'))

def ln(v,n,q):
	''' l_n norm (mod q) '''
	return RR(sum([mod_near(x**n,q) for x in v]))**(1/RR(n))

def l2(v,q):
	return vector([mod_near(x,q) for x in v]).norm()
