import prec_float
from prec_float import prec_float as pf


print 1+1, 2*2
print pf(0,prec=212),pf(1.1,prec=212), pf(1.1,prec=212),pf("1.23",prec=212), pf(1.1,prec=212)
a=pf(0,prec=212)
b = pf(1,prec=212)
c = pf("1.234567891234567891234566789",prec=212)
a=b *c *c *c
e = c**2 - c*c 
print a,e

b = a+a

print b

b = a+1

print b

b = 1+a

print b

print pf(2,prec=212)*a 
print 2*a 
print a*2 
print [a,b]
print [a,b]+[a,b]

