from math import sqrt,log,cos,sin,pi
from random import random

Z = 79
e = 1.602e-19
E = 7.7e6*e
epsilon0 = 8.854e-12
a0 = 5.292e-11
sigma = a0/100
N = 100000
L = 100
l = 1
while l <= L:
    def gaussian():
        r = sqrt(-2*sigma*sigma*log(1-random()))
        theta = 2*pi*random()
        x = r*cos(theta)
        y = r*sin(theta)
        return x,y

    count = 0
    for i in range(N):
        x,y = gaussian()
        b = sqrt(x*x+y*y)
        if b<Z*e*e/(2*pi*epsilon0*E):
            count += 1
    print(count, "particles were reflected out of", N)
    l += 1

