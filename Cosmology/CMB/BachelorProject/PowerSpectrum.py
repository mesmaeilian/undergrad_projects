import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import math as m
from scipy.special import jv, jvp, spherical_jn
from scipy.integrate import quad
from scipy import interpolate



a_rec = 1/1000
a_eq = 1/24000
h0 = 67.5
h = 0.001
OME_M = 0.3
OME_L = 0.7
sigma = 0.03
ome_m = 0.13
ome_b = 0.02
phi = 1
K = 1



alpha = m.sqrt(a_rec/a_eq)
tr = m.sqrt(4*a_rec/(OME_M*(h0)**2))
k_min = 0.01/tr
k_max = 300/tr
print('tr is :',tr,'\nk_min is :',k_min,'\nk_max is :',k_max)
# k = K*tr
etha = tr*h
dp = -2*phi
dc = 3*dp/4
# vp = (-1/4)*(k*dp/etha)
# vc = vp
x_rec = (m.sqrt(1+alpha**2)-1)/alpha
x_s = 0.6*(ome_m**(1/4))*(ome_b**(-1/2))*(a_rec**(3/4))
k_d = 1/m.sqrt((2*x_s**2)+((sigma**2)*(x_rec**2)))
K_D = k_d/tr
print("K_D is :",K_D,"k_d is :",k_d)
k_min = 0.01
k_max = 300
# step = 0.1
K_MIN = k_min/tr
K_MAX = k_max/tr
print(K_MIN, K_MAX)
Powers = []
L = []
print("z is :",-dp/(etha*4))

for l in np.arange(0,100000,100):
	fig = lambda k : l * (l+1) * ( ( ( phi+dp/4 ) * ( jv(l,k) ) + ( (-1/4) *( k*dp/(etha) ) * jvp(l,k,1) ) ) **2 )*( np.exp(-(k/K_D)**2) ) * (k**2) * 2/np.pi
	res = quad(fig , K_MIN, K_MAX)
	print(l,">>>>>>>>",res[0])
	Powers.append(res[0])
	L.append(l)


# P_tck = interpolate.splrep(L, Powers, s=0)
# P_der = interpolate.splev(L, P_tck, der=1)
# plt.plot(L,P_der)
# L.reverse()
plt.plot(L,Powers,'.')
plt.show()

