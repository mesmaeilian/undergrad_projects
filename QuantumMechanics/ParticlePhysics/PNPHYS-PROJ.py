import numpy as np
import scipy.special.sph_harm as sph



c = 3 * 1e8
alpha = 1/137
hbar = 
es = 
a0 = 
uqm = 4*1e-30
dqm = 8*1e-30

def rmass(m1, m2):
	return m1*m2/(m1+m2)

def eng_lvl(n, l, m1, m2):
	rmass = rmass(m1, m2)
	coff = (n - (l+1/2) + np.sqrt((l+1/2)**2 - alpha**2))**(2)
	return rmass * (c)**2 * (1 + ((alpha**2)/coff))**(-1/2)

# def spr_prt(l, m, theta, phi):
# 	return sph(l, m, theta, phi)

# def sigma(l):
# 	coff = l + 1/2 - ((l+1/2)**2 - alpha**2)**(1/2)


# def rad_prt(n,l):
# 	np.exp(-r/((n - sigma(l))*a0))

