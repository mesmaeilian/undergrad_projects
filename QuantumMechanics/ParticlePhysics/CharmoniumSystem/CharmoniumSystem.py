import numpy as np
import matplotlib.pyplot as plt

### defining a function to change eV to Joules
def eVtoJul(n):
	return n * 1.6*1e-19
### defining a function to change Joules to eV
def JultoeV(n):
	return n / (1.6*1e-19)

N = 5000
h_bar = eVtoJul(6.58*1e-16)
c = 3*1e8
L = 1
delta_r = L / N

### Mass of charmonium System
m = eVtoJul(3600*1e7)

### a and b are parameter of charmonium system potential
a=eVtoJul(0.48*1e6)
b=eVtoJul(0.18*1e6)

### Quantum Numbers
n = 3
l = 2

alpha = (h_bar**2)/(2*m)

### definig charmonium system potential
def V(i):
    return -a*h_bar*c/(delta_r*i+1e10) + b*(delta_r*i+1e10)/h_bar/c

### definig a function to solve schrodinger equation (with eigen-value problem)
def SH_eq(N, l):
	H = np.zeros([N,N])
	H[0][0] = alpha/delta_r**2 + V(0) + alpha*l*(l+1)/(delta_r)**2
	H[N-1][N-1] = alpha/delta_r**2 + V(0) + alpha*l*(l+1)/(delta_r)**2
	H[0][1] = -alpha/delta_r**2
	H[N-1][N-2] = -alpha/delta_r**2

	for i in range(1,N-1):
	    H[i][i-1] = -alpha/delta_r**2
	    H[i][i] = 2*alpha/delta_r**2 + V(i) + alpha*l*(l+1)/(delta_r*i)**2
	    H[i][i+1] = -alpha/delta_r**2

	return H
    
### Solving eigen value problem with numpy function
H = SH_eq(5000, l)
eig = np.linalg.eigh(H)

### Getting results
E = eig[0]
psi = eig[1]

print(JultoeV(E[n]))

