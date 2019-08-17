#!/usr/bin/env python
# coding: utf-8

# In[2]:


import camb
from camb import model, initialpower
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib import rc
import re
import matplotlib.pyplot as plt
from scipy.integrate import quad as qd
import numpy as np


# In[7]:


pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

def redshift_val(redshifts_values):
    pars.set_matter_power(redshifts=redshifts_values, kmax=2.0);

def PK(k, As, ns, alpha, beta):
    return As * (((k / 0.05) ** (ns - 1)) +  (alpha * (k / 0.05 )) + (beta * (((k / 0.05)**(-1))))) 

def prim_n(alpha, beta):
    pars.set_initial_power_function(PK, args=(2e-9, 0.96, alpha, beta));
    
def res():
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    return kh, z, pk


# In[8]:


def b(z):
    return np.sqrt(1+z)

def w(z, w_0, w_1):
    return w_0 + w_1*(1/1+z)

def H(z, omega_m_0, omega_x_0):
    func_1 = lambda z_prime_1 : (1 + w(z_prime_1, -1, 0))/(1 + z_prime_1)
    integ_1 = qd(func_1,0,z)[0]
    return H_0 * np.sqrt((omega_m_0*(1+z)**3) + (omega_x_0)*np.exp(3*integ_1))

def D(z, omega_m_0, omega_x_0):
    func_2 = lambda z_prime_2 : 1 / H(z_prime_2, omega_m_0, omega_x_0)
    integ_2 = qd(func_2, 0, z)[0]
    return (c/(1+z)) * integ_2

def Omega_m(z, omega_m_0, omega_lambda_0, omega_k_0, omega_r_0):
    return (omega_m_0 * (1+z)**3) * (omega_lambda_0 + omega_k_0*(1+z)**2 + omega_m_0*(1+z)**3 + omega_r_0*(1+z)**4)**(-1)

def gamma(z):
    return gamma_0 + gamma_1 * (z / 1+z)

def f_g(z):
    return (Omega_m(z, 0.321, 0.677, 0.001, 0.001)**gamma(z)) * (1+etha)

def G(z):
    func_3 = lambda z_prime_3 : f_g(z_prime_3) * (1+z_prime_3)**(-1)
    integ_3 = qd(func_3, 0, z)[0]
    return np.exp(integ_3)
    
def betta(z):
    return f_g(z) / b(z)

def sigma_r(z):
    return sigma_z/H(z, 0.321, 0.679)


# In[9]:


def P_k(z, alpha, beta):
    redshift_val([z])
    prim_n(alpha, beta)
    kh, z, pk = res()
    return kh, pk

def P_obs(k, z, mu, pk):
    pk_o = []
    for i, j in enumerate(k):
        pkval = kh, (b(z)**2) * ((1 + betta(z) * mu**2)**2) * np.exp(-(j**2)*(mu**2)*((sigma_r(z)**2)+(sigma_v**2))) * pk[i]
        pk_o.append(pkval)
    return kh, pk_o


# In[10]:


redshifts_values = np.arange(0, 1, 0.01)
alpha = 0
beta = 0
z = 0.2
mu = 0.2
H_0 = 67.5
c = 2.1
gamma_0 = 0.545
gamma_1 = 0
etha = 0
sigma_z = 0.001
sigma_v = 3
# def P_com():
#     kh,pk = P_k(z, alpha, beta)
#     _,pk_o = P_obs(kh, z, mu)
    
#     plt.loglog(kh, pk_o)
#     plt.show()

# P_com()    
kh, pk = P_k(z, alpha, beta)
kh, pk_o = P_obs(kh, z, mu, pk[0])
# pk = np.array(pk[0])
# pk_o = np.array(pk_o)
# delta = pk - pk_o
print(pk[0], '>>>>>>>>>>><<<<<<<<<<<<',pk_o)


# In[ ]:




