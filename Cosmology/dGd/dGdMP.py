### Muhammad Sadegh Esmaeilian
### August 15, 2019
### Calculating Fisher Matrix for the galaxy power spectrum 

####################

### importing Libraries that we need in this code
import camb
from camb import model, initialpower
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib import rc
import re
import matplotlib.pyplot as plt
from scipy.integrate import quad as qd
import numpy as np
import time as t
# from tqdm import tqdm, tqdm_pandas
# import pandas as pd
# tqdm.pandas()

####################

### initializing some variables that we are going to use them soon
H_0 = 67.5
gamma_0 = 0.545
gamma_1 = 0
etha = 0
sigma_z = 0.001
sigma_v_0 = 300
c = 2.99*100000
f_sky = 0.3636
delta_z = 0.1
z_min = 0.65
z_max = 2.05
z_med = 0.9
z_0 = z_med/1.412
mu = np.arange(-1,1,0.01)

####################

### using CAMB to initialize cosmic parameters
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=1e-10, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results = camb.get_results(pars)

####################

### defining a function that set z values for us
def redshift_val(redshifts_values):
    pars.set_matter_power(redshifts=redshifts_values, kmax=2.0);

####################

### defining a function that create an arbitrary primordial power function with form of PK = As*(k^(ns-1) + a*(k) + b/(k))
def PK(k, As, ns, alpha, beta):
    return As * (((k / 0.05) ** (ns - 1)) +  (alpha * (k / 0.05 )) + (beta * (((k / 0.05)**(-1))))) 

####################

### defining a function that set our primordial power function
def prim_n(alpha, beta):
    pars.set_initial_power_function(PK, args=(1e-10, 0.96, alpha, beta));

####################

### defining a function that get the results (matter power spectrum)
def res():
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    return kh, z, pk

####################

### defining linear bias factor
def b(z):
    return np.sqrt(1+z)

####################

### defining Hubble Parameter, Angular Diameter and Matter Density Parameter that we get from our results (CAMB)
# Hubble Parameter
def H(z):
    res = results.hubble_parameter(z)
    return res

# Angular Diameter
def D(z):
    res = results.angular_diameter_distance(z)
    return res

# Matter Density Parameter
def Omega_m(z):
    a = results.get_Omega('baryon', z)
    b = results.get_Omega('nu', z)
    c = results.get_Omega('cdm', z)
    return a+b+c

####################

### defining γ that used in γ-parameterization
def gamma(z):
    return gamma_0 + gamma_1 * (z / 1+z)

####################

### defining growth rate
def f_g(z):
    return (Omega_m(z)**gamma(z)) * (1+etha)

####################

### defining growth function
def G(z):
    func_3 = lambda z_prime_3 : f_g(z_prime_3) * (1+z_prime_3)**(-1)
    integ_3 = qd(func_3, z, 0)[0]
    return np.exp(integ_3)
    
####################

### defining linear redshift-space distortion parameter
def betta(z):
    return f_g(z) / b(z)

####################

### defining error induced by spectroscopic redshift measurement 
def sigma_r(z):
    return sigma_z*c/H(z)

####################

### defining dispersion of pairwise peculiar velocities
def sigma_v(z):
    return sigma_v_0/H(z)

####################

### defining Galaxy Number Density
def n_z(z):
    return (z**2) * np.exp(-(z/z_0)**3/2)

####################

### defining V Survey
# Radial Coordinate
def r_z(z):
    func_4 = lambda z_prime_4 : H(z_prime_4)**(-1)
    integ_4 = qd(func_4, 0, z)[0]
    return integ_4*c

# V Survey
def vr_z(z):
    return (4*np.pi/3) * (f_sky) * ((r_z(z+(delta_z/2))**3) - (r_z(z-(delta_z/2))**3))

####################

### defining Matter Power Spectrum, Observational Matter Power Spectrum ...
# Matter Power Spectrum
def P_k(z, alpha, beta):
    redshift_val([z])
    prim_n(alpha, beta)
    kh, z, pk = res()
    return kh, pk

# Observational Matter Power Spectrum (TYPE1 : with alpha, beta)
def P_obs(z, mu, alpha, beta):
    pk_o = []
    kh, pk = P_k(z, alpha, beta)
    for i, j in enumerate(kh):
        pkval = (b(z)**2) * ((1 + betta(z) * mu**2)**2) * np.exp(-(j**2)*(mu**2)*((sigma_r(z)**2)+(sigma_v(z)**2))) * pk[0][i]
        pk_o.append(pkval)
    return kh, pk_o

# Observational Matter Power Spectrum (TYPE2 : from P matter)
def P_obs_z(kh, pk, z, mu):
    pk_o = []
    for i,j in enumerate(kh):
        pkval = (b(z)**2) * ((1 + betta(z) * mu**2)**2) * np.exp(-(j**2)*(mu**2)*((sigma_r(z)**2)+(sigma_v(z)**2))) * pk[0][i]
        pk_o.append(pkval)
    return kh, pk_o

####################

### defining a function calculating derivative by Specifying Parameter
def pod(z, coff, alpha, beta):
    epsilon = 0.001
    kh1, pk1 = P_k(z, alpha, beta)
    deriv = []
    if coff == 'alpha':
        alpha += epsilon
        _, pk2 = P_k(z, alpha, beta)
        for i,j in enumerate(kh1):
            p1 = pk1[0][i]
            p2 = pk2[0][i]
            der = (p2 - p1)/epsilon
            deriv.append(der)
    elif coff == 'beta':
        beta += epsilon
        _, pk2 = P_k(z, alpha, beta)
        for i,j in enumerate(kh1):
            p1 = pk1[0][i]
            p2 = pk2[0][i]
            der = (p2 - p1)/epsilon
            deriv.append(der)
            
    return kh1, pk1, deriv

####################

### defining a small function that give us Power Spectrum for each K 
def val_finder(k, kh, lst):
    index = np.where(kh == k)[0][0]
    return lst[index]
    
####################

### defining a function that calculate Fisher matrix elements (Specify Each element with coff Value; 'alpha' or 'beta' for diagonal elements)
def v(z, coff):
    alpha = 0.008
    beta = 0.002
    cof = (8*(np.pi)**2) ** (-1)
    v = vr_z(z)
    dk = 0.04628311744711677
    dmu = 0.05
    
    kh, pk1, derivA = pod(z, 'alpha', 0.008, 0.002)
    kh, pk1, derivB = pod(z, 'beta', 0.008, 0.002)
    
    if coff == "alpha" or coff == "beta":
        if coff == "alpha":
            c_lst = derivA
        else:
            c_lst = derivB
        mu_integ = []
        for i in mu:
            _, pk_o = P_obs_z(kh, pk1, z, i)
            func_var = list(map(lambda k: cof*v*((k**3)*(((n_z(z))**2)/((n_z(z)*val_finder(k,kh,pk_o)+1)**2))*(val_finder(k,kh,c_lst)**2)),kh))
            res = [i*dk for i in func_var]
            integ = np.sum(res)
            mu_integ.append(integ)
        mu_integ = [c*dmu for c in mu_integ]
        return np.sum(mu_integ)
    
    else:
        mu_integ = []
        for i in mu:
            _, pk_o = P_obs_z(kh, pk1, z, i)
            func_var = list(map(lambda k: cof*v*((k**3)*(((n_z(z))**2)/((n_z(z)*val_finder(k,kh,pk_o)+1)**2))*val_finder(k,kh,derivA)*val_finder(k,kh,derivB)),kh))
            res = [i*dk for i in func_var]
            integ = np.sum(res)
            mu_integ.append(integ)
        mu_integ = [c*dmu for c in mu_integ]
        return np.sum(mu_integ)

####################

### defining a function that calculate cov_matrix from Fisher matrix that we calculate from previous function
def cov_sq(cov_matr):
    sq_cov = np.array([[0,0],[0,0]], dtype=np.float64)
    for i,p in enumerate(cov_matr):
        for j,q in enumerate(p):
            sq_cov[i][j] = abs(cov_matr[i][j])**(1/2)                       
    return sq_cov   

####################

### defining a main function that give us the results using Fisher Matrix from previous functions  
def fc_matr(z):
    t_start = t.time()
    a11 = v(z, 'alpha')
    a12 = v(z, 'alphabeta')
    a21 = a12
    a22 = v(z, 'beta')
    fish_matr = np.array([[a11,a12],[a21,a22]])
    cov_matr = np.linalg.inv(fish_matr)
    sq_cov = cov_sq(cov_matr)
    t_finish = t.time()
    print("for z : %s cov_matrix is \n"%str(z), sq_cov)
    print('Run Time :', t_finish-t_start,'Seconds')
    return sq_cov

####################

### Average Run Time 
# μ steps : 0.01 ---> 50 Seconds
# μ steps : 0.02 ---> 35 Seconds
# μ steps : 0.05 ---> 20 Seconds





