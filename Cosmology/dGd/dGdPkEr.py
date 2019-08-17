import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib import rc
import re

cmap = plt.cm.PRGn

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
totCL=powers['total']
ls = np.arange(totCL.shape[0])[5:]


pars1 = camb.CAMBparams()
pars1.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars1.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars1.set_for_lmax(2500, lens_potential_accuracy=0);


Powers = []
Er = []
t = np.arange(0,2.1,0.2)

totCL = [i[0] for i in totCL]
totCL = totCL[5:]
for q in t:

	def PK(k, As, ns):
	    return As * (((k / 0.05) ** (ns - 1)) +  (q * (k / 0.05 )) + (0 * (((k / 0.002)**(-1))))) 
	    # return  As * (k / 0.05) ** (ns - 1) * (1 + 0.1 * np.sin(10 * k))

	pars1.set_initial_power_function(PK, args=(2e-9, 0.96))
	results1 = camb.get_results(pars1)
	powers1 =results1.get_cmb_power_spectra(pars1, CMB_unit='muK')

	totCL_e=powers1['total']
	totCL1 = [i[0] for i in totCL_e]
	totCL1 = totCL1[5:]
	A = np.array(totCL)
	B = np.array(totCL1)
	Er_i = B - A
	ls1 = np.arange(totCL_e.shape[0])[5:]
	# plt.semilogx(ls1,totCL1)
	Powers.append(totCL1)
	Er.append(Er_i)

alpha = 0
for i, j in enumerate(Er):

	a = 'alpha = %s' %alpha
	alpha += 0.2
	plt.plot(ls1, j, label= a)

plt.legend()
plt.show()





