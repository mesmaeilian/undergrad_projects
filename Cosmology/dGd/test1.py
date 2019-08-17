import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower


cmap = plt.cm.PRGn

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
totCL=powers['total']
ls = np.arange(totCL.shape[0])[5:]


plt.plot(ls,totCL,'r')
plt.show()

