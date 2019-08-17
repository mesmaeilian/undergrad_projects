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

a = [1]
pars.set_matter_power(redshifts=a, kmax=2.0)


results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
totCL=powers['total']
ls = np.arange(totCL.shape[0])[5:]


# for i,j in enumerate(a):

kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
# plt.loglog(kh, pk[i,:],'k')

# plt.xlabel('k/h Mpc')
# plt.show()






pars1 = camb.CAMBparams()
pars1.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars1.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars1.set_for_lmax(2500, lens_potential_accuracy=0);
pars1.set_matter_power(redshifts=a, kmax=2.0)


Powers = []
t = np.arange(0,1.1,0.1)

for q in t:

	def PK(k, As, ns):
	    return As * (((k / 0.05) ** (ns - 1)) +  (q * (k / 0.05 )) + (2 * (((k / 0.002)**(-1))))) 
	    # return  As * (k / 0.05) ** (ns - 1) * (1 + 0.1 * np.sin(10 * k))

	pars1.set_initial_power_function(PK, args=(2e-9, 0.96))
	results1 = camb.get_results(pars1)
	# powers1 =results1.get_cmb_power_spectra(pars1, CMB_unit='muK')

	# totCL_e=powers1['total']
	# totCL1 = [i[0] for i in totCL_e]
	# totCL1 = totCL1[5:]
	# ls1 = np.arange(totCL_e.shape[0])[5:]
	kh1, z1, pk1_e = results1.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
	pk1 = pk1_e[0,:]
	Powers.append(pk1)


fig, ax = plt.subplots()
# ax.semilogx(ls,totCL[:,0][5:],'k',linewidth=1.0)
ax.loglog(kh, pk[0,:], 'k', linewidth=1)
# plt.rc('font', family='serif')
# plt.rc('text', usetex=True)
fig.subplots_adjust(top=0.9, bottom=0.1,left=0.13,right=0.93)
ax.set_xlim(0.00005, 5)
# ax.set_ylim(np.min(Powers), np.max(Powers))


# ax.set_yticks([0, 2000, 4000, 6000, 8000, 10000, 12000])


line_seg = LineCollection([np.column_stack((kh1,power)) for power in Powers], linewidths=(1), linestyles='solid')
line_seg.set_array(t)
# line_seg.cmap = plt.cm.RdBu
# line_seg.cmap = plt.cm.PRGn
line_seg.cmap = plt.cm.BrBG
ax.add_collection(line_seg)
# ax.plot()
# ax.set_xscale('log')
ax.loglog()
# ax.set_xticks([0.0001,,0.001, 0.01, 0.1, 1],[1,2,3,4,5]) 
# ax.set_xticks([100,1000, 10000, 100000],[1,2,3,4]) 
# ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# x= np.arange(0,4,1)
# ax.semilogx()

plt.xlabel('k/h Mpc')
# plt.xlabel(r'${l}$')
# plt.ylabel(r'${l(l+1){C_l}^{TT}}$')

plt.text(0.03, 0.28,r'$P(k) = As\, [(\frac{k}{0.05})^{ns-1} + \alpha (\frac{k}{0.05}) + \beta (\frac{k}{0.002})^{-1}]$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
plt.text(0.03, 0.21,r'$As = 2^{-9} $', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
plt.text(0.03, 0.15,r'$ns = 0.96 $', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
plt.text(0.03, 0.09,r'$\beta = 2 $', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
plt.text(0.912, 0.542,r'$\alpha$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
cbaxes = fig.add_axes([0.85, 0.565, 0.02, 0.3]) 
fig.colorbar(line_seg, cax=cbaxes)



# plt.show()
file_name =  "/Users/Muhammad/Desktop/dgdmatterbeta20.png"
fig.savefig(file_name,format='png', dpi=1000)


