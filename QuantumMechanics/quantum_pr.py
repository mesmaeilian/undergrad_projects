import numpy as np 
import matplotlib.pyplot as plt
from scipy.special import spherical_jn as sjn
from scipy.special import sph_harm as sph
from scipy.integrate import quad
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal


h = 1.0545718 * 10e-34
rad = 5.3 * 10e-11
m = 9.11*10e-31


def spbess_zero(l,n,n_digit):
	count = 0
	arg1 = 0
	step = 10**(-n_digit)
	arg2 = step
	while count != n:
		arg1 = arg2
		arg2 += step
		if sjn(l,arg1) * sjn(l,arg2) < 0:
			count += 1
			mean = (arg1+arg2)/2

	return round(mean,n_digit)

def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def eig_fun(n,l,ml,r,t,p):

	zero = spbess_zero(l,n,4)
	func = lambda r: (sjn(l,r*zero/rad)**2)*(r**2)
	norm_fact = quad(func, 0, rad)
	energy = (zero**2)*(h**2)/(2*m*(rad**2))
	k = 2 * m * energy / (h**2)


	psi = (norm_fact[0])*(sjn(l,r*zero/rad))*(sph(ml,l,t,p))
	print(psi)
	# psi = '%.2E' % Decimal(psi)
	return psi, energy

def q_3dplot():
	N = range(1,6)
	sub = str.maketrans("0123456789","₀₁₂₃₄₅₆₇₈₉")
	for n in N:
		for l in range(n):
			for ml in range(-l,l+1):
				if ml < 3 and ml > -3:
					print(n)
					theta = np.arange(0, 2*np.pi, 0.01)
					phi = np.arange(0, np.pi, 0.01)
					fig = plt.figure()
					ax = plt.axes(projection='3d')
					ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
					ax.view_init(45, 115)
					x, y = np.meshgrid(theta, phi)
					z_label = ('\u03A8'+str(n)+str(l)+str(ml)).translate(sub)
					z, energy = eig_fun(n,l,ml,rad/10,x,y)
					ax.contour3D(x, y, z,60, cmap='binary')
					ax.set_xlabel('\u03B8')
					ax.set_ylabel('\u03A6')
					ax.set_zlabel(z_label)
					ax.set_xlim(ax.get_xlim()[::-1])
					# textstr = '\n'.join((
		   #  		r'$\mathrm{n}=%.f$' % (n, ),
		   #  		r'$\mathrm{l}=%.f$' % (l, ),
		   #  		r'$\mathrm{ml}=%.f$' % (ml, ),
		   #  		r'$\mathrm{E}=%.f$' % (energy, )
		   #  		))
					# ax.text(0.88, 0.97, textstr, transform=ax.transAxes, fontsize=8,
		   #      	verticalalignment='top',horizontalalignment='left', bbox=dict(boxstyle='round', facecolor="wheat", alpha=0.5))
					file_name =  "/Users/Muhammad/Desktop/QP3D/Contour3D(BW)/PsiThetaPhi/" + "EigenFunc(nlm=%s).eps"%(str(n)+str(l)+str(ml))
					fig.savefig(file_name,format='eps', dpi=1000)
	for n in N:
		for l in range(n):
			for ml in range(-l,l+1):
				if ml < 3 and ml > -3:
					print(n)
					phi = np.linspace(0, np.pi,200)
					R = np.linspace(0,rad,80)
					fig = plt.figure()
					ax = fig.gca(projection='3d')
					ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
					ax.view_init(45, 115)
					x, y = np.meshgrid(R, phi)
					z_label = ('\u03A8'+str(n)+str(l)+str(ml)).translate(sub)
					z, energy = eig_fun(n,l,ml,x,0.01,y)
					
					# Z = [np.real(h[0]) for h in z]

					# norm = plt.Normalize(np.min(Z), np.max(Z))
					# colors = plt.cm.RdGy(norm(Z))
					# rcount, ccount = colors.shape
					surf = ax.plot_wireframe(x, y, z, rstride=4, cstride=4)
					# surf = ax.plot_surface(x, y, Z, rcount=rcount, ccount=ccount,
					#                        facecolors=colors, shade=False)
					# surf.set_facecolor((0,0,0,0))
					ax.set_xlabel('R')
					ax.set_ylabel('\u03A6')
					ax.set_zlabel(z_label)
					ax.set_xlim(ax.get_xlim()[::-1])
					# ax.gca().invert_xaxis()
					# textstr = '\n'.join((
		   #  		r'$\mathrm{n}=%.f$' % (n, ),
		   #  		r'$\mathrm{l}=%.f$' % (l, ),
		   #  		r'$\mathrm{ml}=%.f$' % (ml, ),
		   #  		r'$\mathrm{E}=%.f$' % (energy, )
		   #  		))
					# ax.text(0.88, 0.97, textstr, transform=ax.transAxes, fontsize=8,
		   #      	verticalalignment='top',horizontalalignment='left', bbox=dict(boxstyle='round', facecolor="wheat", alpha=0.5))
					file_name =  "/Users/Muhammad/Desktop/QP3D/Contour3D(BW)/PsiRPhi/" + "EigenFunc(nlm=%s).eps"%(str(n)+str(l)+str(ml))
					fig.savefig(file_name,format='eps', dpi=1000)
	for n in N:
		for l in range(n):
			for ml in range(-l,l+1):
				if ml < 3 and ml > -3:
					print(n)
					theta = np.linspace(0, 2*np.pi, 200)
					R = np.linspace(0,rad,80)
					fig = plt.figure()
					ax = plt.axes(projection='3d')
					ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
					ax.view_init(45, 115)
					x ,y = np.meshgrid(R, theta)
					z_label = ('\u03A8'+str(n)+str(l)+str(ml)).translate(sub)
					z ,energy= eig_fun(n,l,ml,x,y,0.01)
					ax.plot_wireframe(x, y, z,rstride=4, cstride=4, color='#414141')
					ax.set_xlabel('R')
					ax.set_ylabel('\u03B8')
					ax.set_zlabel(z_label)
					ax.set_xlim(ax.get_xlim()[::-1])
					# textstr = '\n'.join((
		   #  		r'$\mathrm{n}=%.f$' % (n, ),
		   #  		r'$\mathrm{l}=%.f$' % (l, ),
		   #  		r'$\mathrm{ml}=%.f$' % (ml, ),
		   #  		r'$\mathrm{E}=%.f$' % (energy, )
		   #  		))
					# ax.text(0.88, 0.97, textstr, transform=ax.transAxes, fontsize=8,
		   #      	verticalalignment='top',horizontalalignment='left', bbox=dict(boxstyle='round', facecolor="wheat", alpha=0.5))
					file_name =  "/Users/Muhammad/Desktop/QP3D/Contour3D(BW)/PsiRTheta/" + "EigenFunc(nlm=%s).eps"%(str(n)+str(l)+str(ml))
					fig.savefig(file_name,format='eps', dpi=1000)


q_3dplot()


