import numpy as np 
import matplotlib.pyplot as plt
from scipy.special import spherical_jn as sjn
from scipy.special import sph_harm as sph
from scipy.integrate import quad
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D


h = 1.05457 * 10e-34
m = 1
rad = 1

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



def eig_fun(rad,n,l,ml,zero,up_limit,step):
	# h = 1.0545718 * 10e-34
	# m = 9.1093836 * 10e-31
	# rad = 5.3 * 10e-11
	h = 1
	m = 1
	func = lambda r: (sjn(l,r*zero/rad)**2)*(r**2)
	norm_fact = quad(func, 0, rad)
	energy = (zero**2)*(h**2)/(2*m*(rad**2))
	k = 2 * m * energy / (h**2)
	solutions = []
	R = np.arange(0,up_limit,step)
	# theta = np.arange(0, 2*np.pi, 0.01)
	# phi = np.arange(0, np.pi, 0.01)
	for r in R:
		sol = (norm_fact[0])*(sjn(l,r*k))*(sph(ml,l,[0,np.pi/6],[0,np.pi/3]))
		solutions.append(sol)
	return R, solutions, energy


def opt_plt(height_space, width_space, grid_stat, X_label, Y_label, T_X_pos, T_Y_pos, font_size, face_color):

	up_limit = 1
	step = 0.01
	N = range(1,6)
	sub = str.maketrans("0123456789","₀₁₂₃₄₅₆₇₈₉")
	count = 
	fig = plt.figure()
	for n in N:
		i = 0
		for l in range(n):
			for ml in range(-l,l+1):
				if ml < 3 and ml > -3:
					i += 1 
					zero = spbess_zero(l,n,4)
					r, solutions, energy = eig_fun(1,n,l,ml,zero,up_limit,step)
					print(energy)
					fig.subplots_adjust(hspace=height_space,wspace=width_space)
					ax = fig.add_subplot(n,n,i)
					# ax.set_xlabel(X_label)
					# scale_y = 1
					# ticks_y = ticker.FuncFormatter(lambda x,pos: '{0:g}'.format(solutions/scale_y))
					# ax.yaxis.set_major_formatter(ticks_y)
					# y_label = (Y_label+str(n)+str(l)+str(ml)).translate(sub)
					textstr = '\n'.join((
		    		r'$\mathrm{n}=%.f$' % (n, ),
		    		r'$\mathrm{l}=%.f$' % (l, ),
		    		r'$\mathrm{E}=%.f$' % (energy, )
		    		))
					ax.text(T_X_pos, T_Y_pos, textstr, transform=ax.transAxes, fontsize=font_size,
		        	verticalalignment='top', bbox=dict(boxstyle='round', facecolor=face_color, alpha=0.5))
					ax.grid(grid_stat)
					ax.plot(r,solutions,c='#1f77b4')
	file_name =  "/Users/Muhammad/Desktop/radial.eps"
	fig.savefig(file_name,format='eps', dpi=1000)

# , top=0.95, bottom=0.1, left=0.15, right=0.95
opt_plt(0.4, 0.6, False, "r", "\u03C8", 0.88, 0.97, 8, "wheat")


# plt_opt_prob()


