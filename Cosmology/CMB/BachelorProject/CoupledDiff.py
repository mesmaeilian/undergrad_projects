import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math as m
from scipy.integrate import quad
from scipy.special import spherical_jn
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import axes3d
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib import rc
import re

cmap = plt.cm.PRGn

h = 0.6766

OMEGA_M = 0.3089
OMEGA_C = 0.2589
OMEGA_B = 0.0486

omega_m = OMEGA_M * (h**2)
omega_b = OMEGA_B * (h**2)
omega_c = OMEGA_C * (h**2) 

a_rec = 1/1100
a_eq_1 = 24000*omega_m
alpha = m.sqrt(a_rec*a_eq_1)

x_rec = (m.sqrt(1+(alpha**2))-1)/alpha

sigma = 0.03
tr = 284
t_zero = 4200

de = 1.476
x_s = 0.6*(omega_m**(1/4))*(omega_b**(-1/2))*(a_rec**(3/4))
k_d = 1/m.sqrt((2*x_s**2)+((sigma**2)*(x_rec**2)))
K_D = k_d/tr
t0 = te_zero*de



def VacFlucInit(k,x):

	etha = (2*alpha) * ((alpha*x)+1) / (((alpha*x)**2)+(2*alpha*x))
	phi = 1
	d_n = -2*phi
	d_c = 3*d_n/4
	v_n = -k*d_n/(4*etha)
	v_c = v_n

	return [phi,d_c,v_c,d_n,v_n]

def VacFluc(X,x,k):

	PHI = X[0]
	D_C = X[1]
	V_C = X[2]
	D_N = X[3]
	V_N = X[4]

	etha = (2*alpha) * ((alpha*x)+1) / (((alpha*x)**2)+(2*alpha*x))
	y = (2*alpha*x)+((alpha*x)**2)
	y_c = (y*OMEGA_C)/(OMEGA_M)
	y_b = (y*OMEGA_B)/(OMEGA_M)

	dPHIdx = (-etha * PHI) + ((3*(etha**2)/(2*k))*(((V_N*((4/3)+y-y_c))+(V_C*y_c))/(1+y)))
	dD_Cdx = (-k * V_C) + (3 * dPHIdx)
	dV_Cdx = (-etha * V_C) + (k * PHI)
	dD_Ndx = ((-4 * k * V_N)/3) + (4 * dPHIdx)
	dV_Ndx = ((1/(1+(3*y_b/4)))*((-3*y_b*etha*V_N/4)+(k*D_N/4)))+(k*PHI)

	return [dPHIdx, dD_Cdx, dV_Cdx, dD_Ndx, dV_Ndx]

def Variable_cal(k_min, k_max, step):
	x0 = 0.001
	x = np.linspace(x0, x_rec)
	var_x_rec = []

	for i in np.arange(k_min, k_max, step):
		X0 = VacFlucInit(i,x0)
		k = i
		sol = odeint(VacFluc, X0, x, args=(k,))
		var = [k,((sol[:,0][-1])+(sol[:,3][-1]/4)),(sol[:,4][-1])]
		print(var)
		var_x_rec.append(var)
	
	return var_x_rec, step, x
	
def Spectrum(tab_var, x, step):

	Powers = []
	L = np.arange(1, 1500, 25)
	for i in L:

		k = [j[0] for j in tab_var]
		a = [j[1] for j in tab_var]
		b = [j[2] for j in tab_var]
		dk = step

		func_var = list(map(lambda k,a,b:(2)*i*(i+1)*((a*spherical_jn(i,(t0-tr)*k/tr,derivative=False)+b*spherical_jn(i,(t0-tr)*k/tr,1))**2)*np.exp(-(k/(tr*K_D))**2)*(k**(-1.04)), k,a,b))
		integ = np.sum([j*dk for j in func_var])
		print(i,">>>>>>>>>>>>",integ)
		Powers.append(integ)

	return L, Powers

def Smooth_func(x_lst, y_lst, number):

	x_new = []
	y_new = []
	counter = 0

	while counter != len(x_lst) - 1 - int((number-1)/2) :

		x = [i for i in x_lst[counter:counter+number]]
		y = [j for j in y_lst[counter:counter+number]]
		x_mean = np.mean(x)
		y_mean = np.mean(y)
		x_new.append(x_mean)
		y_new.append(y_mean)
		counter += 1

	return x_new, y_new

def interp_func(l, powers):

	f = interp1d(l, powers, kind='cubic')
	
	return l, f


# t = np.arange(0.02,0.06,0.01)
# for i, j in enumerate(t):
# 	OMEGA_B = j
# 	OMEGA_M = 0.3089
# 	OMEGA_C = 0.2589
# 	coup_ans, step, x = Variable_cal(0.01,300,0.05)
# 	L, Powers = Spectrum(coup_ans, x, step)
# 	# fig = plt.figure()	
# 	line_colors = cmap(np.linspace(0.2,1,len(t)))
# 	plt.semilogx(L, Powers, c=line_colors[i])
# plt.show()

power1 = []
t = np.arange(0.04,0.14,0.01)
for i, j in enumerate(t):
	OMEGA_B = j
	coup_ans, step, x = Variable_cal(0.01,300,20)
	L, Powers = Spectrum(coup_ans, x, step)
	power1.append(Powers)



fig, ax = plt.subplots()
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
# fig.subplots_adjust(top=0.9, bottom=0.1,left=0.1,right=0.9)
ax.set_xlim(np.min(L), np.max(L))
ax.set_ylim(np.min(power1), np.max(power1))


ax.set_yticks([0, .2, .4, .6, .8, 1])
line_seg = LineCollection([np.column_stack((L,power)) for power in power1], linewidths=(1), linestyles='solid')
line_seg.set_array(t)
line_seg.cmap = plt.cm.RdGy
# line_seg.cmap = plt.cm.PRGn
ax.add_collection(line_seg)
ax.semilogx()
ax.set_xticks([10,100,1000,10000])
plt.xlabel(r'${l}$')
plt.ylabel(r'${l(l+1){C_l}^{TT}}$')
plt.text(0.066, 0.51,r'$\Omega_b$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
cbaxes = fig.add_axes([0.16, 0.55, 0.02, 0.3]) 
fig.colorbar(line_seg, cax=cbaxes)

# plt.sci(line_seg)
plt.show()
# file_name =  "/Users/Muhammad/Desktop/21.eps"
# fig.savefig(file_name,format='eps', dpi=1000)

