from CoupledDiff import *

a_rec = 1/1100
a_eq = 1/24000
h0 = 67.74

sigma = 0.03

OMEGA_M = 0.3089
OMEGA_C = 0.2589
OMEGA_B = 0.0486

omega_m = 0.13
omega_b = 0.02230
omega_c = 0.1188

alpha = m.sqrt(a_rec/a_eq)
# tr = m.sqrt(4*a_rec/(OMEGA_M*(h0)**2))
tr = 269
t0 = 3.266 * 3 * (10**5) / 70
x_rec = (m.sqrt(1+alpha**2)-1)/alpha
x_s = 0.6*(omega_m**(1/4))*(omega_b**(-1/2))*(a_rec**(3/4))
k_d = 1/m.sqrt((2*x_s**2)+((sigma**2)*(x_rec**2)))
K_D = k_d/tr

x0 = 0.001
etha0 = (2*alpha)*(alpha*x0+1)/(((alpha*x0)**2)+2*alpha*x0)

tab_var, step = Variable_cal(0, 300, 1)

f_var = [i[1] for i in tab_var]
s_var = [j[2]*(3**-1/2) for j in tab_var]
kmpc = [k[0]for k in tab_var]

plt.plot(kmpc, f_var, label='First Component')
plt.plot(kmpc, s_var, label='Second Component')
plt.legend()
plt.show()
