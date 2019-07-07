import numpy as np
import matplotlib.pyplot as plt


b = 7

L = 2
m = 1
k = 0.5
E = -1
e = 1/5

tet = np.arange(0,10*np.pi,0.01)

r = (((b**2)*(L**2))/(m*k))/(1 - e*np.cos(b*tet))

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='polar')
ax.plot(tet, r)
ax.set_yticklabels([])
ax.set_theta_zero_location('N')
f_name = '/Users/Muhammad/Desktop/' + 'baste1' + '.eps'
fig.savefig(f_name,format='eps',dpi=1000)