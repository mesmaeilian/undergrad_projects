import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib import rc
import re


# phi1 = 0
# # phi2 = np.arange(0,np.pi/6,0.1)
# phi2 = 0
# # phi2 = [0.001,0.001]
# # w = [2,4,6]
# # w = [0,2,4,6]
# # w = 2
# # cmap = plt.cm.viridis
# # line_colors = cmap(np.linspace(0,1,len(phi2)))
# fig = plt.figure()
# for i, j in enumerate(phi2):
# 	t = np.arange(0,100,0.01)
# 	x = np.sin(t + phi1)
# 	y = np.sin(w*t + j)

# 	plt.plot(x,y, c=line_colors[i])

# # fig.savefig('/Users/Muhammad/Desktop/d.eps',format='eps', dpi=1000)
# plt.show()
# X1 = [np.sin(t) for t in np.arange(0,100,0.01)]
# X2 = []
# w = np.arange(2,2.6,0.2)
# # line_colors = cmap(np.linspace(0,1,len(w)))
# # fig = plt.figure()
# for q, p in enumerate(w):
# 		t = np.arange(0,100,0.01)
# 		x = np.sin(t + phi1)
# 		y = np.sin(p*t + phi2)
# 		# plt.plot(x,y, c=line_colors[q])
# 		X2.append(y)

# fig.savefig('/Users/Muhammad/Desktop/oddw.eps',format='eps', dpi=1000)

# fig, ax = plt.subplots()
# ax.set_xlim(np.min(X1), np.max(X1))
# ax.set_ylim(np.min(X2), np.max(X2))
# line_seg = LineCollection([np.column_stack((X1,x)) for x in X2], linewidths=(1), linestyles='solid')
# line_seg.set_array(w)
# line_seg.cmap = plt.cm.RdGy
# ax.add_collection(line_seg)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.text(1.04, 0.565,r'$\omega$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
# cbaxes = fig.add_axes([0.92, 0.58, 0.02, 0.3]) 
# fig.colorbar(line_seg, cax=cbaxes)

# fig.savefig('/Users/Muhammad/Desktop/b.eps',format='eps', dpi=1000)



phi1 = 0
# phi2 = 0

# w = 6
phi2 = np.arange(0,2,0.5)
# w = np.arange(2,3,0.2)
w = (3 + 8*np.pi)/(3*np.pi)
t = np.arange(0,20,0.01)

for i, p in enumerate(phi2) :
	fig = plt.figure()
	x = np.sin(t + phi1)
	y = np.sin(w*t + p)
	plt.plot(x,y,c='#D08633')
	f_name = '/Users/Muhammad/Desktop/' + 'd' + str(i) + '.eps'
	fig.savefig(f_name,format='eps',dpi=500)

















