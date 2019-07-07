import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn

tau = [6400,14400]
l = 100

for j in tau:

	k = np.arange(0,1.2,0.002)
	bes = [spherical_jn(l, j*i, derivative = True) for i in k]
	ll = 'tau = ' + str(j)
	plt.plot(k,bes,label = ll)

plt.legend()
plt.show()