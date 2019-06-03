import numpy as np 
import matplotlib.pyplot as plt
from scipy.special import spherical_jn as sjn
from scipy.special import sph_harm as sph
from scipy.integrate import quad
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal



def prob(rad,n,l,ml,zero,up_limit,step):

	h = 1
	m = 1
	func = lambda r: (sjn(l,r*zero/rad)**2)*(r**2)
	norm_fact = quad(func, 0, rad)
	energy = (zero**2)*(h**2)/(2*m*(rad**2))
	k = 2 * m * energy / (h**2)
	probability = []
	R = np.arange(0,up_limit,step)
	theta = np.arange(0, 2*np.pi, 0.1)
	phi = np.arange(0, np.pi, 0.1)
	counter = 0
	for r in R:
		for i in theta:
			for j in phi:
				counter += 1
				spat = (((-1**l) * sph(ml,l,[i],[j]) * sph(-ml,l,[i],[j]))[0]).real
				if spat < 0 :
					spat *= -1
				# print(">>>>>>>>",spat,"r:",r,"theta:",i,"phi:",j)
				radi = (norm_fact[0]) * (sjn(l,r*k))
				psi_2 = (radi**2)*spat
				# print("r,theta,phi,psi2 >>>>>",r,i,j,psi_2)
				probability.append([r,i,j,psi_2])
				# print(probability)
	return probability

def prob_normalization(probability):

	probabilities = []
	for i in probability:
		probabilities.append(i[3])
		# print(i[3])

	tot_probability = sum(probabilities)
	# print(tot_probability)
	n_probability = []
	for j in probabilities:
		n_probability.append(j/tot_probability)
		print(j/tot_probability)


	# return n_probability
	return probabilities

def plt_opt_prob():
	R = 1
	up_limit = 1
	step = 0.02
	# N = range(1,6)
	N = [1]
	sub = str.maketrans("0123456789","₀₁₂₃₄₅₆₇₈₉")
	for n in N:
		for l in range(n):
			for ml in range(-l,l+1):
				zero = spbess_zero(n,l,4)
				probability = prob(R,n,l,ml,zero,up_limit,step)
				n_probability = prob_normalization(probability)
				print(len(probability),len(n_probability))
				# print("Normalized Probability >>>>", n_probability)
				x_coordinates = []
				y_coordinates = []
				z_coordinates = []
				counter = -1
				for i in n_probability:
					counter += 1
					# print(i)
					if i > 9*10e-7:
						x = R * np.sin(probability[counter][3]) * np.cos(probability[counter][2])
						y = R * np.sin(probability[counter][3]) * np.sin(probability[counter][2])
						z = R * np.cos(probability[counter][3])
						x_coordinates.append(x)
						y_coordinates.append(y)
						z_coordinates.append(z)

				# theta, phi = np.arange(0, 2 * np.pi, 0.01), np.arange(0, np.pi, 0.01)
				# THETA, PHI = np.meshgrid(theta, phi)
				# X = R * np.sin(PHI) * np.cos(THETA)
				# Y = R * np.sin(PHI) * np.sin(THETA)
				# Z = R * np.cos(PHI)
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				ax.scatter(x_coordinates, y_coordinates, z_coordinates, c='r', marker='.')
				ax.set_xlabel('X')
				ax.set_ylabel('Y')
				ax.set_zlabel('Z')
				plt.show()
				# file_name =  "/Users/mohammadsadegh/Desktop/QP(3d)/" + "ProbabilityDensity(nlm=%s).eps"%(str(n)+str(l)+str(ml))
				# fig.savefig(file_name,format='eps', dpi=1000)