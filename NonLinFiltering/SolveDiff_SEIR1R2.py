#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy             as np
import matplotlib.pyplot as plt
import random
from scipy.integrate import odeint

from SEIR_Bacaer     import SEIR1R2_Bacaer

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)


class SolveDiff_SEIR1R2:

	def __init__(self, N, dt=1, verbose=1):

		self.verbose = verbose
		
		# Total population, N.
		self.N = N
		# Initial number of infected and recovered individuals, I0 and R0.
		E0, I0, R10, R20 = 0, 13, 0, 0
		# Everyone else, S0, is susceptible to infection initially.
		S0 = self.N - E0 - I0 - R10 - R20
		# Initial conditions vector
		self.y0 = S0, E0, I0, R10, R20

		 # Modele d'eq. diff non linéaires
		self.dt = dt # 1 day
		self.modele = SEIR1R2_Bacaer(self.N, dt=dt)
		
	def __str__(self):
		S  = 'Parameters for solveur:\n  N=' + str(self.N) + '\n  y0='+ str(self.y0) + '\n'
		S += self.modele.__str__()
		return S
		
	def setParamInit(self, E0, I0, R10, R20):
		S0 = self.N - E0 - I0 - R10 - R20
		# Initial conditions vector
		self.y0 = S0, E0, I0, R10, R20

	def getParamInit(self):
		return self.y0

	def solve_SEIR1R2_sol2(self, simulLenght):
		self.solution = []

		# SOLUTION 2 - Integrate the SEIR1R2 equations over the time grid t, one by one.
		self.solution = np.zeros(shape=(simulLenght, 5))
		ret = odeint(self.modele.deriv, self.y0, [0, 0], args=self.modele.getParam())
		self.solution[0, :] = ret.T[:, -1]
		for i in range(simulLenght-1):
			ret = odeint(self.modele.deriv, self.solution[i], [0, 1], args=self.modele.getParam())
			self.solution[i+1, :] = ret.T[:, -1]

		return self.solution

	def solve_SEIR1R2_sol1(self, vectTime):

		# Integrate the SIR equations over the time grid t.
		self.solution = odeint(self.modele.deriv, self.y0, vectTime, args=self.modele.getParam())
		# print(type(self.solution))
		# print('shape=', np.shape(self.solution))
		# input('attente')
		return self.solution


	def simul_SEIR1R2(self, simulLenght):
		self.solution = np.zeros(shape=(simulLenght, 5), dtype=int) # 5 parameters: S, E, I, R1, R2

		# initalisation
		self.solution[0, :] = self.y0
		# print('self.solution[0, :]=', self.solution[0, :])
		# input('pause')


		print('\n')
		for t in range(1, simulLenght):
			print('\rTime ', t, ' over ', simulLenght, '      ', end='', flush = True)

			# compartiment S
			tab = np.random.random_sample(size=(self.solution[t-1, 0],))
			mask = (tab <= self.modele.a).astype(np.int)
			self.solution[t, 1] =  np.sum(mask)
			self.solution[t, 0] =  self.solution[t-1, 0] - self.solution[t, 1]
			# print('self.solution[t, 1]=', self.solution[t, 1])
			# print('np.sum(mask)=', np.sum(mask))
			# print(self.solution[t-1, 0] - self.solution[t, 1])
			# input('apuse')

			tab = np.random.random_sample(size=(self.solution[t-1, 1],))
			mask = (tab <= self.modele.b).astype(np.int)
			self.solution[t, 2]  =  np.sum(mask)
			self.solution[t, 1] +=  self.solution[t-1, 1] - self.solution[t, 2]

			tab = np.random.random_sample(size=(self.solution[t-1, 2],))
			mask = (tab <= self.modele.c).astype(np.int)
			temp = np.sum(mask)
			self.solution[t, 3] =  int(temp * self.modele.f)
			self.solution[t, 4] =  np.sum(mask)-self.solution[t, 3]
			self.solution[t, 2] +=  self.solution[t-1, 2] - temp
			
			# for individu in range(self.solution[t-1, 0]):
			# 	if random.random() <= self.modele.a: # passage de S à E
			# 		self.solution[t, 1] += 1
			# 	else:
			# 		self.solution[t, 0] += 1

			# # compartiment E
			# for individu in range(self.solution[t-1, 1]):
			# 	if random.random() <= self.modele.b: # passage de E à I
			# 		self.solution[t, 2] += 1
			# 	else:
			# 		self.solution[t, 1] += 1

			# # compartiment I
			# for individu in range(self.solution[t-1, 2]):
			# 	if random.random() <= self.modele.c: # passage de I vers R
			# 		if random.random() <= self.modele.f:
			# 			self.solution[t, 3] += 1
			# 		else:
			# 			self.solution[t, 4] += 1
			# 	else:
			# 		self.solution[t, 2] += 1


			# compartiments R1 e R2 : on ajoute le passé (ces compartiments cumulent)
			self.solution[t, 3] += self.solution[t-1, 3]
			self.solution[t, 4] += self.solution[t-1, 4]

			# print('self.solution[t, :]=', self.solution[t, :])
			# print('sum = ', np.sum(self.solution[t, :]))
			# input('pause')


		return self.solution


	def plot_SEIR1R2(self, name, vectTime, withR2=True, zs=None):

		# set color cycle identical to panda's default color cycle
		# ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

		R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]

		#ax.plot(vectTime, self.solution[:, 0], color='black',   alpha=1.0, lw=2, label='$S(t)$')
		ax.plot(vectTime, self.solution[:, 1], color='#1f77b4', alpha=1.0, lw=2, label='$E(t)$')
		ax.plot(vectTime, self.solution[:, 2], color='#ff7f0e', alpha=1.0, lw=2, label='$I(t)$')
		ax.plot(vectTime, self.solution[:, 3], color='#2ca02c', alpha=1.0, lw=2, label='$R^1(t)$')
		if withR2==True:
			ax.plot(vectTime, self.solution[:, 4], color='#d62728', alpha=1.0, lw=2, label='$R^2(t)$')
			ax.plot(vectTime, R2bis,               color='#9467bd', alpha=1.0, lw=2, label='$R^2(t)=N-\sum{SEIR^1}$')
		if zs != None:
			ax.plot(vectTime, zs, color='#8c564b', alpha=1.0, lw=2, label='Observations (R1)', marker='x', ls='')

		ax.set_xlabel('Time (days)')
		ax.set_ylabel('Number (' + str(self.N) + ')')
		#ax.set_ylim(0, self.N*1.01)
		ax.yaxis.set_tick_params(length=0)
		ax.xaxis.set_tick_params(length=0)
		ax.grid(b=True, which='major', c='w', lw=1, ls='-')
		legend = ax.legend()
		legend.get_frame().set_alpha(0.5)
		for spine in ('top', 'right', 'bottom', 'left'):
			ax.spines[spine].set_visible(False)
		plt.savefig(name, dpi=dpi)
		plt.close()



if __name__ == '__main__':

	prefixFig = './figures/'

	# Solveur eq. diff.
	verbose = 1
	plot    = True
	N       = 65.E6
	dt      = 1
	solveur = SolveDiff_SEIR1R2(N, dt, verbose)
	if verbose>0:
		print(solveur)

	# integration time grid
	simulLenght = 600
	vectTime    = np.linspace(0, simulLenght-1, simulLenght)

	# Solveur to get the theoretical behavior
	########################################################################
	solveur.solve_SEIR1R2_sol1(vectTime)
	if plot==True:
		solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model.png', vectTime, withR2= False)
	
	# call main function
	# solveur.solve_SEIR1R2_sol2(simulLenght)
	# if plot==True:
	#	solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model_bis.png', vectTime)

	# Simulateur to get the theoretical behavior
	########################################################################
	# solveur.simul_SEIR1R2(simulLenght)
	# if plot==True:
	# 	solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model_ter.png', vectTime, withR2= False)



