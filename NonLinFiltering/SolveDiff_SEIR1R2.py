#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy             as np
import matplotlib.pyplot as plt
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

	def solve_SEIR1R2_sol1(self, t):

		# Integrate the SIR equations over the time grid t.
		self.solution = odeint(self.modele.deriv, self.y0, t, args=self.modele.getParam())
		return self.solution

	def plot_SEIR1R2(self, name, vectTime, zs=None):

		# set color cycle identical to panda's default color cycle
		# ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

		R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]

		#ax.plot(vectTime, self.solution[:, 0], color='black',   alpha=1.0, lw=2, label='$S(t)$')
		ax.plot(vectTime, self.solution[:, 1], color='#1f77b4', alpha=1.0, lw=2, label='$E(t)$')
		ax.plot(vectTime, self.solution[:, 2], color='#ff7f0e', alpha=1.0, lw=2, label='$I(t)$')
		ax.plot(vectTime, self.solution[:, 3], color='#2ca02c', alpha=1.0, lw=2, label='$R^1(t)$')
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

	# call main function
	solveur.solve_SEIR1R2_sol1(vectTime)
	if plot==True:
		solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model.png', vectTime)
	
	# call main function
	# solveur.solve_SEIR1R2_sol2(simulLenght)
	# if plot==True:
	#	solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model_bis.png', vectTime)
