#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy             as np
import matplotlib.pyplot as plt
import random

from scipy.integrate import odeint
from lmfit           import minimize, Parameters, Parameter, report_fit
from scipy           import signal
from sklearn.metrics import mean_squared_error

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
		S  = 'Parameters for solveur:\n  N=' + str(self.N) + '\n  dt='+str(self.dt) + '\n  y0='+ str(self.y0) + '\n'
		S += self.modele.__str__()
		return S
		
	def setN(self, N):
		self.N = N
		# MAJ des autres paramètres en conséquence
		_, E0, I0, R10, R20 = self.y0
		S0 = self.N - E0 - I0 - R10 - R20
		self.y0 = S0, E0, I0, R10, R20
		# MAJ de N dans le modele
		self.modele.setN(self.N)
		# print('N-sum=', self.N-np.sum(self.y0))
		# input('pause')

	def set_paramf(self, f):
		self.modele.setf(f)

	def setParamInit(self, N, E0, I0, R10, R20):
		self.N = N
		S0 = N - E0 - I0 - R10 - R20
		self.y0 = S0, E0, I0, R10, R20
		# print('N-sum=', self.N-np.sum(self.y0))
		# input('pause')

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

	def plot_SEIR1R2(self, name, vectTime, plot, zs=(None, 0)):

		if len(plot)==0 or plot is None: pass

		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

		for p in plot:
			#print('self.solution[:, p]=', self.solution[:, p])
			ax.plot(vectTime, self.solution[:, p], color=self.modele.getColor(p), alpha=1.0, lw=2, label=self.modele.getString(p))

		# seconde estimation de R2
		# R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]
		# ax.plot(vectTime, R2bis, color=self.modele.getColor(indice), alpha=1.0, lw=2, label='$R^2(t)=N-\sum{SEIR^1}$')

		# les données observées
		if zs[0] != None:
			delta = zs[1]
			ax.plot(vectTime[delta:delta+len(zs[0])], zs[0], color=self.modele.getColor(3), alpha=1.0, lw=2, label='Observations ($R^1(n)$)', marker='x', ls='')

		ax.set_xlabel('Time (days)')
		ax.set_ylabel('Pop. size (' + str(int(self.N)) + ')')
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


	def paramOptimization(self, ts0, zs, vectTime):

		# set parameters including bounds; you can also fix parameters (use vary=False)
		S0, E0, I0, R10, R20 = self.y0 
		_, a0, b0, c0, f0    = self.modele.getParam()
	 
		params = Parameters()
		params.add('N',   value=self.N, vary=False)
		params.add('S0',  value=S0,     vary=False) #, min=1E7,    max=6.5E7)
		params.add('E0',  value=E0,     vary=False) #, min=0.,     max=200)
		params.add('I0',  value=I0,     vary=False) #, min=0.,     max=200)
		params.add('R10', value=R10,    vary=False) #, min=0.,     max=200)
		params.add('R20', value=R20,    vary=False) #, min=0.,     max=200)
		params.add('a',   value=a0,     vary=True, min=0.4*a0, max=2*a0) #min=0.01*a0, max=1.7*a0) #min=0.001, max=0.999)
		params.add('b',   value=b0,     vary=True, min=0.4*b0, max=2*b0) #min=0.01*b0, max=3.5*b0) #min=0.001, max=0.999)
		params.add('c',   value=c0,     vary=True, min=0.5*c0, max=2.3*c0) #min=0.03*c0, max=4.*c0)  #min=0.001, max=0.999)
		params.add('f',   value=f0,     vary=False) #, min=0.8*f0, max=1.4*f0)
		params.add('ts',  value=ts0,    vary=True,  min=0, max=len(vectTime)-len(zs)-1)
		
		# fit model
		result = minimize(residual, params, args=(vectTime, zs, self), method='least_squares')
		if self.verbose>0:
			result.params.pretty_print()

		# with the optimized values
		self.setParamInit(result.params['N'].value, result.params['E0'].value, result.params['I0'].value, result.params['R10'].value, result.params['R20'].value)
		self.modele.setParam(result.params['N'].value, result.params['a'].value, result.params['b'].value, result.params['c'].value, result.params['f'].value)
		if self.verbose>0:
			print('Solver''s state after optimization=', self)

		return result.params['ts'].value


def residual(paras, t, data, solveur):
	"""
	compute the residual between actual data and fitted data
	"""

	# print('current value', paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value, ' -- ', paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
	# print('current ts', int(paras['ts'].value))
	
	solveur.setParamInit(paras['N'].value, paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value)
	solveur.modele.setParam(paras['N'].value, paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
	model = solveur.solve_SEIR1R2_sol1(t)
	# print('model=', model)
	# input('apuse')

	# Only R1 to calculate the residual
	eqm=[]
	# print(paras['ts'].min, paras['ts'].max)
	for t in range(paras['ts'].min, paras['ts'].max):
		eqm.append(mean_squared_error(model[t:t+len(data), 3], data))
	delay = np.argmin(eqm)
	ts    = paras['ts'].min+delay
	paras['ts'].set(ts)
	# print('delay=', delay)
	# print('next ts', int(paras['ts'].value))

	# Only R1 to calculate the residual, only on the window's size of the data
	result = (model[ts:ts+len(data), 3] - data).ravel()
	# paras.pretty_print()
	#print('sum=', np.sum(np.abs(result)))

	return result



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

	# MAJ des parametres
	solveur.setParamInit(N, N*0.01, N*0.01, 23000, 0) # 1% de la population infectée
	solveur.modele.setParam(N=N, a=0.062, b=0.07692, c=0.1856, f=0.064)

	# integration time grid
	simulLenght = 80
	vectTime    = np.linspace(0, simulLenght-1, simulLenght)

	# Solveur to get the theoretical behavior
	########################################################################
	solveur.solve_SEIR1R2_sol1(vectTime)
	if plot==True:
		listePlot=[1,2,3]
		solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model.png', vectTime, plot=listePlot)
	
	# call main function
	# solveur.solve_SEIR1R2_sol2(simulLenght)
	# if plot==True:
	#	solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model_bis.png', vectTime)

	# Simulateur to get the theoretical behavior
	########################################################################
	# solveur.simul_SEIR1R2(simulLenght)
	# istePlot=[1,2,3]
	# if plot==True:
		# listePlot=[1,2,3]
	# 	solveur.plot_SEIR1R2(prefixFig+'simSEIR1R2model_pop.png', vectTime, plot=listePlot)

