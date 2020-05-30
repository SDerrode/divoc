#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator

from scipy.integrate import odeint
from lmfit           import minimize, Parameters, Parameter, report_fit
from scipy           import signal
from sklearn.metrics import mean_squared_error

from SEIR_Bacaer     import SEIR1R2_Bacaer
from common          import addDaystoStrDate

dpi     = 150    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)


class SolveEDO_SEIR1R2:

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

		self.TS = -1 # time shift with respect to some data

		 # Modele d'eq. diff non linéaires
		self.dt = dt # 1 day
		self.modele = SEIR1R2_Bacaer(self.N, dt=dt)
		
	def __str__(self):
		S  = 'Solveur:\n  N=' + str(self.N) + '\n  dt='+str(self.dt) + '\n  y0='+ str(self.y0) + '\n  TS='+ str(self.TS) + '\n'
		S = self.modele.__str__()
		return S

	def getTextParam(self, startDate=None):
		S  = self.modele.getTextParam()
		S += '\nSolveur init:\n'# + '\n  N=' + str(self.N) + '\n  dt='+str(self.dt) \
		S += r'  $S_0=' + str(self.y0[0]) + '$\n'
		S += r'  $E_0=' + str(self.y0[1]) + r', I_0='+str(self.y0[2]) + '$\n'
		S += r'  $R^1_0=' + str(self.y0[3]) + r', R^2_0=' + str(self.y0[4])+ '$'
		if self.TS != -1:
			if startDate==None:
				S += '\n  Time Shift=$'+str(self.TS)+'$'
			else:
				dateI0 = addDaystoStrDate(startDate, -self.TS)
				S += '\n  Start date:'+dateI0
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

	def getParamInit(self):
		return self.y0

	def compute_tsfromEQM(self, data, T):
		
		dataLength = len(data)
		eqm = np.zeros(shape=(T-dataLength))
		for t in range(T-dataLength):
			eqm[t] = mean_squared_error(self.solution[t:t+dataLength, 3], data)
		self.TS = np.argmin(eqm)
		return self.TS

	def solve_SEIR1R2_withSwitch(self, T, timeswitch=-1):
		'''
		Integrate the SEIR1R2 equations over the time grid t, one by one; to dela with parameter switching
		'''
		#print('timeswitch=', timeswitch)

		self.solution = np.zeros(shape=(T, 5))
		ret = odeint(self.modele.deriv, self.y0, [0, 0], args=self.modele.getParam())
		self.solution[0, :] = ret.T[:, -1]
		for i in range(1,T):
			if i == timeswitch:
				self.modele.seta(0.)
				self.modele.setb(0.)
			ret = odeint(self.modele.deriv, self.solution[i-1], [0, 1], args=self.modele.getParam())
			self.solution[i, :] = ret.T[:, -1]

		return self.solution

	def solve_SEIR1R2(self, time):
		# Integrate the SIR equations over the time grid t.
		self.solution = odeint(self.modele.deriv, self.y0, time, args=self.modele.getParam())
		# print(type(self.solution))
		# print('shape=', np.shape(self.solution))
		# input('attente')
		return self.solution

	def paramOptimization(self, data, time):

		# set parameters including bounds; you can also fix parameters (use vary=False)
		S0, E0, I0, R10, R20 = self.y0 
		_, a0, b0, c0, f0    = self.modele.getParam()
	 
		params = Parameters()
		params.add('N',   value=self.N,  vary=False)
		params.add('E0',  value=E0,      vary=False) #True, min=0 , max=2000)
		params.add('I0',  value=I0,      vary=False) #True, min=0 , max=2000)
		params.add('R10', value=R10,     vary=False) #True, min=0 , max=2000)
		params.add('R20', value=R20,     vary=False) #True, min=0 , max=2000)
		params.add('ts',  value=self.TS, vary=True, min=0,     max=len(time)-len(data)-2)
		params.add('a',   value=a0,      vary=True, min=0.0001, max=0.85) 
		params.add('b',   value=b0,      vary=True, min=0.100, max=0.6) 
		params.add('c',   value=c0,      vary=True, min=0.010, max=0.2) 
		params.add('f',   value=f0,      vary=True, min=0.001, max=0.4)
		# params.add('a',   value=a0,      vary=True, min=0.3*a0, max=1./0.3*a0)
		# params.add('b',   value=b0,      vary=True, min=0.5*b0, max=1./0.5*b0)
		# params.add('c',   value=c0,      vary=True, min=0.5*c0, max=1./0.5*c0)
		# params.add('f',   value=f0,      vary=True, min=0.4*f0, max=1./0.4*f0)
		
		# fit model
		result = minimize(residual, params, args=(time, data, self), method='powell') #powell, least_squares
		if self.verbose>0:
			result.params.pretty_print()

		# with the optimized values
		self.setParamInit(result.params['N'].value, result.params['E0'].value, result.params['I0'].value, result.params['R10'].value, result.params['R20'].value)
		self.modele.setParam(result.params['N'].value, result.params['a'].value, result.params['b'].value, result.params['c'].value, result.params['f'].value)
		
		self.TS = result.params['ts'].value


	# def simul_SEIR1R2(self, T):
	# 	self.solution = np.zeros(shape=(T, 5), dtype=int) # 5 parameters: S, E, I, R1, R2

	# 	# initalisation
	# 	self.solution[0, :] = self.y0
	# 	# print('self.solution[0, :]=', self.solution[0, :])
	# 	# input('pause')

	# 	print('\n')
	# 	for t in range(1, T):
	# 		print('\rTime ', t, ' over ', T, '      ', end='', flush = True)

	# 		# compartiment S
	# 		tab = np.random.random_sample(size=(self.solution[t-1, 0],))
	# 		mask = (tab <= self.modele.a).astype(np.int)
	# 		self.solution[t, 1] =  np.sum(mask)
	# 		self.solution[t, 0] =  self.solution[t-1, 0] - self.solution[t, 1]
	# 		# print('self.solution[t, 1]=', self.solution[t, 1])
	# 		# print('np.sum(mask)=', np.sum(mask))
	# 		# print(self.solution[t-1, 0] - self.solution[t, 1])
	# 		# input('apuse')

	# 		tab = np.random.random_sample(size=(self.solution[t-1, 1],))
	# 		mask = (tab <= self.modele.b).astype(np.int)
	# 		self.solution[t, 2]  =  np.sum(mask)
	# 		self.solution[t, 1] +=  self.solution[t-1, 1] - self.solution[t, 2]

	# 		tab = np.random.random_sample(size=(self.solution[t-1, 2],))
	# 		mask = (tab <= self.modele.c).astype(np.int)
	# 		temp = np.sum(mask)
	# 		self.solution[t, 3] =  int(temp * self.modele.f)
	# 		self.solution[t, 4] =  np.sum(mask)-self.solution[t, 3]
	# 		self.solution[t, 2] +=  self.solution[t-1, 2] - temp
			
	# 		# for individu in range(self.solution[t-1, 0]):
	# 		# 	if random.random() <= self.modele.a: # passage de S à E
	# 		# 		self.solution[t, 1] += 1
	# 		# 	else:
	# 		# 		self.solution[t, 0] += 1

	# 		# # compartiment E
	# 		# for individu in range(self.solution[t-1, 1]):
	# 		# 	if random.random() <= self.modele.b: # passage de E à I
	# 		# 		self.solution[t, 2] += 1
	# 		# 	else:
	# 		# 		self.solution[t, 1] += 1

	# 		# # compartiment I
	# 		# for individu in range(self.solution[t-1, 2]):
	# 		# 	if random.random() <= self.modele.c: # passage de I vers R
	# 		# 		if random.random() <= self.modele.f:
	# 		# 			self.solution[t, 3] += 1
	# 		# 		else:
	# 		# 			self.solution[t, 4] += 1
	# 		# 	else:
	# 		# 		self.solution[t, 2] += 1


	# 		# compartiments R1 e R2 : on ajoute le passé (ces compartiments cumulent)
	# 		self.solution[t, 3] += self.solution[t-1, 3]
	# 		self.solution[t, 4] += self.solution[t-1, 4]

	# 		# print('self.solution[t, :]=', self.solution[t, :])
	# 		# print('sum = ', np.sum(self.solution[t, :]))
	# 		# input('pause')

	# 	return self.solution

	def plot_SEIR1R2(self, name, title, sliceedo, slicedata, plot, data='', text=''):

		if len(plot)==0 or plot is None: pass

		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

		a = slicedata.start+sliceedo.stop-sliceedo.start-1
		timeedo = np.linspace(slicedata.start, a, a-slicedata.start+1)
		for p in plot:
			ax.plot(timeedo, self.solution[sliceedo.start:min(sliceedo.stop, np.shape(self.solution)[0]), p], color=self.modele.getColor(p), alpha=1.0, lw=2, label=self.modele.getString(p))

		# seconde estimation de R2
		# R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]
		# ax.plot(timeedo, R2bis, color=self.modele.getColor(indice), alpha=1.0, lw=2, label='$R^2(t)=N-\sum{SEIR^1}$')

		# les données observées
		time = np.linspace(slicedata.start, slicedata.stop-1, slicedata.stop-slicedata.start)
		if len(data) != 0:
			ax.plot(time, data[slicedata], color=self.modele.getColor(3), alpha=1.0, lw=2, label='Data $R^1$', marker='x', ls='')

		ax.set_xlabel('Time (days)')
		ax.set_ylabel('Pop. size (' + str(int(self.N)) + ')')
		ax.yaxis.set_tick_params(length=0)
		ax.xaxis.set_tick_params(length=0)
		ax.grid(b=True, which='major', c='w', lw=1, ls='-')
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		
		legend = ax.legend()
		legend.get_frame().set_alpha(0.5)
		for spine in ('top', 'right', 'bottom', 'left'):
			ax.spines[spine].set_visible(False)

		# ajout d'un text d'annotation
		if text != '':
			bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
			ax.annotate(text, xy=(time[0], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=8, bbox=bbox, ha="left", va="center") 
		plt.title(title)
		plt.savefig(name, dpi=dpi)
		plt.close()


	def plot_dR1(self, name, title, sliceedoderiv, slicedataderiv, deriv_data='', text=''):

		time = np.linspace(slicedataderiv.start, slicedataderiv.stop-1, slicedataderiv.stop-slicedataderiv.start)

		deriv_sol3 = (self.solution[1:, 3] - self.solution[:-1, 3]) / self.dt

		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
		
		ax.plot(time, deriv_sol3[sliceedoderiv.start:min(sliceedoderiv.stop, np.shape(deriv_sol3)[0])], color=self.modele.getColor(3), alpha=1.0, lw=2, label=r'$\frac{\partial R^1(t)}{\partial t}$')

		# seconde estimation de R2
		# R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]
		# ax.plot(timeedo, R2bis, color=self.modele.getColor(indice), alpha=1.0, lw=2, label='$R^2(t)=N-\sum{SEIR^1}$')

		# les données observées
		if len(deriv_data) != 0:
			ax.plot(time, deriv_data, color=self.modele.getColor(3), alpha=1.0, lw=0.5, label=r'$\frac{\partial R^1(n)}{\partial n}$', marker='x')

		ax.set_xlabel('Time (days)')
		ax.yaxis.set_tick_params(length=0)
		ax.xaxis.set_tick_params(length=0)
		ax.grid(b=True, which='major', c='w', lw=1, ls='-')
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		
		legend = ax.legend()
		legend.get_frame().set_alpha(0.5)
		for spine in ('top', 'right', 'bottom', 'left'):
			ax.spines[spine].set_visible(False)

		# ajout d'un text d'annotation
		if text != '':
			bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
			ax.annotate(text, xy=(time[0], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=8, bbox=bbox, ha="left", va="center") 
		plt.title(title)
		plt.savefig(name, dpi=dpi)
		plt.close()


def residual(paras, t, data, solveur):
	"""
	compute the residual between actual data and fitted data
	"""

	# print('current value', paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value, ' -- ', paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
	# print('current ts', int(paras['ts'].value))
	
	solveur.setParamInit(paras['N'].value, paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value)
	solveur.modele.setParam(paras['N'].value, paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
	model = solveur.solve_SEIR1R2(t)

	# Only R1 to calculate the residual
	eqm=[]
	for t in range(paras['ts'].min, paras['ts'].max):
		eqm.append(mean_squared_error(model[t:t+len(data), 3], data))
	delay = np.argmin(eqm)
	ts    = paras['ts'].min+delay
	paras['ts'].set(ts)
	# print('delay=', delay)
	# print('next ts', int(paras['ts'].value))

	# print('current value', paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value, ' -- ', paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
	# print('current ts', int(paras['ts'].value))

	# print(data)
	# print(model[ts:ts+len(data), 3])
	#print(model[:, 3])
	#input('attente')

	# Dessin des courbes théoriques
	# fig = plt.figure(facecolor='w',figsize=figsize)
	# ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
	# ax.plot(model[ts:ts+len(data), 3], color='red',   label='model')
	# ax.plot(data,                 color='green', label='data', marker='x', ls='')
	# plt.legend()
	# plt.show()

	# Only R1 to calculate the residual, only on the window's size of the data
	result = model[ts:ts+len(data), 3] - data
	# print('result=', result)
	#print('sum=', np.sum(np.multiply(result, result)))
	# input('pause')
	# paras.pretty_print()
	#print('sum=', np.sum(np.abs(result)))

	return result



if __name__ == '__main__':

	# Constantes
	import os
	repertoire = './figures/simulSEIR1R2'
	if not os.path.exists(repertoire):
		os.makedirs(repertoire)
	prefFig = repertoire + '/'

	# Solveur eq. diff.
	verbose = 1
	plot    = True
	N       = 66987244 # Population de la France
	dt      = 1
	solveur = SolveEDO_SEIR1R2(N, dt, verbose)
	if verbose>0:
		print(solveur)

	# MAJ des parametres
	
	E0, I0, R10, R20 = 0, 1, 0, 0
	a, b, c, f       = 0.11, 0.24, 0.060, 0.20

	solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20)
	solveur.modele.setParam(N=N,  a=a,   b=b,    c=c,     f=f)

	# integration time grid
	T    = 750
	time = np.linspace(0, T-1, T)

	# Solveur to get the theoretical behavior
	########################################################################
	result = solveur.solve_SEIR1R2(time)
	if plot==True:
		
		sliceedo = slice(0, 0+T)
		
		listePlot = [3]
		filename  = prefFig + 'SEIR1R2model_' + ''.join(map(str, listePlot)) + '.png'
		solveur.plot_SEIR1R2(filename, '', sliceedo, sliceedo, plot=listePlot, data='', text=solveur.getTextParam())
		listePlot = [1,2,3]
		filename  = prefFig + 'SEIR1R2model_' + ''.join(map(str, listePlot)) + '.png'
		solveur.plot_SEIR1R2(filename, '', sliceedo, sliceedo, plot=listePlot, data='', text=solveur.getTextParam())
		listePlot = [0,1,2,3,4]
		filename  = prefFig + 'SEIR1R2model_' + ''.join(map(str, listePlot)) + '.png'
		solveur.plot_SEIR1R2(filename, '', sliceedo, sliceedo, plot=listePlot, data='', text=solveur.getTextParam())


	# call main function
	# solveur.solve_SEIR1R2_withSwitch(T, xx)
	# if plot==True:
	#	solveur.plot_SEIR1R2(prefFig+'simSEIR1R2model_bis.png', time)

	# Simulateur to get the theoretical behavior
	########################################################################
	# solveur.simul_SEIR1R2(T)
	# istePlot=[1,2,3]
	# if plot==True:
		# listePlot=[1,2,3]
	# 	solveur.plot_SEIR1R2(prefFig+'simSEIR1R2model_pop.png', time, plot=listePlot)
