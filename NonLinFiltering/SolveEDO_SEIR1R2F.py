#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator

from scipy.integrate import odeint
from lmfit           import minimize, Parameters, Parameter, report_fit
from sklearn.metrics import mean_squared_error

from SolveEDO        import SolveEDO
from SEIR1R2F        import SEIR1R2F
from common          import addDaystoStrDate, getRepertoire

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)


class SolveEDO_SEIR1R2F(SolveEDO):

	def __init__(self, N, dt=1, verbose=1):

		# appel constructeur classe mère
		super().__init__(N, dt, verbose)

		# Initial number of infected and recovered individuals, I0 and R0.
		E0, I0, R10, R20, F0 = 0, 1, 0, 0, 0
		# Everyone else, S0, is susceptible to infection initially.
		S0 = self.N - E0 - I0 - R10 - R20 - F0
		# Initial conditions vector
		self.y0 = S0, E0, I0, R10, R20, F0
		self.nbparam = len(self.y0)

		self.indexdata=[3,5]

		# Modèle d'eq. diff non linéaires
		self.modele = SEIR1R2F(self.N, dt=self.dt)

	def getTextParam(self, startDate=None):
		S = super().getTextParam(startDate)
		S += '\n' + r'  $F_0=' + str(self.y0[5]) + '$'
		return S

	def getTextParamWeak(self, startDate=None):
		S = super().getTextParamWeak(startDate)
		return S

	def setN(self, N):
		self.N = N
		# MAJ des autres paramètres en conséquence
		_, E0, I0, R10, R20, F0 = self.y0
		S0 = self.N - E0 - I0 - R10 - R20 - F0
		self.y0 = S0, E0, I0, R10, R20, F0
		# MAJ de N dans le modele
		self.modele.setN(self.N)

	def setParamInit(self, N, E0, I0, R10, R20, F0):
		self.N = N
		S0 = N - E0 - I0 - R10 - R20 - F0
		self.y0 = S0, E0, I0, R10, R20, F0
		# MAJ de N dans le modele
		self.modele.setN(self.N)

	def paramOptimization(self, data, time, ts=None):

		# set parameters including bounds; you can also fix parameters (use vary=False)
		S0, E0, I0, R10, R20, F0     = self.y0 
		_, a0, b0, c0, f0, muI0, xi0 = self.modele.getParam()

		params = Parameters()
		params.add('N',   value=self.N,  vary=False)
		params.add('E0',  value=E0,      vary=False)
		params.add('I0',  value=I0,      vary=False)
		params.add('R10', value=R10,     vary=False)
		params.add('R20', value=R20,     vary=False)
		params.add('F0',  value=F0,      vary=False)
		
		params.add('a',   value=a0,      vary=True, min=0.0, max=0.99) 
		params.add('b',   value=b0,      vary=True, min=0.0, max=0.99) 
		params.add('c',   value=c0,      vary=True, min=0.0, max=0.99) 
		params.add('f',   value=f0,      vary=True, min=0.0, max=0.99)
		params.add('muI', value=muI0,    vary=True, min=0.00001, max=0.01)
		params.add('xi',  value=xi0,     vary=False)
		
		if ts != None:
			params.add('ts', value=ts, vary=False)
		else:
			params.add('ts', value=self.TS, vary=True, min=0, max=len(time)-np.shape(data)[0]-2)
		
		# fit model
		result = minimize(residual, params, args=(time, data, self, self.indexdata), method='powell') #powell, least_squares
		if self.verbose>0:
			result.params.pretty_print()
			# display fitted statistics
			report_fit(result)

		# Update params with the optimized values
		self.setParamInit(result.params['N'].value, result.params['E0'].value, result.params['I0'].value, result.params['R10'].value, result.params['R20'].value, result.params['F0'].value)
		self.modele.setParam(result.params['N'].value, result.params['a'].value, result.params['b'].value, result.params['c'].value, result.params['f'].value, result.params['muI'].value, result.params['xi'].value)
		self.TS = result.params['ts'].value


def residual(paras, t, data, solveur, indexdata):
	"""
	compute the residual between actual data and fitted data
	"""

	solveur.setParamInit(paras['N'].value, paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value, paras['F0'].value)
	solveur.modele.setParam(paras['N'].value, paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value, paras['muI'].value, paras['xi'].value)
	model = solveur.solveEDO(t)

	# Only R1 and F to calculate the residual
	ts = paras['ts'].value
	if paras['ts'].vary == True:
		eqm=[]
		for t in range(paras['ts'].min, paras['ts'].max):
			eqm.append(mean_squared_error(model[t:t+np.shape(data)[0], indexdata], data))
		delay = np.argmin(eqm)
		ts    = paras['ts'].min + delay
		paras['ts'].set(ts)

	# Dessin des courbes théoriques
	# fig = plt.figure(facecolor='w',figsize=figsize)
	# ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
	# ax.plot(model[ts:ts+len(data), 3], color='red',   label='model')
	# ax.plot(data,                 color='green', label='data', marker='x', ls='')
	# plt.legend()
	# plt.show()

	# Only R1 and F to calculate the residual, only on the window's size of the data
	result = model[ts:ts+np.shape(data)[0], indexdata] - data
	# result = (model[ts:ts+np.shape(data)[0], indexdata] - data).ravel()
	return result



if __name__ == '__main__':

	# Repertoire des figures
	repertoire = getRepertoire(True, './figures/simul/simulSEIR1R2F')
	prefFig    = repertoire + '/'

	# Solveur eq. diff.
	verbose = 1
	plot    = True
	N       = 66987244 # Population de la France
	dt      = 1
	solveur = SolveEDO_SEIR1R2F(N, dt, verbose)
	if verbose>0:
		print(solveur)

	# MAJ des parametres
	E0, I0, R10, R20, F0 = 0, 1, 0, 0, 0
	a, b, c, f, muI, xi  = 0.11, 0.24, 0.060, 0.05, 0.001, 0.0

	solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20,  F0=F0)
	solveur.modele.setParam(N=N,  a=a,   b=b,    c=c,     f=f  , muI=muI, xi=xi)

	# integration time grid
	T    = 750
	time = np.linspace(0, T-1, T)

	# Solveur to get the theoretical behavior
	########################################################################
	result = solveur.solveEDO(time)
	if plot==True:
		
		sliceedo = slice(0, 0+T)
		
		listePlot = [3,5]
		filename  = prefFig + 'SEIRFmodel_' + ''.join(map(str, listePlot)) + '.png'
		solveur.plotEDO(filename, '', sliceedo, sliceedo, plot=listePlot, data='', text=solveur.getTextParam())
		listePlot = [1,2,3,5]
		filename  = prefFig + 'SEIRFmodel_' + ''.join(map(str, listePlot)) + '.png'
		solveur.plotEDO(filename, '', sliceedo, sliceedo, plot=listePlot, data='', text=solveur.getTextParam())
		listePlot = [0,1,2,3,4,5]
		filename  = prefFig + 'SEIRFmodel_' + ''.join(map(str, listePlot)) + '.png'
		solveur.plotEDO(filename, '', sliceedo, sliceedo, plot=listePlot, data='', text=solveur.getTextParam())
