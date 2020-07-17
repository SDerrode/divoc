#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator

from scipy.integrate import odeint
from lmfit           import minimize, Parameters, Parameter, report_fit
from sklearn.metrics import mean_squared_error

from SEIR1R2         import SEIR1R2
from common          import addDaystoStrDate

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)


class SolveEDO:

	def __init__(self, N, dt=1, verbose=1):

		self.verbose = verbose
		
		# Total population
		self.N = N
		# Time delay
		self.TS = -1

		self.dt = dt

	def __str__(self):
		S  = 'Solver:\n  N=' + str(self.N) + '\n  dt='+str(self.dt) + '\n  y0='+ str(self.y0) + '\n  TS='+ str(self.TS) + '\n'
		S = self.modele.__str__()
		return S

	def getTextParam(self, startDate=None, ROsignificatif=True, Degenerate_case=False, Period=1):
		if Degenerate_case==False:
			S = self.modele.getTextParam(ROsignificatif, Period=Period)
		else:
			S = 'Degenerate case - Period ' + str(Period)
		S += '\nSolveur init:\n'# + '\n  N=' + str(self.N) + '\n  dt='+str(self.dt) \
		S += r'  $S_0=' + str(self.y0[0]) + '$\n'
		S += r'  $E_0=' + str(self.y0[1]) + r', I_0='+str(self.y0[2]) + '$\n'
		S += r'  $R^1_0=' + str(self.y0[3]) + r', R^2_0=' + str(self.y0[4])+ '$'
		if startDate!=None:
			dateI0 = addDaystoStrDate(startDate, -self.TS)
			S += '\n  Date ' + r'$d_{I=1}$' + ':'+dateI0
		return S

	def getTextParamWeak(self, startDate=None, ROsignificatif=True, Period=1):
		S  = self.modele.getTextParam(ROsignificatif, Period=Period)
		if startDate!=None:
			dateI0 = addDaystoStrDate(startDate, -self.TS)
			S += '\n  Date ' + r'$d_{I=1}$' + ':'+dateI0
		return S

	def getParamInit(self):
		return self.y0

	def compute_tsfromEQM(self, data, T, indexdata):
		dataLength = np.shape(data)[0]
		eqm = np.zeros(shape=(T-dataLength))
		for t in range(T-dataLength):
			eqm[t] = mean_squared_error(self.solution[t:t+dataLength, indexdata], data)
		self.TS = np.argmin(eqm)
		return self.TS

	def solveEDO_withSwitch(self, T, timeswitch=-1):
		'''
		Integrate the EDO equations over the time grid t, one by one; to deal with parameter switching
		'''

		self.solution = np.zeros(shape=(T, self.nbparam))
		ret = odeint(self.modele.deriv, self.y0, [0, 0], args=self.modele.getParam())
		self.solution[0, :] = ret.T[:, -1]
		for i in range(1, T):
			if i == timeswitch:
				self.modele.seta(0.)
				self.modele.setb(0.)
			ret = odeint(self.modele.deriv, self.solution[i-1], [0, self.dt], args=self.modele.getParam())
			self.solution[i, :] = ret.T[:, -1]

		return self.solution

	def solveEDO(self, time):
		# Integrate the EDO equations over the time grid t.
		self.solution = odeint(self.modele.deriv, self.y0, time, args=self.modele.getParam())
		return self.solution


	def plotEDO(self, name, title, sliceedo, slicedata, plot, data='', text=''):

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
			ax.plot(time, data[slicedata, 0], color=self.modele.getColor(3), alpha=1.0, lw=2, label='Data $R^1$', marker='x', ls='')
			if np.shape(data)[1]>1:
				ax.plot(time, data[slicedata, 1], color=self.modele.getColor(5), alpha=1.0, lw=2, label='Data $F$', marker='+', ls='')

		ax.set_xlabel('Time (days)')
		ax.set_ylabel('Pop. size (' + str(int(self.N)) + ')')
		ax.yaxis.set_tick_params(length=0)
		ax.xaxis.set_tick_params(length=0)
		ax.grid(b=True, which='major', c='w', lw=1, ls='-')
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))
		
		legend = ax.legend()
		legend.get_frame().set_alpha(0.5)
		for spine in ('top', 'right', 'bottom', 'left'):
			ax.spines[spine].set_visible(False)

		# ajout d'un text d'annotation
		if text != '':
			bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
			ax.annotate(text, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=8, bbox=bbox, ha="left", va="center") 
		
		plt.tight_layout(rect=(0, 0, 1., 0.95))
		plt.title(title)
		plt.savefig(name, dpi=dpi)
		plt.close()


	def plotEDO_deriv(self, name, title, sliceedoderiv, slicedataderiv, deriv_data='', indexdata=None, text=''):

		time = np.linspace(slicedataderiv.start, slicedataderiv.stop-1, slicedataderiv.stop-slicedataderiv.start)

		deriv_sol3 = (self.solution[1:, indexdata] - self.solution[:-1, indexdata]) / self.dt

		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
		
		labels=[r'$\frac{\partial R^1(t)}{\partial t}$', r'$\frac{\partial D(t)}{\partial t}$']
		q=0
		for k in indexdata:
			ax.plot(time, deriv_sol3[sliceedoderiv.start:min(sliceedoderiv.stop, np.shape(deriv_sol3)[0]), q], color=self.modele.getColor(int(indexdata[q])), alpha=1.0, lw=2, label=labels[q])
			q += 1
		# seconde estimation de R2
		# R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]
		# ax.plot(timeedo, R2bis, color=self.modele.getColor(indice), alpha=1.0, lw=2, label='$R^2(t)=N-\sum{SEIR^1}$')

		# les données observées
		labels=[r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial D(n)}{\partial n}$']
		markers=['x', '+']
		if len(deriv_data) != 0:
			q=0
			for k in indexdata:
				ax.plot(time, deriv_data[:, q], color=self.modele.getColor(int(indexdata[q])), alpha=1.0, lw=0.5, label=labels[q], marker=markers[q])
				q += 1

		ax.set_xlabel('Time (days)')
		ax.yaxis.set_tick_params(length=0)
		ax.xaxis.set_tick_params(length=0)
		ax.grid(b=True, which='major', c='w', lw=1, ls='-')
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))
		
		legend = ax.legend()
		legend.get_frame().set_alpha(0.5)
		for spine in ('top', 'right', 'bottom', 'left'):
			ax.spines[spine].set_visible(False)

		# ajout d'un text d'annotation
		if text != '':
			bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
			ax.annotate(text, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=8, bbox=bbox, ha="left", va="center") 
		
		plt.tight_layout(rect=(0, 0, 1., 0.95))
		plt.title(title)
		plt.savefig(name, dpi=dpi)
		plt.close()

