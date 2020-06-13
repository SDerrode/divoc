#!/usr/bin/env python3
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
		
		# Total population, N.
		self.N = N
		# time shift with respect to some data
		self.TS = -1

		# Modèle d'eq. diff non linéaires
		self.dt = dt # 1 day

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

	def getParamInit(self):
		return self.y0

	def compute_tsfromEQM(self, data, T, indexdata):
		
		dataLength = np.shape(data)[0]
		eqm = np.zeros(shape=(T-dataLength))
		for t in range(T-dataLength):
			# print(np.shape(self.solution[t:t+dataLength, indexdata]))
			# print(np.shape(data))
			# input('apuse')
			eqm[t] = mean_squared_error(self.solution[t:t+dataLength, indexdata], data)
			# print('eqm[t]=', eqm[t])
			# input('pause')
		self.TS = np.argmin(eqm)
		# print('self.TS=', self.TS)
		# input('apsue')
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


	# def simulEDO(self, T):
	# 	self.solution = np.zeros(shape=(T, self.nbparam), dtype=int) # 5 parameters: S, E, I, R1, R2

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
		
		labels=[r'$\frac{\partial R^1(t)}{\partial t}$', r'$\frac{\partial F(t)}{\partial t}$']
		q=0
		for k in indexdata:
			ax.plot(time, deriv_sol3[sliceedoderiv.start:min(sliceedoderiv.stop, np.shape(deriv_sol3)[0]), q], color=self.modele.getColor(int(indexdata[q])), alpha=1.0, lw=2, label=labels[q])
			q += 1
		# seconde estimation de R2
		# R2bis = self.N-self.solution[:, 0]-self.solution[:, 1]-self.solution[:, 2]-self.solution[:, 3]
		# ax.plot(timeedo, R2bis, color=self.modele.getColor(indice), alpha=1.0, lw=2, label='$R^2(t)=N-\sum{SEIR^1}$')

		# les données observées
		labels=[r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial F(n)}{\partial n}$']
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

