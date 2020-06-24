#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class SEIR1R2:
	''' Modèle SEIR1R2 du papier suivant 
		"Un modèle mathématique des débuts de l'épidémie de coronavirus en France", Nicolas Bacaër
	'''

	def __init__(self, N, dt=1, verbose=1):

		self._verbose=verbose

		self.dimState       = 5   # Les 5 variables du modèle SEIR1R2 (S, E, ...)
		self.dimObs         = 1   # On observe uniquement R1 (nombre sumulé de cas confirmés)
		self.modelName      = 'SEIR1R2 (N. Bacaer)'  # Long name for the model
		self.modelShortName = 'SEIR1R2'              # Short name for the model

		# Default values for parameters - from the paper
		self.N  = N              # taille de la population étudiée
		self.a  = 0.155			 # tau de contact effectif
		self.b  = 1./5.2         # phase de latence de 5.2 jours
		self.c  = 1./12.39       # durée moyenne dans le compartiment I
		self.f  = 0.003          # france 0.0032 , RU : 0.007      # fraction d'individus infectieux qui sont comptabilisés parmi les cas confirmés au moment de l'isolement
		self.dt = dt

		if 'SEIR1R2.SEIR1R2' in str(type(self)):
			self.setR0()       # MAJ de R0
		
		self.colorCycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

	def __str__(self):
		S  = 'Model ' + self.modelName + ':' + ' \n  N=' + str(self.N) 
		S += '\n  a='+ str(np.round(self.a, decimals=3))
		S += '\n  b=' + str(np.round(self.b, decimals=3)) + '\n  c=' + str(np.round(self.c, decimals=3)) + '\n  f=' + str(np.round(self.f, decimals=3))
		return S

	def getTextParam(self, ROsignificatif=True):
		S  = 'Model '  + self.modelShortName + ':'
		S += '\n' + r'  $a='  + str(np.round(self.a, decimals=4)) + r', b=' + str(np.round(self.b, decimals=4)) + '$'
		S += '\n' + r'  $c='  + str(np.round(self.c, decimals=4)) + r', f=' + str(np.round(self.f, decimals=4)) + '$'
		if self.c!= 0. and ROsignificatif==True:
			S += '\n' + r'  $R_0=' + str(np.round(self.R0, decimals=2)) + '$'
		if ROsignificatif==False:
			S += '\n' + r'  $R_0$ non significatif'
		return S

	def setR0(self):
		if self.c != 0.:
			self.R0 = self.a/self.c
		else:
			self.R0 = -1.

	def setParam(self, N, a, b, c, f):
		self.N, self.a, self.b, self.c, self.f = N, a, b, c, f
		self.setR0()       # MAJ de R0

	def setN(self, N):
		self.N = N
	def seta(self, a):
		self.a = a
		self.setR0()       # MAJ de R0
	def setb(self, b):
		self.b = b
	def setc(self, c):
		self.c = c
		self.setR0()       # MAJ de R0
	def setf(self, f):
		self.f = f

	def getParam(self):
		return (self.N, self.a, self.b, self.c, self.f)
	def getR0(self):
		return self.R0

	# The SEIR1R2 model differential equations.
	def deriv(self, y, t, N, a, b, c, f):
		S, E, I, R1, R2 = y
		dSdt  = -a * S * I / N
		dEdt  = a * S * I / N - b * E
		dIdt  = b * E - c * I
		dR1dt = f * c * I
		dR2dt = (1.-f) * c * I
		return dSdt, dEdt, dIdt, dR1dt, dR2dt

	def getString(self, indice):
		if indice == 0: return r'$S(t)$'
		if indice == 1: return r'$E(t)$'
		if indice == 2: return r'$I(t)$'
		if indice == 3: return r'$R^1(t)$'
		if indice == 4: return r'$R^2(t)$'
		return indice

	def getColor(self, indice):
		if indice >= 0 and indice<5: return self.colorCycle[indice]
		print('PB getColor - indice =', indice, ' does not exist!')
		exit(1)

	def getColorFromString(self, string):
		if string == r'$S(t)$'  : return self.colorCycle[0] 
		if string == r'$E(t)$'  : return self.colorCycle[1] 
		if string == r'$I(t)$'  : return self.colorCycle[2] 
		if string == r'$R^1(t)$': return self.colorCycle[3] 
		if string == r'$R^2(t)$': return self.colorCycle[4] 
		return string
