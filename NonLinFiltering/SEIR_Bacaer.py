#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class SEIR1R2_Bacaer:
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
		self.f  = 0.003    # france 0.0032 , RU : 0.007      # fraction d'individus infectieux qui sont comptabilisés parmi les cas confirmés au moment de l'isolement
		self.dt = dt

		self.colorCycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
		
		#a, b, c, f = 0.9025, 0.68, 0.43, 0.4
		#a, b, c, f = 0.155, 1./5.2, 1./12.39, 0.05

	def __str__(self):
		S  = 'Model ' + self.modelName + ':' + ' \n  N=' + str(self.N) 
		S += '\n  a='+ str(np.round(self.a, decimals=3))
		S += '\n  b=' + str(np.round(self.b, decimals=3)) + '\n  c=' + str(np.round(self.c, decimals=3)) + '\n  f=' + str(np.round(self.f, decimals=3))
		return S

	def getTextParam(self):
		S  = 'Model '  + self.modelShortName + ':'
		S += '\n  a='  + str(np.round(self.a, decimals=3)) + ', b=' + str(np.round(self.b, decimals=3))
		S += '\n  c='  + str(np.round(self.c, decimals=3)) + ', f=' + str(np.round(self.f, decimals=3))
		if self.c!= 0.:
			S += '\n  R0=' + str(np.round(self.a/self.c, decimals=3))
		return S

	def setParam(self, N, a, b, c, f):
		self.N, self.a, self.b, self.c, self.f = N, a, b, c, f

	def setN(self, N):
		self.N = N
	def seta(self, a):
		self.a = a
	def setb(self, b):
		self.b = b
	def setc(self, c):
		self.c = c
	def setf(self, f):
		self.f = f

	def getParam(self):
		return (self.N, self.a, self.b, self.c, self.f)

	def h_bacaer(self, x):
		return x[[3]] # on renvoie R1 (4ieme élément dans le vecteur SEIR1R2)

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
		if indice == 0: return '$S(t)$'
		if indice == 1: return '$E(t)$'
		if indice == 2: return '$I(t)$'
		if indice == 3: return '$R^1(t)$'
		if indice == 4: return '$R^2(t)$'
		print('PB getString - indice =', indice, ' does not exist!')
		exit(1)

	def getColor(self, indice):
		if indice >= 0 and indice<5: return self.colorCycle[indice]
		print('PB getColor - indice =', indice, ' does not exist!')
		exit(1)

	def getColorFromString(self, string):
		if string == '$S(t)$'  : return self.colorCycle[0] 
		if string == '$E(t)$'  : return self.colorCycle[1] 
		if string == '$I(t)$'  : return self.colorCycle[2] 
		if string == '$R^1(t)$': return self.colorCycle[3] 
		if string == '$R^2(t)$': return self.colorCycle[4] 

		return 'black'

	def f_bacaer(self, x, dt):
		'''State transition function for Bacaer's model SEIR1R2'''

		# petite normalisation pour véiter des dérives
		x /= abs(np.sum(x))/self.N
		ret = odeint(self.deriv, x, [0, dt], args=self.getParam())
		return ret.T[:, -1]
