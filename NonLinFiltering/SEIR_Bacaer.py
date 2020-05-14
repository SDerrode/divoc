#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import odeint

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
		self.a  = 0.155
		self.b  = 1./5.2         # phase de latence de 4 jours
		self.c  = 1./12.39       # durée moyenne dans le compartiment I
		self.f  = 0.2            # fraction d'individus infectieux qui sont comptabilisés parmi les cas confirmés au moment de l'isolement
		self.dt = dt

		#a, b, c, f = 0.9025, 0.68, 0.43, 0.4
		#a, b, c, f = 0.155, 1./5.2, 1./12.39, 0.05

	def __str__(self):
		S  = 'Parameters for ' + self.modelName + ' model: \n  N=' + str(self.N) + '\n  a='+ str(self.a)
		S += '\n  b=' + str(self.b) + '\n  c=' + str(self.c) + '\n  f=' + str(self.f)
		return S

	def setParam(self, a, b, c, f):
		self.a, self.b, self.c, self.f = a, b, c, f

	def getParam(self):
		return (self.a, self.b, self.c, self.f)

	def h_bacaer(self, x):
		return x[[3]] # on renvoie R1 (4ieme élément dans le vecteur SEIR1R2)

	# The SEIR1R2 model differential equations.
	def deriv(self, y, t, a, b, c, f):
	    S, E, I, R1, R2 = y
	    dSdt  = -a * S * I / self.N
	    dEdt  = a * S * I / self.N - b * E
	    dIdt  = b * E - c * I
	    dR1dt = f * c * I
	    dR2dt = (1.-f) * c * I
	    return dSdt, dEdt, dIdt, dR1dt, dR2dt

	def f_bacaer(self, x, dt):
		'''State transition function for Bacaer's model SEIR1R2'''

		# petite normalisation pour véiter des dérives
		x /= abs(np.sum(x))/self.N
		ret = odeint(self.deriv, x, [0, dt], args=(self.N, self.a, self.b, self.c, self.f))
		#xout = ret.T[:, -1]
		# print('xout=', xout, 'shpe=', np.shape(xout))
		#input('attente')
		return ret.T[:, -1]
