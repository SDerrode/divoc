#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class SEIR1R2_Bacaer:
	''' Modèle SEIR1R2 du papier suivant 
		"Un modèle mathématique des débuts de l'épidémie de coronavirus en France", Nicolas Bacaër
	'''

	def __init__(self, dt=1):

		self.dimState       = 10        # Les 5 variables du modèle SEIR1R2 et leur dérivée, intercalées (S, \dot{S}, E, \dot{E}, ...)
		self.dimObs         = 1         # On observe uniquement R1 (nombre sumulé de cas confirmés)
		self.modelName      = 'SEIR1R2 (N. Bacaer)'  # Long name for the model
		self.modelShortName = 'SEIR1R2'              # Short name for the model

		# Default values for parameters - from the paper
		self.l  = 0.225          # coefficient de croissance exponentielle: exp(l * t)
		                         # Ce paramètre est estimé par la pente de la croube sur le début de l'épidémie
		self.N  = 65.E6          # taille de la population étudiée
		self.b  = 1./4.0         # phase de latence de 4 jours
		self.c  = 1./4.0         # durée moyenne dans le compartiment I
		self.f  = 0.2            # fraction d'individus infectieux qui sont comptabilisés parmi les cas confirmés au moment de l'isolement
		self.dt = dt
		# parameter 'a' is computed from some of the others
		self.set_a()

	def __str__(self):
		S  = 'Parameters: \n  l=' + str(self.l) + '\n  N='+ str(self.N) + '\n  b=' + str(self.b)
		S  +=  '\n  c=' + str(self.c) + '\n  f=' + str(self.f) + '\n  a=' + str(self.a)
		return S

	def setParam(self, l, N, b, c, f):
		self.l, self.N, self.b, self.c, self.f = l, N, b, c, f
		self.set_a()

	def set_a(self):
		# taux de contact effectif
		self.a = (self.l+self.c)*(1.+self.l/self.b) 

	def f_bacaer(self, x, dt):
		''' state transition function for Bacaer's model SEIR1R2'''

		# bidouille pour etre sur que la somme des valeurs fasse bien N -> à améliorer?
		x[[0,2,4,6,8]] /= abs(np.sum(x[[0,2,4,6,8]]))/self.N
		xout = np.empty_like(x)

		# param S, \dot{S}
		xout[0] = x[0]+dt*(-self.a/self.N * x[0] * x[4])
		xout[1] = x[1]
		# param E, \dot{E}
		xout[2] = x[2]+dt*(self.a/self.N * x[0] * x[4] - self.b * x[2])
		xout[3] = x[3]
		# param I, \dot{I}
		xout[4] = x[4]+dt*(self.b * x[2] - self.c * x[4])
		xout[5] = x[5]
		# param R^1, \dot{R}^1
		xout[6] = x[6]+dt*(self.f * self.c * x[4])
		xout[7] = x[7]
		# param R^2, \dot{R}^2
		xout[8] = x[8]+dt*((1.-self.f) * self.c * x[4])
		xout[9] = x[9]
		return xout

	def h_bacaer(self, x):
		return x[[6]] # on renvoie R1 (7ieme élément dans le vecteur x)

