#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from SEIR1R2 import SEIR1R2

class SEIR1R2F(SEIR1R2):
	''' Modèle SEIR1R2F du papier suivant 
		"Universal masking is urgent in he covid-19 pandemic...", 
		De Kai, Guy-Philippe Goldstein, et al. ArXiv
	'''

	def __init__(self, N, dt=1, verbose=1):

		super().__init__(N, dt, verbose)

		# Correction of some varaible specific
		self.dimState       = 6          # Les 6 variables du modèle SEIR1R2F (S, E, ...)
		self.dimObs         = 2          # On observe R1 et F 
		self.modelName      = 'SEIRF'    # Long name for the model
		self.modelShortName = 'SEIR1R2F' # Short name for the model

		# Add suplementary parameters - from the paper
		self.muI = 0.001         # taux de mortalité
		self.xi  = 0.0           # rate of re-susceptibility

	def __str__(self):
		S  = super().__str__()
		S += '\n  muI=' + str(np.round(self.muI, decimals=3))
		S += '\n  xi =' + str(np.round(self.xi, decimals=3))
		return S

	def getTextParam(self):
		S  = 'Model '  + self.modelShortName + ':'
		S += '\n' + r'  $a='  + str(np.round(self.a, decimals=4)) + r', b=' + str(np.round(self.b, decimals=4)) + '$'
		S += '\n' + r'  $c='  + str(np.round(self.c, decimals=4)) + r', f=' + str(np.round(self.f, decimals=4)) + '$'
		S += '\n' + r'  $\mu_I='  + str(np.round(self.muI, decimals=5)) + r', \xi=' + str(np.round(self.xi, decimals=5)) + '$'
		if self.c!= 0.:
		 	S += '\n' + r'  $R_0=' + str(np.round(self.a/self.c, decimals=2)) + '$'
		return S

	def setParam(self, N, a, b, c, f, muI, xi):
		super().setParam(N, a, b, c, f)
		self.muI, self.xi = muI, xi
	
	def setf(self, muI):
		self.muI = muI
	def setf(self, xi):
		self.xi = xi

	def getParam(self):
		return (self.N, self.a, self.b, self.c, self.f, self.muI, self.xi)

	# The SEIR1R2F model differential equations.
	def deriv(self, y, t, N, a, b, c, f, muI, xi):
		S, E, I, R1, R2, F = y
		R = R1+R2
		dSdt  = -a * S * I / N + xi * R
		dEdt  = a * S * I / N - b * E
		dIdt  = b * E - c * I - muI * I
		dR1dt = f * (c * I - xi * R)
		dR2dt = (1.-f) * (c * I - xi * R)
		dFdt  = muI * I
		return dSdt, dEdt, dIdt, dR1dt, dR2dt, dFdt

	def getString(self, indice):
		string = super().getString(indice)
		if indice == 5: return r'$F(t)$'
		return string

	def getColor(self, indice):
		if indice >= 0 and indice<6: return self.colorCycle[indice]
		print('PB getColor - indice =', indice, ' does not exist!')
		exit(1)

	def getColorFromString(self, string):
		col = super().getColorFromString(string)
		if string == r'$F(t)$'  : return self.colorCycle[5] 
		return col

	# def fx(self, x, dt):
	# 	'''State transition function for Bacaer's model SEIR1R2'''

	# 	# petite normalisation pour éviter des dérives
	# 	input('verifier si cette deriv est reelle - SEIR1R2F')
	# 	x /= abs(np.sum(x))/self.N
	# 	ret = odeint(self.deriv, x, [0, dt], args=self.getParam())
	# 	return ret.T[:, -1]

	# def hx(self, x):
	# 	return x[[3, 5]] # on renvoie R1 et F (4ieme et 6ième éléments dans le vecteur SEIR1R2F)

