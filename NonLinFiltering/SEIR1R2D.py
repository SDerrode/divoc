#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from SEIR1R2 import SEIR1R2

class SEIR1R2D(SEIR1R2):
	''' Modèle SEIR1R2D du papier suivant 
		"Universal masking is urgent in the covid-19 pandemic...", 
		De Kai, Guy-Philippe Goldstein, et al. ArXiv
	'''

	def __init__(self, N, dt=1, verbose=1):

		super().__init__(N, dt, verbose)

		# Correction of some varaible specific
		self.dimState       = 6          # Les 6 variables du modèle SEIR1R2D (S, E, ...)
		self.dimObs         = 2          # On observe R1 et F 
		self.modelName      = 'SEIR1R2D' # Long name for the model
		self.modelShortName = 'SEIR1R2D' # Short name for the model

		# Add suplementary parameters - from the paper
		self.mu = 0.001         # taux de mortalité spécifique à la pandémie
		self.xi = 0.001         # rate of re-susceptibility

		self.setR0()       # MAJ de R0

	def __str__(self):
		S  = super().__str__()
		S += '\n  mu=' + str(np.round(self.mu, decimals=3))
		S += '\n  xi=' + str(np.round(self.xi, decimals=3))
		return S

	def setR0(self):
		if self.c+self.mu != 0.:
			self.R0 = (self.a+self.xi)/(self.c+self.mu)
		else:
			self.R0 = -1.

	def getTextParam(self, ROsignificatif=True, Period=1):
		S  = 'Model '  + self.modelShortName +  ' - Period ' + str(Period) + ':'
		S += '\n' + r'  $a='    + str(np.round(self.a, decimals=4))  + r', b='   + str(np.round(self.b, decimals=4))  + '$'
		S += '\n' + r'  $c='    + str(np.round(self.c, decimals=4))  + r', f='   + str(np.round(self.f, decimals=4))  + '$'
		S += '\n' + r'  $\mu='  + str(np.round(self.mu, decimals=5)) + r', \xi=' + str(np.round(self.xi, decimals=5)) + '$'
		if self.c!= 0. and ROsignificatif==True:
			S += '\n' + r'  $R_0=' + str(np.round(self.R0, decimals=2)) + '$'
		if ROsignificatif==False:
			S += '\n' + r'  $R_0$ non significatif'
		return S

	def setParam(self, N, a, b, c, f, mu, xi):
		super().setParam(N, a, b, c, f)
		self.mu, self.xi = mu, xi
		self.setR0()       # MAJ de R0
	
	def setf(self, mu):
		self.mu = mu
		self.setR0()       # MAJ de R0
	def setf(self, xi):
		self.xi = xi
		self.setR0()       # MAJ de R0

	def getParam(self):
		return (self.N, self.a, self.b, self.c, self.f, self.mu, self.xi)
	def getR0(self):
		return self.R0

	# The SEIR1R2D model differential equations.
	def deriv(self, y, t, N, a, b, c, f, mu, xi):
		S, E, I, R1, R2, D = y
		R = R1+R2
		dSdt  = -a * S * I / N + xi * R
		dEdt  = a * S * I / N - b * E
		dIdt  = b * E - c * I - mu * I
		dR1dt = f      * (c * I - xi * R)
		dR2dt = (1.-f) * (c * I - xi * R)
		dDdt  = mu * I
		return dSdt, dEdt, dIdt, dR1dt, dR2dt, dDdt

	def getString(self, indice):
		string = super().getString(indice)
		if indice == 5: return r'$D(t)$'
		return string

	def getColor(self, indice):
		if indice >= 0 and indice<6: return self.colorCycle[indice]
		print('PB getColor - indice =', indice, ' does not exist!')
		exit(1)

	def getColorFromString(self, string):
		col = super().getColorFromString(string)
		if string == r'$D(t)$'  : return self.colorCycle[5] 
		return col
