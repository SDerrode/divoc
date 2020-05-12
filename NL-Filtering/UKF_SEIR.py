#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from filterpy.kalman import UnscentedKalmanFilter as UKF
from filterpy.kalman import JulierSigmaPoints, MerweScaledSigmaPoints, rts_smoother
from filterpy.common import Q_discrete_white_noise
from common import readDataGouvFr, readDataEurope, drawAnnotation, Plot

from SEIR_Bacaer        import SEIR1R2_Bacaer
from Covid_SpecialDates import Covid_SpecialDates

# constante
verbose       = 1      # verbose level during run-time (0/1/2)
fileLocalCopy = False   # if we upload the file from the url (to get latest results) or from a local copy file


def getDates(country='France', verbose=False):

	Dates = None
	if country == 'France':
		Dates = Covid_SpecialDates(country=country, verbose=verbose)
		Dates.addConfDates('2020-03-16')
		Dates.addDeconfDates('2020-05-11')
		Dates.setListOtherDates(['2020-03-06'])
	if country == 'Germany':
		Dates = Covid_SpecialDates(country=country, verbose=verbose)
		Dates.addConfDates('2020-03-22')
		Dates.addDeconfDates('2020-04-20')

	if verbose>1:
			print('Interesting dates for ', country, '=', Dates)
	return Dates


if __name__ == '__main__':

	# Choose the country by uncommenting the country name you are interested in
	#############################################################################
	country = 'France'
	#country = 'Germany'

	# These are the date of confinement and deconfinement + other. See function getDates to add or delete dates to put the focus on
	Dates = getDates(country=country, verbose=verbose)
	
	# Lecture des données et copy of the observation
	#############################################################################

	startDate='2020-02-25'
	pd, z_observ = \
		readDataEurope(country=country, dateMin=startDate, dateMax=None, Plot=False, fileLocalCopy=fileLocalCopy, verbose=verbose)
	zs = []
	for z in pd[z_observ[0]]:
		zs.append(np.array([z]))

	# modele
	dt        = 1. # 1 day
	modele    = SEIR1R2_Bacaer(dt=dt)
	prefixFig = './figures/' + modele.modelShortName + '_' + country
	if modele.l == None:
		# estimation de lambda
		# correction de a # a = (l+c)*(1.+l/b)
		if verbose>1:
			print('lambda apres=', l, ' - a apres=', a)
			input('attente')
	if verbose>0:
		print('Modele', modele.modelName)
		print(modele)

	# Unscented KF
	#############################################################################

	sigmas = MerweScaledSigmaPoints(n=modele.dimState, alpha=.5, beta=2., kappa=modele.dimState-3.)
	ukf    = UKF(dim_x=modele.dimState, dim_z=modele.dimObs, fx=modele.f_bacaer, hx=modele.h_bacaer, dt=modele.dt, points=sigmas)
	# Filter init
	ukf.x[0] = modele.N   # on initialise S(0) à la population totale, le reste étant à 0 par défaut
	# Filter noise, starting with the measurment noise (R), and then the process noise (Q)
	ukf.Q[0:2,  0:2]  = Q_discrete_white_noise(2, dt=dt, var=1000)
	ukf.Q[2:4,  2:4]  = Q_discrete_white_noise(2, dt=dt, var=100)
	ukf.Q[4:6,  4:6]  = Q_discrete_white_noise(2, dt=dt, var=50)
	ukf.Q[6:8,  6:8]  = Q_discrete_white_noise(2, dt=dt, var=10)
	ukf.Q[8:10, 8:10] = Q_discrete_white_noise(2, dt=dt, var=10)
	ukf.R             = np.diag([10.])
	if verbose>1:
		print('ukf.R=', ukf.R)
		print('ukf.Q=', ukf.Q)
		input('pause')

	# UKF filtering, day by day
	# xMean_est = []
	# for z in zs:
	# 	ukf.predict()
	# 	ukf.update(z)
	# 	xMean_est.append(ukf.x.copy())
	# xMean_est = np.array(xMean_est)

	# UKF filtering and smoothing, batch mode
	xMean_est, _ = ukf.batch_filter(zs)
	# UKF filtering and smoothing, batch mode -- DOESN'T WORK!!!
	# uxs, ucovs = ukf.batch_filter(zs)
	# xMean_est, _, _ = ukf.rts_smoother(uxs, ucovs)
	
	# Panda dataframe preparation for plotting
	columns = ['$S(t)$', '$\dot{S}(t)$', '$E(t)$', '$\dot{E}(t)$', '$I(t)$', '$\dot{I}(t)$', '$R^1(t)$', \
				'$\dot{R}^1(t)$', '$R^2(t)$', '$\dot{R}^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$', '$R(t)=R^1(t)+R^2(t)$']
	for j in range(10):
		pd[columns[j]]=xMean_est[:, j]
	pd[columns[10]] = modele.N-np.sum(xMean_est[:, [0, 2, 4, 6]], axis=1)
	pd[columns[11]] = np.sum(xMean_est[:, [6, 8]], axis=1)
	if verbose>1:
		print(pd.tail())
	
	# save the genjerated data in csv form
	pd.to_csv (prefixFig + '_all.csv', header=True, index=False)

	# Graphiques
	#############################################################################
	
	# Plot de E, I, R^1, R^2
	Plot(pd, ['$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$'], country, prefixFig+'_EIR1R2.png', Dates, z_observ)

	# Plot de S
	Plot(pd, ['$S(t)$'], country, prefixFig+'_S.png', Dates)

	# Plot dérivées de de E, I, R^1, R
	Plot(pd, ['$\dot{E}(t)$', '$\dot{I}(t)$', '$\dot{R}^1(t)$', '$\dot{R}^2(t)$'], country, prefixFig+'_EIR1R2_dot.png', Dates)

	# Plot dérivée de S
	Plot(pd, ['$\dot{S}(t)$'], country, prefixFig+'_S_dot.png', Dates)

