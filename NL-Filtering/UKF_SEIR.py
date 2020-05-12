#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from filterpy.kalman    import UnscentedKalmanFilter as UKF
from filterpy.kalman    import JulierSigmaPoints, MerweScaledSigmaPoints, rts_smoother
from filterpy.common    import Q_discrete_white_noise

from common             import readDataGouvFr, readDataEurope, drawAnnotation, Plot
from SEIR_Bacaer        import SEIR1R2_Bacaer
from Covid_SpecialDates import Covid_SpecialDates

# constante
fileLocalCopy = False         # if we upload the file from the url (to get latest results) or from a local copy file
startDate     = '2020-02-25'  # whatever the upload data starts, this set the satrt date to be processed


def main():
	"""
		Program to perform UKF filtering on Covid Data.
 
		:Example:

		>> python3 UKF_SEIR.py 
		>> python3 UKF_SEIR.py Italy
		>> python3 UKF_SEIR.py Italy 1 1
		>> python3 UKF_SEIR.py France,Germany 1 1 # Shortcut for the two countries successively
		>> python3 UKF_SEIR.py Europe 1 1         # Sum of all coutries in Europe (no time-shift of data)
		>> python3 UKF_SEIR.py World 1 1         # Sum of all coutries in Europe (no time-shift of data)

		argv[1] : Country (or list separeted by ','), or Continent, or 'World'. Default: France 
		argv[2] : Verbose level (debug: 3, ..., almost mute: 0). Default: 1
		argv[3] : Plot graphique (0/1).                          Default: 1
	"""

	print('Command line : ', sys.argv, flush=True)
	if len(sys.argv) > 4:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	listcountries = ['France']
	verbose       = 1
	plot          = True

	# Parameters from argv
	if len(sys.argv)>1: listcountries = list(sys.argv[1].split(','))
	if len(sys.argv)>2: verbose = int(sys.argv[2])
	if len(sys.argv)>3 and int(sys.argv[3])==0: plot = False

	# Modele d'eq. diff non lineaires
	#############################################################################
	dt     = 1. # 1 day
	modele = SEIR1R2_Bacaer(dt=dt)
	if modele.l == None:
		# estimation de lambda
		# correction de a # a = (l+c)*(1.+l/b)
		if verbose>1:
			print('lambda apres=', l, ' - a apres=', a)
			input('attente')
	if verbose>0:
		print('Modele', modele.modelName)
		print(modele)

	
	for country in listcountries:

		print('PROCESSING of ', country, ' in ', listcountries)

		prefixFig = './figures/' + modele.modelShortName + '_' + country

		# Lecture des données et copy of the observation
		#############################################################################

		pd, z_observ = readDataEurope(country=country, dateMin=startDate, dateMax=None, \
								Plot=plot, fileLocalCopy=fileLocalCopy, verbose=verbose)
		if pd.empty==True:
			continue

		zs = []
		for z in pd[z_observ[0]]:
			zs.append(np.array([z]))

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
		
		# Panda dataframe
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

		if plot == True:

			# These are the date of confinement and deconfinement + other. See function getDates to add or delete dates to put the focus on
			Dates = getDates(country, verbose)
			print('Dates=', Dates)

			# Plot de E, I, R^1, R^2
			Plot(pd, ['$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$'], country, prefixFig+'_EIR1R2.png', Dates, z_observ)

			# Plot de S
			Plot(pd, ['$S(t)$'], country, prefixFig+'_S.png', Dates)

			# Plot dérivées de de E, I, R^1, R
			Plot(pd, ['$\dot{E}(t)$', '$\dot{I}(t)$', '$\dot{R}^1(t)$', '$\dot{R}^2(t)$'], country, prefixFig+'_EIR1R2_dot.png', Dates)

			# Plot dérivée de S
			Plot(pd, ['$\dot{S}(t)$'], country, prefixFig+'_S_dot.png', Dates)


def getDates(country, verbose):

	Dates = None
	if country == 'France':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates('2020-03-16')
		Dates.addDeconfDates('2020-05-11')
		Dates.setListOtherDates(['2020-03-06'])
	if country == 'Germany':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates('2020-03-22')
		Dates.addDeconfDates('2020-04-20')
	if country == 'Italy':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates('2020-03-09')
		Dates.addDeconfDates('2020-05-04')
	if country == 'Spain':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates('2020-03-15')
		Dates.addDeconfDates('2020-05-10')

	if verbose>1:
		print(Dates)
	
	return Dates


if __name__ == '__main__':
	main()
