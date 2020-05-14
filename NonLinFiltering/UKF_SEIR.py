#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from filterpy.kalman    import UnscentedKalmanFilter as UKF
from filterpy.kalman    import JulierSigmaPoints, MerweScaledSigmaPoints, rts_smoother
from filterpy.common    import Q_discrete_white_noise

from common             import readDataEurope, getDates, Plot
from SolveDiff_SEIR1R2  import SolveDiff_SEIR1R2

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
		>> python3 UKF_SEIR.py France,Germany 1 1 # Shortcut for processing the two countries successively
		>> python3 UKF_SEIR.py Europe 1 1         # Sum of all coutries in Europe (no time-shift of data)
		>> python3 UKF_SEIR.py World 1 1          # Sum of all coutries in Europe (no time-shift of data)

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
	if len(sys.argv)>2: verbose       = int(sys.argv[2])
	if len(sys.argv)>3 and int(sys.argv[3])==0: plot = False

	# Modele d'eq. diff non lineaires
	#############################################################################
	# Solveur eq. diff.
	N       = 65.E6
	dt      = 1
	solveur = SolveDiff_SEIR1R2(N, dt, verbose)
	if verbose>0:
		print(solveur)

	for country in listcountries:

		print('PROCESSING of ', country, ' in ', listcountries)
		prefixFig = './figures/' + solveur.modele.modelShortName + '_' + country

		# Lecture des données et copy of the observation
		#############################################################################

		pd, z_observ = readDataEurope(country=country, dateMin=startDate, dateMax=None, \
								plot=plot, fileLocalCopy=fileLocalCopy, verbose=verbose)
		if pd.empty==True:
			continue # on passe au pays suivant si celui-ci n'existe pas

		zs = []
		for z in pd[z_observ[0]]:
			zs.append(np.array([z]))

		# Unscented KF
		#############################################################################
		sigmas = MerweScaledSigmaPoints(n=solveur.modele.dimState, alpha=.5, beta=2., kappa=solveur.modele.dimState-3.)
		ukf    = UKF(dim_x=solveur.modele.dimState, dim_z=solveur.modele.dimObs, fx=solveur.modele.f_bacaer, hx=solveur.modele.h_bacaer, dt=solveur.modele.dt, points=sigmas)
		# Filter init
		ukf.x[1] = zs[0]
		ukf.x[2] = zs[0]
		ukf.x[3] = zs[0]
		ukf.x[0] = solveur.modele.N-np.sum(ukf.x[1:])
		# Filter noise, starting with the measurment noise (R), and then the process noise (Q)
		# ukf.Q[0:1, 0:1] = 1000 #Q_discrete_white_noise(1, dt=dt, var=1000)
		# ukf.Q[1:2, 1:2] = 50   #Q_discrete_white_noise(1, dt=dt, var=100)
		# ukf.Q[2:3, 2:3] = 50   #Q_discrete_white_noise(1, dt=dt, var=50)
		# ukf.Q[3:4, 3:4] = 10   #Q_discrete_white_noise(1, dt=dt, var=10)
		# ukf.Q[4:5, 4:5] = 10   #Q_discrete_white_noise(1, dt=dt, var=10)
		ukf.Q           = np.diag([1., 1., 1., 1., 1.])
		ukf.R           = np.diag([1.])
		if verbose>1:
			print('ukf.x[1:]=', ukf.x[1:])
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
		columns = ['$S(t)$', '$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$', '$R(t)=R^1(t)+R^2(t)$']
		for j in range(5):
			pd[columns[j]]=xMean_est[:, j]
		pd[columns[5]] = solveur.modele.N-np.sum(xMean_est[:, 0:4], axis=1)
		pd[columns[6]] = np.sum(xMean_est[:, 3:5], axis=1)
		if verbose>1:
			print(pd.tail())
		
		# save the genjerated data in csv form
		pd.to_csv (prefixFig + '_all.csv', header=True, index=False)

		if plot == True:

			# These are the date of confinement and deconfinement + other. See function getDates to add or delete dates to put the focus on
			Dates = getDates(country, verbose)
			#print('Dates=', Dates)

			titre = country + ' - ' + solveur.modele.modelName

			# Plot de E, I, R^1, R^2
			Plot(pd, ['$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$'], titre , prefixFig+'_EIR1R2.png', Dates, z_observ)

			# Plot de S
			Plot(pd, ['$S(t)$'], titre, prefixFig+'_S.png', Dates)

			# Plot dérivées de de E, I, R^1, R^2
			# Plot(pd, ['$\dot{E}(t)$', '$\dot{I}(t)$', '$\dot{R}^1(t)$', '$\dot{R}^2(t)$'], titre, prefixFig+'_EIR1R2_dot.png', Dates)

			# Plot dérivée de S
			# Plot(pd, ['$\dot{S}(t)$'], titre, prefixFig+'_S_dot.png', Dates)





if __name__ == '__main__':
	main()
