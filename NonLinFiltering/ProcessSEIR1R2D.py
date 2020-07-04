#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from filterpy.kalman   import UnscentedKalmanFilter as UKF
from filterpy.kalman   import JulierSigmaPoints, MerweScaledSigmaPoints, rts_smoother
from filterpy.common   import Q_discrete_white_noise

from datetime          import datetime, timedelta
from sklearn.metrics   import mean_squared_error

from common            import readDataEurope, readDataFrance, readDates, addDaystoStrDate, getRepertoire
from common            import getNbDaysBetweenDateFromString, GetPairListDates
from France            import getPlace
from SolveEDO_SEIR1R2D import SolveEDO_SEIR1R2D

strDate = "%Y-%m-%d"

def fR1D(R1D, dt):
	return R1D # on renvoie R1 et F

def hR1D(R1D):
	return R1D # on renvoie R1 et F


def fit(sysargv):
	"""
		Program to process Covid Data.
 
		:Example:

		For countries (European database)
		>> python ProcessSEIR1R2D.py 
		>> python ProcessSEIR1R2D.py France 0 1 0 0 1 1
		>> python ProcessSEIR1R2D.py France 2 3 8 0 1 1          # 3 périodes pour les femmes en France avec un décalage de 11 jours
		>> python ProcessSEIR1R2D.py France,Germany 1 1 0 0 1 1 # 1 période pour les hommes francais et les hommes allemands 

		For French Region (French database)
		>> python ProcessSEIR1R2D.py FRANCE,D69         0 -1 11 0 1 1 # Code Insee Dpt 69 (Rhône)
		>> python ProcessSEIR1R2D.py FRANCE,R84         0 -1 11 0 1 1 # Tous les dpts de la Région dont le code Insee est 
		>> python ProcessSEIR1R2D.py FRANCE,R32+        0 -1 11 0 1 1 # Somme de tous les dpts de la Région 32 (Hauts-de-France)
		>> python ProcessSEIR1R2D.py FRANCE,MetropoleD  0 -1 11 0 1 1 # Tous les départements de la France métropolitaine
		>> python ProcessSEIR1R2D.py FRANCE,MetropoleD+ 0 -1 11 0 1 1 # Toute la France métropolitaine (en sommant les dpts)
		>> python ProcessSEIR1R2D.py FRANCE,MetropoleR+ 0 -1 11 0 1 1 # Somme des dpts de toutes les régions françaises
		Toute combinaison est possible de lieu : exemple FRANCE,R32+,D05,R84
		
		argv[1] : List of countries (ex. France,Germany,Italy), or see above.          Default: France 
		argv[2] : Sex (male:1, female:2, male+female:0). Only for french database             Default: 0 
		argv[3] : Periods ('1' -> 1 period ('all-in-on'), '!=1' -> severall periods).  Default: -1
		argv[4] : Shift of periods (in days).                                          Default: 11
		argv[5] : UKF filtering of data (0/1).                                         Default: 0
		argv[6] : Verbose level (debug: 3, ..., almost mute: 0).                       Default: 1
		argv[7] : Plot graphique (0/1).                                                Default: 1
	"""

	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	# Il y a 18 pays
	
	if len(sysargv) > 7:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Constantes
	######################################################@
	fileLocalCopy    = True  # if we upload the file from the url (to get latest data) or from a local copy file
	readStartDateStr = "2020-03-01" # "2020-03-01" Le 8 mars, pour inclure un grand nombre de pays européens dont la date de premier était postérieur au 1er mars
	readStopDateStr  = None
	recouvrement     = -1
	dt               = 1
	France           = 'France'
	thresholdSignif  = 1.5E-6

	# Interpetation of arguments - reparation
	######################################################@

	# Default value for parameters
	listplaces    = ['France']
	sexe, sexestr = 0, 'male+female'
	nbperiodes    = -1
	decalage      = 11
	UKF_filt      = False
	verbose       = 1
	plot          = True

	# Parameters from argv
	if len(sysargv)>0: liste = list(sysargv[0].split(','))
	if len(sysargv)>1: sexe = int(sysargv[1])
	if len(sysargv)>2: nbperiodes = int(sysargv[2])
	if len(sysargv)>3: decalage   = int(sysargv[3])
	if len(sysargv)>4 and int(sysargv[4])==1: UKF_filt = True
	if len(sysargv)>5: verbose    = int(sysargv[5])
	if len(sysargv)>6 and int(sysargv[6])==0: plot     = False
	if nbperiodes==1:       decalage = 0  # nécessairement pas de décalage (on compense le recouvrement)
	if sexe not in [0,1,2]:	sexe, sexestr = 0, 'male+female'      # sexe indiférencié
	if sexe == 1: 		          sexestr =    'male'
	if sexe == 2:                 sexestr =    'female'
	
	listplaces = []
	listnames  = []
	if liste[0]=='FRANCE':
		FrDatabase = True
		liste = liste[1:]
		for el in liste:
			l,n=getPlace(el)
			if el=='MetropoleR+':
				for l1,n1 in zip(l,n):
					listplaces.extend(l1)
					listnames.extend([n1])
			else:
				listplaces.extend(l)
				listnames.extend(n)
	else:
		listplaces = liste[:]
		FrDatabase = False

	if verbose>0:
		print('  Full command line : '+sysargv[0]+' '+str(nbperiodes)+' '+str(decalage)+' '+str(UKF_filt)+' '+str(verbose)+' '+str(plot), flush=True)

	# Data reading to get first and last date available in the data set
	######################################################@
	if FrDatabase == True: 
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr, _ = readDataFrance(['D69'],   readStartDateStr, readStopDateStr, fileLocalCopy, sexe, verbose=0)
	else:
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr, _ = readDataEurope(France, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
	dataLength = pd_exerpt.shape[0]

	readStartDate = datetime.strptime(readStartDateStr, strDate)
	if readStartDate<pd_exerpt.index[0]:
		readStartDate    = pd_exerpt.index[0]
		readStartDateStr = pd_exerpt.index[0].strftime(strDate)
	readStopDate = datetime.strptime(readStopDateStr, strDate)
	if readStopDate<pd_exerpt.index[-1]:
		readStopDate     = pd_exerpt.index[-1]
		readStopDateStr  = pd_exerpt.index[-1].strftime(strDate)

	dataLength = pd_exerpt.shape[0]
	if verbose>1:
		print('readStartDateStr=', readStartDateStr, ', readStopDateStr=', readStopDateStr)
		print('readStartDate   =', readStartDate,    ', readStopDate   =', readStopDate)
		print('dataLength      =', dataLength)
		#input('pause')

	# Collections of data return by this function
	modelSEIR1R2D = np.zeros(shape=(len(listplaces), dataLength, 6))
	data_deriv    = np.zeros(shape=(len(listplaces), dataLength, 2))
	modelR1_deriv = np.zeros(shape=(len(listplaces), dataLength, 2))
	data_all      = np.zeros(shape=(len(listplaces), dataLength, 2))
	modelR1_all   = np.zeros(shape=(len(listplaces), dataLength, 2))
	Listepd            = []
	ListetabParamModel = []

	# data observed
	data = np.zeros(shape=(dataLength, 2))

	# Paramètres sous forme de chaines de caractères
	ListeTextParam = []
	ListeDateI0    = []

	# Loop on the places to process
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			placefull   = 'France-' + listnames[indexplace][0]
			DatesString = readDates(France, verbose)
		else:
			placefull   = place
			DatesString = readDates(place, verbose)

		print('PROCESSING of', placefull, 'in', listnames)


		# data reading of the observations
		#############################################################################
		if FrDatabase == True:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr, dateFirstNonZeroStr = readDataFrance(place, readStartDateStr, readStopDateStr, fileLocalCopy, sexe, verbose=0)
		else:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr, dateFirstNonZeroStr = readDataEurope(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)

		shift_0value = getNbDaysBetweenDateFromString(readStartDateStr, dateFirstNonZeroStr)

		# UKF Filtering ?
		if UKF_filt == True:
			data2Filt = pd_exerpt[[HeadData[0], HeadData[1]]].to_numpy(copy=True)
			sigmas = MerweScaledSigmaPoints(n=2, alpha=.5, beta=2., kappa=1.) #1-3.)
			ukf    = UKF(dim_x=2, dim_z=2, fx=fR1D, hx=hR1D, dt=dt, points=sigmas)
			# Filter init
			ukf.x[:] = data2Filt[0, :]
			ukf.Q    = np.diag([30., 15.])
			ukf.R    = np.diag([170., 100.])
			if verbose>1:
				print('ukf.x[:]=', ukf.x[:])
				print('ukf.R   =', ukf.R)
				print('ukf.Q   =', ukf.Q)
			
			# UKF filtering and smoothing, batch mode
			R1Ffilt, _ = ukf.batch_filter(data2Filt)
			HeadData[0] += ' filt'
			HeadData[1] += ' filt'
			pd_exerpt[HeadData[0]] = R1Ffilt[:, 0]
			pd_exerpt[HeadData[1]] = R1Ffilt[:, 1]
			

		# Get the list of dates to process
		ListDates, ListDatesStr = GetPairListDates(readStartDate, readStopDate, DatesString, decalage+recouvrement, nbperiodes, recouvrement)
		if verbose>1:
			#print('ListDates   =', ListDates)
			print('ListDatesStr=', ListDatesStr)
			#input('pause')
		
		# Solveur edo
		solveur   = SolveEDO_SEIR1R2D(N, dt, verbose)
		indexdata = solveur.indexdata
		E0, I0, R10, R20, D0 = 0, 1, 0, 0, 0
		
		# Repertoire des figures
		if plot==True:
			repertoire = getRepertoire(UKF_filt, './figures/SEIR1R2D_UKFilt/'+placefull+'/sexe_'+str(sexe)+'_shift_'+str(decalage), './figures/SEIR1R2D/'+placefull+'/sexe_'+str(sexe)+'_shift_'+str(decalage))
			prefFig    = repertoire+'/Process_'
		
		# Remise à 0 des données
		data.fill(0.)


		# Boucle pour traiter successivement les différentes fenêtres
		###############################################################
		
		ListeTextParamPlace     = []
		ListetabParamModelPlace = []
		ListeEQM                = []

		DEGENERATE_CASE = False

		for i in range(len(ListDatesStr)):

			# dates of the current period
			fitStartDate,    fitStopDate    = ListDates[i]
			fitStartDateStr, fitStopDateStr = ListDatesStr[i]

			# Est-on dans un CAS degénéré?
			# print(getNbDaysBetweenDateFromString(dateFirstNonZeroStr, fitStopDateStr))
			if getNbDaysBetweenDateFromString(dateFirstNonZeroStr, fitStopDateStr)<5: # Il faut au moins 5 données pour fitter
				DEGENERATE_CASE = True

			if i>0:
				DatesString.addOtherDates(fitStartDateStr)

			# Récupérations des données observées
			dataLengthPeriod = 0
			indMinPeriod     = (fitStartDate-readStartDate).days

			for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), HeadData[0]]):
				data[indMinPeriod+j, 0] = z
				dataLengthPeriod +=1
			for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), HeadData[1]]):
				data[indMinPeriod+j, 1] = z
			slicedata      = slice(indMinPeriod, indMinPeriod+dataLengthPeriod)
			slicedataderiv = slice(slicedata.start+1, slicedata.stop)
			if verbose>0:
				print('  dataLength      =', dataLength)
				print('  indMinPeriod    =', indMinPeriod)
				print('  dataLengthPeriod=', dataLengthPeriod)
				print('  fitStartDateStr =', fitStartDateStr)
				print('  fitStopDateStr  =', fitStopDateStr)
				#input('attente')
			

			# Set initialisation data for the solveur
			############################################################################

			# paramètres initiaux à optimiser
			if i==0:
				datelegend = fitStartDateStr
				# ts=getNbDaysBetweenDateFromString(DatesString.listFirstCaseDates[0], readStartDateStr)
				# En premiere approximation, on prend la date du premier cas estimé pour le pays (même si c'est faux pour les régions et dpts)
				ts=getNbDaysBetweenDateFromString(DatesString.listFirstCaseDates[0], dateFirstNonZeroStr)
				if ts<0:
					continue # On passe au pays suivant
				if nbperiodes!=1: # pour plusieurs périodes
					#l, b0, c0, f0 = 0.255, 1./5.2, 1./12, 0.08
					#a0 = (l+c0)*(1.+l/b0)
					#a0, b0, c0, f0, mu0, xi0 = 0.55, 0.34, 0.12, 0.25, 0.0005, 0.0001
					a0, b0, c0, f0, mu0, xi0 = 0.60, 0.55, 0.30, 0.50, 0.0005, 0.0001
					T  = 150
				else: # pour une période
					#a0, b0, c0, f0, mu0, xi0  = 0.10, 0.29, 0.10, 0.0022, 0.00004, 0.
					a0, b0, c0, f0, mu0, xi0  = 0.70, 0.25, 0.05, 0.003, 0.0005, 0.0001
					T = 350
				
			if i==1 or i==2:
				datelegend = None

				_, a0, b0, c0, f0, mu0, xi0 = solveur.modele.getParam()
				R10 = int(data[indMinPeriod, 0]) # on corrige R1 à la valeur numérique 
				F0  = int(data[indMinPeriod, 1]) # on corrige F à la valeur numérique 
				if i == 1:
					a0 /= 4. # le confinement réduit drastiquement (pour aider l'optimisation) 
				T  = 120
				ts = 0

			time = np.linspace(0, T-1, T)

			solveur.modele.setParam(N=N,  a=a0,  b=b0,   c=c0,    f=f0, mu=mu0, xi=xi0)
			solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20, D0=D0)

			# Before optimization
			###############################

			# Solve ode avant optimization
			sol_ode = solveur.solveEDO(time)
			# calcul time shift initial (ts) with respect to data
			if i==0:
				ts = solveur.compute_tsfromEQM(data[slicedata, :], T, indexdata)
			else:
				solveur.TS = ts = 0
			sliceedo = slice(ts, min(ts+dataLengthPeriod, T))
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))

			# plot
			if plot==True and DEGENERATE_CASE==False:
				titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - [' + fitStartDateStr + '\u2192' + addDaystoStrDate(fitStopDateStr, -1) + '] (Sex=', + sexestr + 'Shift=' + str(decalage) + ')'
				listePlot = indexdata
				filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Init.png'
				solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(datelegend, Period=i))

			# Parameters optimization
			############################################################################

			solveur.paramOptimization(data[slicedata, :], time) # version lorsque ts est calculé automatiquement
			#solveur.paramOptimization(data[slicedata, :], time, ts) # version lorsque l'on veut fixer ts
			_, a1, b1, c1, f1, mu1, xi1 = solveur.modele.getParam()
			R0 = solveur.modele.getR0()
			if verbose>0:
				print('Solver''s state after optimization=', solveur)
				print('  Reproductivité après: ', R0)

			# After optimization
			###############################
			
			# Solve ode avant optimization
			sol_ode = solveur.solveEDO(time)
			# calcul time shift with respect to data
			if i==0:
				ts = solveur.compute_tsfromEQM(data[slicedata, :], T, indexdata)
			else:
				solveur.TS = ts = 0
			sliceedo      = slice(ts, min(ts+dataLengthPeriod, T))
			sliceedoderiv = slice(sliceedo.start+1, sliceedo.stop)
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))
			if i==0: # on se souvient de la date du premier infesté
				dateI0 = addDaystoStrDate(fitStartDateStr, -ts+shift_0value)
				if verbose>2:
					print('dateI0=', dateI0)
					input('attente')

			# sauvegarde des param (tableau et texte)
			seuil = (data[slicedata.stop-1, 0]-data[slicedata.start, 0])/getNbDaysBetweenDateFromString(fitStartDateStr, fitStopDateStr)/N
			#print('seuil=', seuil)
			#print('DEGENERATE_CASE=', DEGENERATE_CASE)
			if DEGENERATE_CASE==True:
				ROsignificatif = False
				ListetabParamModelPlace.append([-1., -1., -1., -1., -1., -1., -1.])
			else:
				if seuil<thresholdSignif:
					ROsignificatif = False
					ListetabParamModelPlace.append([a1, b1, c1, f1, mu1, xi1, -1.])
				else:
					ROsignificatif = True
					ListetabParamModelPlace.append([a1, b1, c1, f1, mu1, xi1, R0])
				# print('seuil=', seuil)
				# print('ROsignificatif=', ROsignificatif)
				# print('R0=', R0)
				# input('pause')
			
			ListeTextParamPlace.append(solveur.getTextParamWeak(datelegend, ROsignificatif, Period=i))
			
			data_deriv_period    = (data[slicedataderiv, :]           - data   [slicedataderiv.start-1:slicedataderiv.stop-1, :])        / dt
			modelR1_deriv_period = (sol_ode[sliceedoderiv, indexdata] - sol_ode[sliceedoderiv.start-1 :sliceedoderiv.stop-1, indexdata]) / dt
			data_all_period      = data[slicedataderiv, :]
			modelR1_all_period   = sol_ode[sliceedoderiv, indexdata]

			if plot==True and DEGENERATE_CASE==False:
				titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - [' + fitStartDateStr + '\u2192' + addDaystoStrDate(fitStopDateStr, -1) + '] (Sex=', + sexestr + 'Shift=' + str(decalage) + ')'

				# listePlot = [0,1,2,3,4,5]
				# filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + '.png'
				# solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(datelegend, ROsignificatif, Period=i))
				listePlot = [1,2,3,5]
				filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Final.png'
				solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(datelegend, ROsignificatif, Period=i))
				listePlot = indexdata
				filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Final.png'
				solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(datelegend, ROsignificatif, Period=i))

				# dérivée  numérique de R1 et F
				filename = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Deriv.png'
				solveur.plotEDO_deriv(filename, titre, sliceedoderiv, slicedataderiv, data_deriv_period, indexdata, text=solveur.getTextParam(datelegend, ROsignificatif, Period=i))

			# sol_ode_withSwitch = solveur.solveEDO_withSwitch(T, timeswitch=ts+dataLengthPeriod)

			# ajout des données dérivées
			data_all     [indexplace, slicedataderiv, :] = data_all_period
			modelR1_all  [indexplace, slicedataderiv, :] = modelR1_all_period
			data_deriv   [indexplace, slicedataderiv, :] = data_deriv_period
			modelR1_deriv[indexplace, slicedataderiv, :] = modelR1_deriv_period

			# ajout des SEIR1R2D
			modelSEIR1R2D[indexplace, slicedata.start:slicedata.stop, :] = sol_ode[ts:ts+sliceedo.stop-sliceedo.start, :]

			# preparation for next iteration
			_, E0, I0, R10, R20, D0 = map(int, sol_ode[ts+dataLengthPeriod+recouvrement, :])
			#print('A LA FIN : E0, I0, R10, R20, D0=', E0, I0, R10, R20, D0)

			if verbose>1:
				input('next step')
		
		Listepd.append(pd_exerpt)
		ListeDateI0.append(dateI0)

		# calcul de l'EQM sur les données (et non sur les dérivées des données)
		#EQM = mean_squared_error(data_deriv[indexplace, :], modelR1_deriv[indexplace, :])
		EQM = mean_squared_error(data_all[indexplace, :], modelR1_all[indexplace, :])
		ListeEQM.append(EQM)

		# udpate des listes pour transmission
		ListeTextParam.append(ListeTextParamPlace)
		ListetabParamModel.append(ListetabParamModelPlace)

	return modelSEIR1R2D, ListeTextParam, Listepd, data_deriv, modelR1_deriv, ListetabParamModel, ListeEQM, ListeDateI0


if __name__ == '__main__':
	modelSEIR1R2D, ListeTextParam, Listepd, data_deriv, modelR1_deriv, ListetabParamModel, ListeEQM, ListeDateI0 = fit(sys.argv[1:])
