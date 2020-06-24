#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from filterpy.kalman  import UnscentedKalmanFilter as UKF
from filterpy.kalman  import JulierSigmaPoints, MerweScaledSigmaPoints, rts_smoother
from filterpy.common  import Q_discrete_white_noise

from datetime         import datetime, timedelta
from sklearn.metrics  import mean_squared_error

from common           import readDataEurope, readDataFrance, readDates, addDaystoStrDate, getRepertoire
from common           import getNbDaysBetweenDateFromString, GetPairListDates, getListPlaces
from SolveEDO_SEIR1R2 import SolveEDO_SEIR1R2

strDate = "%Y-%m-%d"

def fR1(r1, dt):
	return r1 # on renvoie R1

def hR1(r1):
	return r1 # on renvoie R1


def fit(sysargv):
	"""
		Program to process Covid Data.
 
		:Example:

		For countries (European database)
		>> python3 ProcessSEIR1R2.py 
		>> python3 ProcessSEIR1R2.py France 1 0 0 1 1
		>> python3 ProcessSEIR1R2.py France 3 8 0 1 1
		>> python3 ProcessSEIR1R2.py France,Germany 1 0 0 1 1

		For French Region (French database)
		>> python3 ProcessSEIR1R2.py France,69        3 8 0 1 1 # Dpt 69 (Rhône)
		>> python3 ProcessSEIR1R2.py France,all       3 8 0 1 1 # Tous les dpts francais
		>> python3 ProcessSEIR1R2.py France,metropole 3 8 0 1 1 # Tous les dpts francais de la métropole uniquement
		>> python3 ProcessSEIR1R2.py France,69,01     3 8 1 1 1 # Dpt 69 (Rhône) + Dpt 01 (Ain) with UKF filtering
		
		argv[1] : Country (or list separeted by ','), or 'France' followed by a list of dpts. Default: France 
		argv[2] : Periods ('1' -> 1 period ('all-in-on'), '!=1' -> severall periods).         Default: 1
		argv[3] : Shift of periods (in days).                                                 Default: 0
		argv[4] : UKF filtering of data (0/1).                                                Default: 0
		argv[5] : Verbose level (debug: 3, ..., almost mute: 0).                              Default: 1
		argv[6] : Plot graphique (0/1).                                                       Default: 1
	"""
	
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	# Il y a 18 pays

	if len(sysargv) > 6:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Constantes
	######################################################@
	fileLocalCopy    = True  # if we upload the file from the url (to get latest data) or from a local copy file
	readStartDateStr = "2020-03-08" # "2020-03-01" Le 8 mars, pour inclure un grand nombre de pays européens dont la date de premier était postérieur au 1er mars
	readStopDateStr  = None
	recouvrement     = -1
	dt               = 1
	France           = 'France'

	# Interpetation of arguments - reparation
	######################################################@

	# Default value for parameters
	listplaces = ['France']
	nbperiodes = 1
	decalage   = 2
	UKF_filt   = False
	verbose    = 1
	plot       = True

	# Parameters from argv
	######################################@

	if len(sysargv)>0: listplaces = list(sysargv[0].split(','))
	FrDatabase, listplaces = getListPlaces(listplaces)
	if len(sysargv)>1: nbperiodes    = int(sysargv[1])
	if len(sysargv)>2: decalage      = int(sysargv[2])
	if len(sysargv)>3 and int(sysargv[3])==1: UKF_filt = True
	if len(sysargv)>4: verbose       = int(sysargv[4])
	if len(sysargv)>5 and int(sysargv[5])==0: plot     = False
	if nbperiodes==1: 
		decalage = 0  # nécessairement pas de décalage (on compense le recouvrement)

	if verbose>0:
		print('  Full command line : '+sysargv[0]+' '+str(nbperiodes)+' '+str(decalage)+' '+str(UKF_filt)+' '+str(verbose)+' '+str(plot), flush=True)
	

	# Data reading to get first and last date available in the data set
	######################################################
	if FrDatabase == True: 
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance('69',   readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
	else:
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(France, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
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
		print('dataLength=', dataLength)
		#input('pause')

	# Collections of data return by this function
	modelSEIR1R2  = np.zeros(shape=(len(listplaces), dataLength, 5))
	data_deriv    = np.zeros(shape=(len(listplaces), dataLength, 1))
	modelR1_deriv = np.zeros(shape=(len(listplaces), dataLength, 1))
	data_all      = np.zeros(shape=(len(listplaces), dataLength, 1))
	modelR1_all   = np.zeros(shape=(len(listplaces), dataLength, 1))
	Listepd            = []
	ListetabParamModel = []

	# data observed
	data = np.zeros(shape=(dataLength, 1))

	# Paramètres sous forme de chaines de caractères
	ListeTextParam = [] 
	ListeDateI0    = []

	# Loop on the places to process
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			placefull   = 'France-' + place
			DatesString = readDates(France, verbose)
		else:
			placefull   = place
			DatesString = readDates(place, verbose)

		print('PROCESSING of', placefull, 'in', listplaces)
		

		# data reading of the observations
		#############################################################################
		if FrDatabase == True:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
		else:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)

		# UKF Filtering ?
		if UKF_filt == True:
			data2filter = pd_exerpt[HeadData[2]].tolist()
			sigmas = MerweScaledSigmaPoints(n=1, alpha=.5, beta=2., kappa=0.)
			ukf    = UKF(dim_x=1, dim_z=1, fx=fR1, hx=hR1, dt=dt, points=sigmas)
			# Filter init
			ukf.x[0] = data2filter[0]
			ukf.Q    = np.diag([30.])
			ukf.R    = np.diag([170.])
			if verbose>1:
				print('ukf.x[0]=', ukf.x[0])
				print('ukf.R   =', ukf.R)
				print('ukf.Q   =', ukf.Q)
			
			# UKF filtering and smoothing, batch mode
			R1filt, _ = ukf.batch_filter(data2filter)
			HeadData[2] = HeadData[2] + ' filt'
			pd_exerpt[HeadData[2]] = R1filt

		# Get the list of dates to process
		ListDates, ListDatesStr = GetPairListDates(readStartDate, readStopDate, DatesString, decalage+recouvrement, nbperiodes, recouvrement)
		if verbose>1:
			#print('ListDates   =', ListDates)
			print('ListDatesStr=', ListDatesStr)
			#input('pause')
		
		# Solveur edo
		solveur   = SolveEDO_SEIR1R2(N, dt, verbose)
		indexdata = solveur.indexdata
		E0, I0, R10, R20 = 0, 1, 0, 0

		# Repertoire des figures
		if plot==True:
			repertoire = getRepertoire(UKF_filt, './figures/SEIR1R2_UKFilt/'+placefull + '/' + str(decalage), './figures/SEIR1R2/' + placefull+ '/' + str(decalage))
			prefFig    = repertoire  + '/Process_'
		
		# Remise à 0 des données
		data.fill(0.)


		# Boucle pour traiter successivement les différentes fenêtres
		###############################################################
		
		ListeTextParamPlace     = []
		ListetabParamModelPlace = []
		ListeEQM                = []

		for i in range(len(ListDatesStr)):

			# dates of the current period
			fitStartDate,    fitStopDate    = ListDates[i]
			fitStartDateStr, fitStopDateStr = ListDatesStr[i]

			if i>0:
				DatesString.addOtherDates(fitStartDateStr)

			# Récupérations des données observées
			dataLengthPeriod = 0
			indMinPeriod     = (fitStartDate-readStartDate).days

			for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), (HeadData[2])]):
				data[indMinPeriod+j, 0] = z
				dataLengthPeriod +=1
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
				ts=getNbDaysBetweenDateFromString(DatesString.listFirstCaseDates[0], readStartDateStr)
				if ts<0:
					continue # On passe au pays suivant
				if nbperiodes!=1: # pour plusieurs périodes
					#l, b0, c0, f0 = 0.255, 1./5.2, 1./12, 0.08
					#a0 = (l+c0)*(1.+l/b0)
					#a0, b0, c0, f0 = 0.55, 0.34, 0.12, 0.25
					a0, b0, c0, f0 = 0.60, 0.55, 0.30, 0.50
					T  = 150
				else: # pour une période
					#a0, b0, c0, f0  = 0.35, 0.29, 0.075, 0.0022
					a0, b0, c0, f0  = 0.70, 0.25, 0.05, 0.003
					T = 350
				# R20 = int((1.-f0)*(R10/f0))
				
			if i==1 or i==2:
				_, a0, b0, c0, f0 = solveur.modele.getParam()
				R10 = int(data[indMinPeriod, 0]) # on corrige R1 à la valeur numérique 
				if i == 1:
					a0 /= 4. # le confinement réduit drastiquement (pour aider l'optimisation) 
				T  = 120
				ts = 0

			time = np.linspace(0, T-1, T)

			solveur.modele.setParam(N=N,  a=a0,  b=b0,   c=c0,    f=f0)
			solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20)

			# Before optimization
			###############################

			# Solve ode avant optimization
			sol_ode = solveur.solveEDO(time)
			# calcul time shift initial (ts) with respect to data
			if i==0:
				ts         = solveur.compute_tsfromEQM(data[slicedata, :], T, indexdata)
			else:
				solveur.TS = ts
			sliceedo = slice(ts, min(ts+dataLengthPeriod, T))
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))
			
			# plot
			if plot==True:
				titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - [' + fitStartDateStr + '\u2192' + addDaystoStrDate(fitStopDateStr, -1) + '] (Shift=' + str(decalage) + ')'
				listePlot = [3]
				filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Init.png'
				solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(fitStartDateStr))

			# Parameters optimization
			############################################################################

			if i==0:
				solveur.paramOptimization(data[slicedata, :], time) # version lorsque ts est calculé automatiquement
			else:
				solveur.paramOptimization(data[slicedata, :], time, ts) # version lorsque l'on veut fixer ts
			_, a1, b1, c1, f1 = solveur.modele.getParam()
			R0=solveur.modele.getR0()
			if verbose>0:
				print('Solver''s state after optimization=', solveur)
				print('  Reproductivité après: ', R0)

			# After optimization
			###############################
			
			# Solve ode apres optimization
			sol_ode = solveur.solveEDO(time)
			# calcul time shift with respect to data
			if i==0:
				ts            = solveur.compute_tsfromEQM(data[slicedata, :], T, indexdata)
			else:
				ts = solveur.TS
			#print('ts=', ts)
			sliceedo      = slice(ts, min(ts+dataLengthPeriod, T))
			sliceedoderiv = slice(sliceedo.start+1, sliceedo.stop)
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))
			if i==0: # on se souvient de la date du premier infesté
				dateI0 = addDaystoStrDate(fitStartDateStr, -ts)
				if verbose>2:
					print('dateI0=', dateI0)
					input('attente')

			# Sauvegarde des param (tableau et texte)
			ROsignificatif=True
			threshold = (data[slicedata.stop-1, 0]-data[slicedata.start, 0])/getNbDaysBetweenDateFromString(fitStartDateStr, fitStopDateStr)
			if  threshold <1.0: # moins de 1 cas détecté par jour sur la période 3
				ROsignificatif = False
				ListetabParamModelPlace.append([a1, b1, c1, f1, -1.])
			else:
				ListetabParamModelPlace.append([a1, b1, c1, f1, R0])
	
			ListeTextParamPlace.append(solveur.getTextParamWeak(fitStartDateStr, ROsignificatif))
			
			data_deriv_period    = (data[slicedataderiv, :]           - data   [slicedataderiv.start-1:slicedataderiv.stop-1, :])        / dt
			modelR1_deriv_period = (sol_ode[sliceedoderiv, indexdata] - sol_ode[sliceedoderiv.start-1 :sliceedoderiv.stop-1, indexdata]) / dt
			data_all_period      = data[slicedataderiv, :]
			modelR1_all_period   = sol_ode[sliceedoderiv, indexdata]

			if plot==True:
				titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - [' + fitStartDateStr + '\u2192' + addDaystoStrDate(fitStopDateStr, -1) +'] (Shift=' + str(decalage) + ')'
				
				# listePlot =[0,1,2,3,4]
				# filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + '.png'
				# solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(fitStartDateStr, ROsignificatif))
				listePlot =[1,2,3]
				filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Final.png'
				solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(fitStartDateStr, ROsignificatif))
				listePlot =[3]
				filename  = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Final.png'
				solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(fitStartDateStr, ROsignificatif))

				# dérivée  numérique de R1
				filename = prefFig + str(decalage) + '_Period' + str(i) + '_' + ''.join(map(str, listePlot)) + 'Deriv.png'
				solveur.plotEDO_deriv(filename, titre, sliceedoderiv, slicedataderiv, data_deriv_period, indexdata, text=solveur.getTextParam(fitStartDateStr, ROsignificatif))

			# sol_ode_withSwitch = solveur.solveEDO_withSwitch(T, timeswitch=ts+dataLengthPeriod)

			# ajout des données et des données dérivées
			data_all     [indexplace, slicedataderiv, :] = data_all_period
			modelR1_all  [indexplace, slicedataderiv, :] = modelR1_all_period
			data_deriv   [indexplace, slicedataderiv, :] = data_deriv_period
			modelR1_deriv[indexplace, slicedataderiv, :] = modelR1_deriv_period

			# ajout des SEIR1R2
			modelSEIR1R2[indexplace, slicedata.start:slicedata.stop, :] = sol_ode[ts:ts+sliceedo.stop-sliceedo.start, :]

			# preparation for next iteration
			_, E0, I0, R10, R20 = map(int, sol_ode[ts+dataLengthPeriod+recouvrement, :])
			#print('A LA FIN : E0, I0, R10, R20=', E0, I0, R10, R20)

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

	return modelSEIR1R2, ListeTextParam, Listepd, data_deriv, modelR1_deriv, ListetabParamModel, ListeEQM, ListeDateI0


if __name__ == '__main__':
	modelSEIR1R2, ListeTextParam, Listepd, data_deriv, modelR1_deriv, ListetabParamModel, ListeEQM, ListeDateI0 = fit(sys.argv[1:])
