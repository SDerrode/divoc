#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt
from datetime          import datetime, timedelta
from sklearn.metrics   import mean_squared_error

from Common            import readDataEurope, readDataFrance, getDates, Plot, addDaystoStrDate
from Common            import getLowerDateFromString, getNbDaysBetweenDateFromString
from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2

strDate = "%Y-%m-%d"

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
		>> python3 ProcessSEIR1R2.py France,69    3 8 0 1 1 # Dpt 69 (Rhône)
		>> python3 ProcessSEIR1R2.py France,69,01 3 8 0 1 1 # Dpt 69 (Rhône) + Dpt 01 (Ain)
		
		argv[1] : Country (or list separeted by ','), or 'France' followed by a list of departments. Default: France 
		argv[2] : Periods (1 -> 1 period ('all-in-on'), 'x!=1' -> severall periods).                 Default: 1
		argv[3] : Shift of periods (in days).                                                        Default: 0
		argv[4] : Compensation strategy (0/1).                                                       Default: 0
		argv[5] : Verbose level (debug: 3, ..., almost mute: 0).                                     Default: 1
		argv[6] : Plot graphique (0/1).                                                              Default: 1
	"""

	# Interpetation of arguments - reparation
	######################################################@

	print('Command line : ', sysargv, flush=True)
	if len(sysargv) > 6:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	listplaces = ['France']
	nbperiodes = 1
	decalage   = 0
	surplus    = False
	verbose    = 1
	plot       = True

	# Parameters from argv
	if len(sysargv)>1: listplaces = list(sysargv[0].split(','))
	FrDatabase = False
	if listplaces[0]=='France' and len(listplaces)>1:
		try:
			int(listplaces[1]) #If this is a number, then it is a french dpt
		except Exception as e:
			FrDatabase=False
		else:
			FrDatabase = True
			listplaces = listplaces[1:]

	if len(sysargv)>1: nbperiodes    = int(sysargv[1])
	if len(sysargv)>2: decalage      = int(sysargv[2])
	if len(sysargv)>3 and int(sysargv[3])==1: surplus = True
	if len(sysargv)>4: verbose       = int(sysargv[4])
	if len(sysargv)>5 and int(sysargv[5])==0: plot = False
	if nbperiodes==1: 
		decalage   = 0  # nécessairement pas de décalage
	if surplus==True: 
		decalage   = 0  # nécessairement pas de décalage
		nbpériodes = 2  # necessairement plusieurs périodes

	# Constantes
	######################################################@
	fileLocalCopy    = True  # if we upload the file from the url (to get latest data) or from a local copy file
	readStartDateStr = "2020-02-01" #"2020-03-01" 
	readStopDateStr  = None
	recouvrement     = -3
	dt               = 1

	# Data reading to get first and last date available in the data set
	######################################################@
	if FrDatabase == True: 
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance('69', readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
	else:
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope('France', readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
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
		#input('pause')



	# Collections of data return by this function
	modelSEIR1R2  = np.zeros(shape=(len(listplaces), dataLength, 5))
	data_deriv    = np.zeros(shape=(len(listplaces), dataLength))
	modelR1_deriv = np.zeros(shape=(len(listplaces), dataLength))
	Listepd            = []
	ListetabParamModel = []

	# Surplus and correction from previou period
	data_surplus   = np.zeros(shape=(dataLength))
	data_corrected = np.zeros(shape=(dataLength))

	# Paramètres sou form de chaines
	ListeTextParam = [] 

	# Loop on the place to process
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			placefull   = 'France-' + place
			DatesString = getDates('France', verbose)
		else:
			placefull   = place
			DatesString = getDates(place, verbose)

		print('PROCESSING of', placefull, 'in', listplaces)
		

		# data reading of the observations
		#############################################################################
		if FrDatabase == True:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
		else:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)

		# Get the list of dates to process
		ListDates, ListDatesStr = GetListDates(readStartDate, readStopDate, DatesString, decalage, nbperiodes, recouvrement)
		if verbose>1:
			print('ListDates   =', ListDates)
			print('ListDatesStr=', ListDatesStr)
			#input('pause')
		
		# Solveur edo
		solveur = SolveEDO_SEIR1R2(N, dt, verbose)

		# Constantes
		import os
		repertoire = './figures/'+ placefull
		if not os.path.exists(repertoire):
			os.makedirs(repertoire)
		prefFig = repertoire + '/' + solveur.modele.modelShortName + '_' + placefull
		
		columns = [r'$S(t)$', r'$E(t)$', r'$I(t)$', r'$R^1(t)$', r'$R^2(t)$', r'$R^2(t)=N-\sum(SEIR^1)$', r'$R(t)=R^1(t)+R^2(t)$']
		E0, I0, R10, R20 = 0, 1, 0, 0
		
		# Remise à 0 du surplus et de la correction
		data_surplus.fill(0.)
		data_corrected.fill(0.)


		# Boucle pour traiter successivement les différentes fenêtres
		###############################################################
		
		ListeTextParamPlace     = []
		ListetabParamModelPlace = []
		ListeEQM = []

		for i in range(len(ListDatesStr)):

			# dates of the current period
			fitStartDate,    fitStopDate    = ListDates[i]
			fitStartDateStr, fitStopDateStr = ListDatesStr[i]

			if i>0:
				DatesString.addOtherDates(fitStartDateStr)

			# Récupérations des données observées
			dataLengthPeriod = 0
			indMinPeriod     = (fitStartDate-readStartDate).days

			for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), (HeadData[0])]):
				data_corrected[indMinPeriod+j] = z - data_surplus[indMinPeriod+j]
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
				if nbperiodes!=1: # pour plusieurs périodes
					l, b0, c0, f0 = 0.255, 1./5.2, 1./12, 0.08
					a0 = (l+c0)*(1.+l/b0)
					T  = 150
				else: # pour une période
					a0, b0, c0, f0  = 0.35, 0.29, 0.075, 0.0022
					T  = 350

			if i==1 or i==2:
				R10 = int(data_corrected[indMinPeriod]) # on corrige R1 à la valeur numérique 
				_, a0, b0, c0, f0   = solveur.modele.getParam()
				if i == 1:
					a0 /= 3 # le confinement réduit drastiquement a (pour aider l'optimisation) 
				T = 120

			time = np.linspace(0, T-1, T)

			solveur.modele.setParam(N=N,  a=a0,  b=b0,   c=c0,    f=f0)
			solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20)

			# Before optimization
			###############################

			# Solve ode avant optimization
			sol_ode = solveur.solve_SEIR1R2(time)
			# calcul time shift initial (ts) with respect to data
			ts       = solveur.compute_tsfromEQM(data_corrected[slicedata], T)
			sliceedo = slice(ts, min(ts+dataLengthPeriod, T))
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))
			
			# plot
			if plot==True:
				listePlot = [3]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fitinit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, placefull, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))

			# Parameters optimization
			############################################################################

			solveur.paramOptimization(data_corrected[slicedata], time)
			_, a1, b1, c1, f1 = solveur.modele.getParam()
			if c1 != 0.:
				R0=a1/c1
			else:
				R0 = -1.
			if verbose>0:
				print('Solver''s state after optimization=', solveur)
				if c1 != 0.:
					print('  Reproductivité après: ', R0)

			# After optimization
			###############################
			
			# Solve ode avant optimization
			sol_ode = solveur.solve_SEIR1R2(time)
			# calcul time shift with respect to data
			ts            = solveur.compute_tsfromEQM(data_corrected[slicedata], T)
			sliceedo      = slice(ts, min(ts+dataLengthPeriod, T))
			sliceedoderiv = slice(sliceedo.start+1, sliceedo.stop)
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))

			# sauvegarde des param (tableau et texte)
			ListetabParamModelPlace.append([a1, b1, c1, f1, R0])
			ListeTextParamPlace    .append(solveur.getTextParam(readStartDateStr))
			
			data_deriv_period    = (data_corrected[slicedataderiv] - data_corrected[slicedataderiv.start-1:slicedataderiv.stop-1  ]) / dt
			modelR1_deriv_period = (sol_ode[sliceedoderiv, 3]      - sol_ode       [sliceedoderiv.start-1 :sliceedoderiv.stop-1, 3]) / dt
			# print('len(data_deriv_period)   =', len(data_deriv_period))
			# print('len(modelR1_deriv_period)=', len(modelR1_deriv_period))
			#input('pause')

			if plot==True:
				titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - Shift=' + str(decalage)
				
				listePlot =[0,1,2,3,4]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, titre, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
				listePlot =[1,2,3]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, titre, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
				listePlot =[3]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, titre, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))

				# dérivée  numérique de R1
				filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_3Deriv.png'
				solveur.plot_dR1(filename, titre, sliceedoderiv, slicedataderiv, data_deriv_period, text=solveur.getTextParam(readStartDateStr))

			# sol_ode_withSwitch = solveur.solve_SEIR1R2_withSwitch(T, timeswitch=ts+dataLengthPeriod)
			# if plot==True:
			#     titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - Shift=' + str(decalage)
			#
			#     listePlot=[0,1,2,3,4]
			#     filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
			#     solveur.plot_SEIR1R2(filename, titre, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
			#     listePlot=[1,2,3]
			#     filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
			#     solveur.plot_SEIR1R2(filename, titre, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
			#     listePlot=[3]
			#     filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
			#     solveur.plot_SEIR1R2(filename, titre, sliceedo, slicedata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))

			#     begDate = addDaystoStrDate(fitStartDate, -int(ts))
			#     #print('begDate=', begDate)
			#     for j in range(5):
			#         pd_exerpt.loc[begDate:, (columns[j])]=sol_ode_withSwitch[0:dataLengthPeriod+ts, j]
			#     pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(sol_ode_withSwitch[0:dataLengthPeriod+ts, 0:4], axis=1)
			#     pd_exerpt.loc[begDate:, (columns[6])] = np.sum(sol_ode_withSwitch[0:dataLengthPeriod+ts, 3:5], axis=1)
			#     if verbose>1:
			#         print(pd_exerpt.tail())

			#     titre = placefull + ' - ' + solveur.modele.modelName
			#     # Plot de E, I, R^1, R^2
			#     Plot(pd_exerpt, titre, prefFig+'_Fit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)
			
			# if surplus == True:
			#     # get the residuals
			#     bout  = min(T, ts+dataLengthPeriod +(dataLength-indMinPeriod+dataLengthPeriod))
			#     debut = ts+dataLengthPeriod
			#     # print('bout=', bout)
			#     # print('debut=', debut)
			#     # print('data_corrected[indMinPeriod+dataLengthPeriod-1]=', data_corrected[indMinPeriod+dataLengthPeriod-1])
			#     data_surplus[indMinPeriod+dataLengthPeriod:indMinPeriod+dataLengthPeriod+(bout-debut)] += sol_ode_withSwitch[debut:bout, 3]-data_corrected[indMinPeriod+dataLengthPeriod-1]
			
			# if surplus == True:
			#     fig = plt.figure(facecolor='w', figsize=(8, 4))
			#     ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
			#     ax.plot(data_surplus, label='data_surplus')
			#     plt.legend()
			#     #plt.show()
			#     plt.savefig(prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_surplus.png')
			#     #plt.savefig(prefFig + 'surplus.png')
			#     plt.close()

			# ajout des données dérivées
			# print('len(modelR1_deriv[indexplace, slicedataderiv])=', len(modelR1_deriv[indexplace, slicedataderiv]))
			# print('len(data_deriv    [indexplace, slicedataderiv])=', len(data_deriv    [indexplace, slicedataderiv]))
			data_deriv   [indexplace, slicedataderiv] = data_deriv_period
			modelR1_deriv[indexplace, slicedataderiv] = modelR1_deriv_period

			# ajout des SEIR1R2
			# print(np.shape(modelSEIR1R2[indexplace, slicedata, :]))
			# print(np.shape(sol_ode[indMinPeriod+sliceedo.start:indMinPeriod+sliceedo.stop, :]))
			# print(slicedata.start, slicedata.stop)
			# print(ts, ts+sliceedo.stop-sliceedo.start)
			# print(np.shape(modelSEIR1R2[indexplace, slicedata.start:slicedata.stop, :]))
			# print(np.shape(sol_ode))
			# print(np.shape(sol_ode[ts:ts+sliceedo.stop-sliceedo.start, :]))
			# input('pause')
			modelSEIR1R2[indexplace, slicedata.start:slicedata.stop, :] = sol_ode[ts:ts+sliceedo.stop-sliceedo.start, :]

			# preparation for next iteration
			_, E0, I0, R10, R20 = map(int, sol_ode[ts+dataLengthPeriod+recouvrement, :])

		Listepd.append(pd_exerpt)

		# calcul de l'EQM
		EQM = mean_squared_error(data_deriv[indexplace, :], modelR1_deriv[indexplace, :])
		ListeEQM.append(EQM)

		# udpate des listes pour transmission
		ListeTextParam.append(ListeTextParamPlace)
		ListetabParamModel.append(ListetabParamModelPlace)

	return modelSEIR1R2, ListeTextParam, Listepd, data_deriv, modelR1_deriv, decalage+recouvrement, ListetabParamModel, ListeEQM




def GetListDates(readStartDate, readStopDate, DatesString, decalage, nbperiodes, recouvrement):

	ListDates = [readStartDate, readStopDate]
	if nbperiodes!=1:
		confin_decalage = None
		if DatesString.listConfDates != []:
			confin_decalage = datetime.strptime(addDaystoStrDate(DatesString.listConfDates[0], decalage), strDate)
			if confin_decalage>ListDates[-2] and confin_decalage<ListDates[-1]:
				ListDates.insert(-1, confin_decalage)
		
		deconfin_decalage = None
		if DatesString.listDeconfDates != []:
			deconfin_decalage = datetime.strptime(addDaystoStrDate(DatesString.listDeconfDates[0], decalage), strDate)
			if deconfin_decalage>ListDates[-2] and deconfin_decalage<ListDates[-1]:
				ListDates.insert(-1, deconfin_decalage)

	ListePairDates = []
	for i in range(len(ListDates)-1):
		if i>0:
			ListePairDates.append((ListDates[i]+timedelta(recouvrement), ListDates[i+1]))
		else:
			ListePairDates.append((ListDates[i], ListDates[i+1]))

	# Conversion date chaine
	ListePairDatesStr = [(date1.strftime(strDate), date2.strftime(strDate)) for date1, date2 in ListePairDates]
	
	return ListePairDates, ListePairDatesStr


if __name__ == '__main__':
	modelSEIR1R2, ListeTextParam, Listepd, data_deriv, modelR1_deriv, decalage_corrige, ListetabParamModel, ListeEQM = fit(sys.argv[1:])
