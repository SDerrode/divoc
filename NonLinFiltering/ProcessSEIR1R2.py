#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt
from datetime          import datetime, timedelta

from common            import readDataEurope, readDataFrance, getDates, Plot, addDaystoStrDate
from common            import getLowerDateFromString, getNbDaysBetweenDateFromString
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
		argv[2] : Number of periods (1, or any number which means automatic).                        Default: 1 
		argv[3] : Shift of periods (in days).                                                        Default: 0
		argv[4] : Compensation strategy (0/1).                                                       Default: 0
		argv[5] : Verbose level (debug: 3, ..., almost mute: 0).                                     Default: 1
		argv[6] : Plot graphique (0/1).                                                              Default: 1
	"""

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
			int(listplaces[1])
		except Exception as e:
			FrDatabase=False
		else:
			FrDatabase = True
			listplaces = listplaces[1:]
			France     = 'France'
	
	if len(sysargv)>1: nbperiodes    = int(sysargv[1])
	if len(sysargv)>2: decalage      = int(sysargv[2])
	if len(sysargv)>3 and int(sysargv[3])==1: surplus = True
	if len(sysargv)>4: verbose       = int(sysargv[4])
	if len(sysargv)>5 and int(sysargv[5])==0: plot = False
	if nbperiodes==1: decalage=0

	# constante
	fileLocalCopy    = False  # if we upload the file from the url (to get latest data) or from a local copy file
	readStartDateStr = "2020-03-01" 
	readStopDateStr  = None

	# Lecture des données pour connaitre le nombre de données
	if FrDatabase == True: 
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance('69', readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
	else:
		pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope('France', readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
	dataLength = pd_exerpt.shape[0]

	# For the derivatives
	data_deriv     = np.zeros(shape=(len(listplaces), dataLength))
	moldelR1_deriv = np.zeros(shape=(len(listplaces), dataLength))
	listepd        = []
	for indexplace, place in enumerate(listplaces):

		placefull = place
		# These are the date of confinement and deconfinement + other. 
		# See function getDates on how to add or delete dates to put the focus on
		DatesString = getDates(place, verbose)
		if FrDatabase == True: 
			placefull = France + '-' + place
			DatesString = getDates(France, verbose)

		print('PROCESSING of', placefull, 'in', listplaces)
		prefFig = './figures/Data_' + placefull

		# Lecture des données et copy of the observations
		#############################################################################
		if FrDatabase == True:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
		else:
			pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)

		readStartDate = datetime.strptime(readStartDateStr, "%Y-%m-%d")
		if readStartDate<pd_exerpt.index[0]:
			readStartDate    = pd_exerpt.index[0]
			readStartDateStr = pd_exerpt.index[0].strftime("%Y-%m-%d")
		readStopDate = datetime.strptime(readStopDateStr, "%Y-%m-%d")
		if readStopDate<pd_exerpt.index[-1]:
			readStopDate    = pd_exerpt.index[-1]
			readStopDateStr = pd_exerpt.index[-1].strftime("%Y-%m-%d")

		dataLength = pd_exerpt.shape[0]
		if verbose>1:
			print('readStartDateStr=', readStartDateStr, ', readStopDateStr=', readStopDateStr)
			print('readStartDate   =', readStartDate,    ', readStopDate   =', readStopDate)
			#input('pause')

		# Get the list od dates to process
		recouvrement = -1
		ListDates, ListDatesStr = GetListDates(readStartDate, readStopDate, DatesString, decalage, nbperiodes, recouvrement)
		if verbose>1:
			print('ListDates   =', ListDates)
			print('ListDatesStr=', ListDatesStr)
			#input('pause')
		
		# Modele d'edo non lineaires
		#############################################################################
		# Solveur edo
		dt      = 1
		solveur = SolveEDO_SEIR1R2(N, dt, verbose)

		# Constantes
		prefFig = './figures/' + solveur.modele.modelShortName + '_' + placefull
		columns = [r'$S(t)$', r'$E(t)$', r'$I(t)$', r'$R^1(t)$', r'$R^2(t)$', r'$R^2(t)=N-\sum(SEIR^1)$', r'$R(t)=R^1(t)+R^2(t)$']

		# Surplus de la précédente période (nulle pour le début)
		data_surplus   = np.zeros(shape=(dataLength))
		data_corrected = np.zeros(shape=(dataLength))

		# Boucle pour traiter successivement les différentes fenêtres
		#############################################################################
		for i in range(len(ListDatesStr)):

			# paramètres initiaux
			if i==0:
				E0, I0, R10, R20 = 0, 1, 0, 0
				# pour une seule période
				a0, b0, c0, f0   = 0.35, 0.29, 0.075, 0.0022
				T  = 350

				# pour plusieurs périodes
				if nbperiodes!=1:
					l, b0, c0, f0 = 0.255, 1./5.2, 1./12, 0.08
					a0 = (l+c0)*(1.+l/b0)
					T  = 150

			if i==1:
				_, E0, I0, R10, R20 = map(int, sol_ode[ts+dataLengthPeriod-1, :]) # Le -1 vient du fait qu'on fait un petit recouvrement

				if surplus==True:
					a0, b0, c0, f0   = 0.100, 0.18, 0.09, 0.14 # en enlevant le surplus
				else:
					a0, b0, c0, f0   = 0.012, 0.25, 0.073, 0.057 # sans enlever le surplus
				T = 120

			if i==2:
				_, E0, I0, R10, R20 = map(int, sol_ode[ts+dataLengthPeriod-1, :]) # Le -1 vient du fait qu'on fait un petit recouvrement
				_, a0, b0, c0, f0 = solveur.modele.getParam()
				T = 50
			

			time = np.linspace(0, T-1, T)

			# date 
			fitStartDate,   fitStopDate    = ListDates[i]
			fitStartDateStr,fitStopDateStr = ListDatesStr[i]

			# Récupérations des données observées
			dataLengthPeriod = 0
			indMinPeriod     = getNbDaysBetweenDateFromString(readStartDateStr, fitStartDateStr)
			for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), (HeadData[0])]):
				data_corrected[indMinPeriod+j] = z - data_surplus[indMinPeriod+j]
				dataLengthPeriod +=1
			timefocusdata      = slice(indMinPeriod, indMinPeriod+dataLengthPeriod)
			timefocusdataderiv = slice(timefocusdata.start+1, timefocusdata.stop)
			if verbose>1:
				print('  dataLength      =', dataLength)
				print('  indMinPeriod    =', indMinPeriod)
				print('  dataLengthPeriod=', dataLengthPeriod)
				print('  fitStartDateStr =', fitStartDateStr)
				print('  fitStopDateStr  =', fitStopDateStr)
				#input('attente')

			# Set initialisation data for the solveur
			solveur.modele.setParam(N=N, a=a0,  b=b0,  c=c0,    f=f0)
			solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20)

			# Before optimization
			###############################

			# Solve ode avant optimization
			sol_ode = solveur.solve_SEIR1R2(time)
			# calcul time shift initial (ts0) with respect to data
			ts0 = solveur.compute_tsfromEQM(data_corrected[timefocusdata], T)
			if verbose>0:
				print(solveur)
				print('  ts0='+str(ts0))
			
			# plot
			if plot==True:
				listePlot    = [3]
				filename     = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fitinit_' + ''.join(map(str, listePlot)) + '.png'
				timefocusedo = slice(ts0, min(ts0+dataLengthPeriod, T))
				solveur.plot_SEIR1R2(filename, placefull, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))

			# Parameters optimization
			############################################################################

			solveur.paramOptimization(data_corrected[timefocusdata], time)
			if verbose>0:
				print('Solver''s state after optimization=', solveur)
				_, a1, b1, c1, f1 = solveur.modele.getParam()
				print('  Reproductivité après: ', a1/c1)


			# After optimization
			###############################
			
			# Solve equa diff apres optimization
			sol_ode = solveur.solve_SEIR1R2(time)
			# calcul time shift (ts) with respect to data
			ts = solveur.compute_tsfromEQM(data_corrected[timefocusdata], T)
			if verbose>0:
				print(solveur)
				print('  ts='+str(ts))

			timefocusedo      = slice(ts, min(ts+dataLengthPeriod, T))
			timefocusedoderiv = slice(timefocusedo.start+1, timefocusedo.stop)
			data_deriv_period    = (data_corrected[timefocusdataderiv] - data_corrected[timefocusdataderiv.start-1:timefocusdataderiv.stop-1 ]) / dt
			modelR1_deriv_period = (sol_ode[timefocusedoderiv, 3]      - sol_ode       [timefocusedoderiv.start-1 :timefocusedoderiv.stop-1,3]) / dt
			# print('len(data_deriv_period)=', len(data_deriv_period))
			# print('len(modelR1_deriv_period)=', len(modelR1_deriv_period))
			# input('pause')

			if plot==True:
				titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - Shift=' + str(decalage)
				
				listePlot =[0,1,2,3,4]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, titre, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
				listePlot =[1,2,3]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, titre, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
				listePlot =[3]
				filename  = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
				solveur.plot_SEIR1R2(filename, titre, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))

				# dérivée  numérique de R1
				filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fit_3Deriv.png'
				solveur.plot_dR1(filename, titre, timefocusedoderiv, timefocusdataderiv, data_deriv_period, text=solveur.getTextParam(readStartDateStr))

			# sol_ode_withSwitch = solveur.solve_SEIR1R2_withSwitch(T, timeswitch=ts+dataLengthPeriod)
			# if plot==True:
			#     titre = placefull + '- Period ' + str(i) + '\\' + str(len(ListDatesStr)-1) + ' - Shift=' + str(decalage)
			#
			#     listePlot=[0,1,2,3,4]
			#     filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
			#     solveur.plot_SEIR1R2(filename, titre, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
			#     listePlot=[1,2,3]
			#     filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
			#     solveur.plot_SEIR1R2(filename, titre, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))
			#     listePlot=[3]
			#     filename = prefFig + '_Period' + str(i) + '_' + str(len(ListDatesStr)-1) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
			#     solveur.plot_SEIR1R2(filename, titre, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam(readStartDateStr))

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
			data_deriv    [indexplace, timefocusdata.start+1:timefocusdata.stop] = data_deriv_period
			moldelR1_deriv[indexplace, timefocusdata.start+1:timefocusdata.stop] = modelR1_deriv_period

		listepd.append(pd_exerpt)

	return listepd, data_deriv, moldelR1_deriv


def GetListDates(readStartDate, readStopDate, DatesString, decalage, nbperiodes, recouvrement):

	ListDates = [readStartDate, readStopDate]
	if nbperiodes!=1:
		confin_decalage = None
		if DatesString.listConfDates != []:
			confin_decalage = datetime.strptime(addDaystoStrDate(DatesString.listConfDates[0], decalage), "%Y-%m-%d")
			if confin_decalage>ListDates[-2] and confin_decalage<ListDates[-1]:
				ListDates.insert(-1, confin_decalage)
		
		deconfin_decalage = None
		if DatesString.listDeconfDates != []:
			deconfin_decalage = datetime.strptime(addDaystoStrDate(DatesString.listDeconfDates[0], decalage), "%Y-%m-%d")
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
	fit(sys.argv[1:])
