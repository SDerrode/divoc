#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
from datetime          import datetime, timedelta
  
from common            import getDates, addDaystoStrDate, get_WE_indice, drawAnnotation
from common            import getLowerDateFromString, getNbDaysBetweenDateFromString, getRepertoire
from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2
from SolveEDO_SEIR1R2D import SolveEDO_SEIR1R2D
from ProcessSEIR1R2D   import fit as fitProcessSEIR1R2D
from ProcessSEIR1R2    import fit as fitProcessSEIR1R2


dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)


def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python3 PlotTimeShift.py 
		>> python3 PlotTimeShift.py France      SEIR1R2  0 2,14 1 1
		>> python3 PlotTimeShift.py Italy,Spain SEIR1R2D 1 2,14 0 1 # Italy and Spain, with UKF filtering

		For French Region (French database)
		>> python3 PlotTimeShift.py France,69    SEIR1R2  0 2,14 0 1 # Dpt 69 (Rhône)
		>> python3 PlotTimeShift.py France,69,01 SEIR1R2D 1 2,14 1 1 # Dpt 69 (Rhône) + Dpt 01 (Ain) with UKF filtering
		
		argv[1] : Country (or list separeted by ','), or 'France' followed by a list of dpts. Default: France 
		argv[2] : EDO model (SEIR1R2 or SEIR1R2D)                                             Default: SEIR2R2         
		argv[3] : UKF filtering of data (0/1).                                                Default: 0
		argv[4] : min shift, max shift, ex. 2,10                                              Default: 2,10
		argv[5] : Verbose level (debug: 3, ..., almost mute: 0).                              Default: 1
		argv[6] : Plot graphique (0/1).                                                       Default: 1
	"""

	# Interpetation of arguments - reparation
	######################################################@

	print('Command line : ', sysargv, flush=True)
	if len(sysargv) > 7:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	places                 = 'France'
	listplaces             = list(places.split(','))
	modeleString           = 'SEIR1R2'
	UKF_filt,   UKF_filt01 = False, 0  #True, 1
	shift_mini, shift_maxi = 2, 10     #4, 18
	verbose                = 1
	plot                   = True
	
	# Parameters from argv
	if len(sysargv)>1: places, listplaces = sysargv[1], list(sysargv[1].split(','))
	FrDatabase = False
	if listplaces[0]=='France' and len(listplaces)>1:
		try:
			int(listplaces[1]) #If this is a number, then it is a french dpt
		except Exception as e:
			FrDatabase=False
		else:
			FrDatabase = True
			listplaces = listplaces[1:]

	if len(sysargv)>2: modeleString  = sysargv[2]
	if len(sysargv)>3 and int(sysargv[3])==1: UKF_filt, UKF_filt01 = True, 1
	if len(sysargv)>4: shift_mini, shift_maxi = map(int, sysargv[4].split(','))
	if len(sysargv)>5: verbose       = int(sysargv[5])
	if len(sysargv)>6 and int(sysargv[6])==0: plot     = False

	# le modèle à traiter (SEIR1R2 or SEIR1R2D)
	if modeleString == 'SEIR1R2':
		fit = fitProcessSEIR1R2
	elif modeleString == 'SEIR1R2D':
		fit = fitProcessSEIR1R2D
	else:
		print('Wrong EDO model, only SEIR1R2 or SEIR1R2D available!')
		exit(1)


	# fit avec 3 périodes + décalage
	##################################
	nbperiodes = -1
	
	TAB_decalage_corrige = []
	TAB_param_model      = []
	TAB_ListeEQM         = []

	for decalage in range(shift_mini, shift_maxi):

		print('TIME-SHIFT', str(decalage), 'OVER', str(shift_maxi))
	
		_, _, _, _, _, decalage_corrige, tabParamModel, ListeEQM = fit([places, nbperiodes, decalage, UKF_filt, 0, 0])

		TAB_decalage_corrige.append(float(decalage_corrige))
		TAB_param_model.append(tabParamModel)
		TAB_ListeEQM.append(ListeEQM)

	#TAB_param_model[decalage][place][nbperiodes]
	
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			placefull   = 'France-' + place
		else:
			placefull   = place

		# Repertoire des figures
		repertoire = getRepertoire(UKF_filt, './figures/' + modeleString + '_UKFilt/'+placefull, './figures/' + modeleString + '/' + placefull)
		prefFig    = repertoire + '/'

		nbperiodes = len(TAB_param_model[0][indexplace][:])
		X = np.linspace(shift_mini, shift_maxi-1, shift_maxi-shift_mini)
		if modeleString == 'SEIR1R2':
			labelsparam  = [r'$a$', r'$b$', r'$c$', r'$f$', r'$R_O$']
		else:
			labelsparam  = [r'$a$', r'$b$', r'$c$', r'$f$', r'$\mu$', r'$\xi$', r'$R_O$']
		labelsperiod = ['Period 1', 'Period 2', 'Period 3']

		# plot pour les 3 périodes
		##########################################@
		Y1 = np.zeros(shape=(shift_maxi-shift_mini, len(labelsparam)))
		for period in range(nbperiodes):

			for decalage in range(shift_maxi-shift_mini):
				try:
					Y1[decalage, :] = TAB_param_model[decalage][indexplace][period][:]
				except IndexError:
					Y1[decalage, :] = 0.

			texte = list(map( lambda s: s.replace('$', '').replace('\\', '').replace('_', ''), labelsparam[:-1]))
			titre    = placefull + ' - ' + modeleString + ' parameters evolution for ' + labelsperiod[period]
			filename = prefFig   + 'ParamEvol_Period' + str(period) + '_' + ''.join(texte) + '.png'
			plotData(TAB_decalage_corrige, Y1[:, :-1], titre, filename, labelsparam[:-1])
			titre    = placefull + ' - ' + modeleString + ' parameters evolution for ' + labelsperiod[period]
			filename = prefFig   + 'ParamEvol_Period' + str(period) + '_R0.png'
			plotData(TAB_decalage_corrige, Y1[:, -1].reshape(shift_maxi-shift_mini, 1), titre, filename, [labelsparam[-1]])

		# plot pour les paramètres
		##########################################
		Y2 = np.zeros(shape=(shift_maxi-shift_mini, len(labelsperiod)))
		for param in range(len(labelsparam)):

			for decalage in range(shift_maxi-shift_mini):
				for period in range(len(labelsperiod)):
					try:
						Y2[decalage, period] = TAB_param_model[decalage][indexplace][period][param]
					except IndexError:
						Y2[decalage, period] = 0.
			titre    = placefull + ' - ' + modeleString + ' periods evolution for param ' + labelsparam[param]
			filename = prefFig   + 'PeriodEvol_Param' + labelsparam[param].replace('$', '') + '.png'
			plotData(TAB_decalage_corrige, Y2, titre, filename, labelsperiod)

		# Write parameters in a file
		if verbose>0:
			with open(prefFig+'Params.txt', 'w') as text_file:
				for param in range(len(labelsparam)):
					text_file.write('\n\nParam: %s' % labelsparam[param].replace('$', '').replace('\\', ''))
					for period in range(len(labelsperiod)):
						#text_file.write('\n  -->%s:\n' % labelsperiod[period])
						np.savetxt(text_file, Y2[:, period], delimiter=', ', newline=', ', fmt='%.4f', header='\n  -->'+labelsperiod[period]+': ')


		# plot de l'EQM
		##########################################
		fig = plt.figure(facecolor='w', figsize=figsize)
		ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
		
		Y3 = np.zeros(shape=(len(X)))
		for k in range(len(Y3)):
			try:
				Y3[k] = TAB_ListeEQM[k][indexplace]
			except IndexError:
				Y3[k] = 0.
		ax.plot(TAB_decalage_corrige, Y3, alpha=1.0, lw=2, label='EQM')

		ax.set_xlabel('Time shift (days)')
		ax.yaxis.set_tick_params(length=0)
		ax.xaxis.set_tick_params(length=0)
		ax.grid(b=True, which='major', c='w', lw=1, ls='-')
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		
		legend = ax.legend()
		legend.get_frame().set_alpha(0.5)
		for spine in ('top', 'right', 'bottom', 'left'):
			ax.spines[spine].set_visible(False)

		plt.xlim([TAB_decalage_corrige[0], TAB_decalage_corrige[-1]])
		#plt.ylim([0, 1.0])

		# ajout d'un text d'annotation
		plt.title(placefull + ' - ' + modeleString + ', EQM on the number of detected cases' )
		plt.savefig(prefFig + 'EQM_Deriv_' + modeleString + '.png', dpi=dpi)
		plt.close()


def plotData(X, Y, titre, filename, labels):
	
	fig = plt.figure(facecolor='w', figsize=figsize)
	ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
	
	for i in range(len(labels)):
		ax.plot(X, Y[:, i], alpha=1.0, lw=2, label=labels[i])

	ax.set_xlabel('Time shift (days)')
	ax.yaxis.set_tick_params(length=0)
	ax.xaxis.set_tick_params(length=0)
	ax.grid(b=True, which='major', c='w', lw=1, ls='-')
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	
	legend = ax.legend()
	legend.get_frame().set_alpha(0.5)
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)

	plt.xlim([X[0], X[-1]])
	# if np.max(np.max(Y)) < 1.:
	# 	plt.ylim([0, 1.])

	# ajout d'un text d'annotation
	plt.title(titre)
	plt.savefig(filename, dpi=dpi)
	plt.close()


if __name__ == '__main__':
	main(sys.argv)


