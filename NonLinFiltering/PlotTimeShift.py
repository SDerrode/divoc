#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
from datetime          import datetime, timedelta
  
from common            import getRepertoire
from France            import getPlace
from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2
from SolveEDO_SEIR1R2D import SolveEDO_SEIR1R2D
from ProcessSEIR1R2D   import fit as fitProcessSEIR1R2D
from ProcessSEIR1R2    import fit as fitProcessSEIR1R2

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

# blackstart    = 2.0

def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python PlotTimeShift.py 
		>> python PlotTimeShift.py France      SEIR1R2  0 11,12 1 1
		>> python PlotTimeShift.py Italy,Spain SEIR1R2D 1 11,12 0 1 # Italy and Spain, with UKF filtering

		For French Region (French database)
		>> python PlotTimeShift.py FRANCE,D69         SEIR1R2  0 11,12 0 1 # Code Insee Dpt 69 (Rhône)
		>> python PlotTimeShift.py FRANCE,R84         SEIR1R2  0 11,12 0 1 # Tous les dpts de la Région dont le code Insee est  regions)
		>> python PlotTimeShift.py FRANCE,R32+        SEIR1R2  0 11,12 0 1 # Somme de tous les dpts de la Région 32 (Hauts de F French regions)
		>> python PlotTimeShift.py FRANCE,MetropoleD  SEIR1R2  0 11,12 0 1 # Tous les départements de la France métropolitaies
		>> python PlotTimeShift.py FRANCE,MetropoleR+ SEIR1R2  0 11,12 0 1 # Somme des dpts de toutes les régions françaises
		Toute combinaison est possible : exemple FRANCE,R32+,D05,R84
		
		argv[1] : List of countries (ex. France,Germany,Italy), or see above.  Default: France 
		argv[2] : EDO model (SEIR1R2 or SEIR1R2D).                             Default: SEIR2R2         
		argv[3] : UKF filtering of data (0/1).                                 Default: 0
		argv[4] : min shift, max shift, ex. 2,10.                              Default: 11,12
		argv[5] : Verbose level (debug: 3, ..., almost mute: 0).               Default: 1
		argv[6] : Plot graphique (0/1).                                        Default: 1
	"""

	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	# Il y a 18 pays
	#Région Auvergne Rhone-Alpes
	# Ain (01), Allier (03), Ardèche (07), Cantal (15), Drôme (26), Isère (38), Loire (42), Haute-Loire (43), Puy-de-Dôme (63), Rhône (69), Savoie (73), Haute-Savoie (74)
	# 01,03,07,15,26,38,42,43,63,69,73,74
	
	# Interpretation of arguments - reparation
	######################################################@

	if len(sysargv) > 7:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	places                 = 'France'
	listplaces             = list(places.split(','))
	modeleString           = 'SEIR1R2'
	UKF_filt, UKF_filt01   = False, 0  #True, 1
	shift_mini, shift_maxi = 11,12     #4,18
	verbose                = 1
	plot                   = True
	
	# Parameters from argv
	if len(sysargv)>1: places, liste = sysargv[1], list(sysargv[1].split(','))
	if len(sysargv)>2: modeleString = sysargv[2]
	if len(sysargv)>3 and int(sysargv[3])==1: UKF_filt, UKF_filt01 = True, 1
	if len(sysargv)>4: shift_mini, shift_maxi = map(int, sysargv[4].split(','))
	if len(sysargv)>5: verbose = int(sysargv[5])
	if len(sysargv)>6 and int(sysargv[6])==0: plot = False
	if shift_maxi-shift_mini==1:
		plot = False # Pas de plot possible s'il n'y a qu'une seule données

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
		places = [el[0] for el in listnames]
		places = 'FRANCE,'+','.join(places)
	else:
		listplaces = liste[:]
		FrDatabase = False

	# le modèle à traiter (SEIR1R2 or SEIR1R2D)
	if modeleString == 'SEIR1R2':
		fit = fitProcessSEIR1R2
	elif modeleString == 'SEIR1R2D':
		fit = fitProcessSEIR1R2D
	else:
		print('Wrong EDO model, only SEIR1R2 or SEIR1R2D available!')
		exit(1)

	if verbose>0:
		print('  Full command line : '+sysargv[0]+' '+places+' '+modeleString+' '+str(UKF_filt)+' '+str(shift_mini)+','+str(shift_maxi)+' '+str(verbose)+' '+str(plot), flush=True)
	

	# fit avec 3 périodes + décalage
	##################################
	nbperiodes = -1
	
	TAB_decalage    = []
	TAB_param_model = []
	TAB_ListeEQM    = []

	for decalage in range(shift_mini, shift_maxi):

		print('TIME-SHIFT', str(decalage), 'OVER', str(shift_maxi))
	
		_, _, _, _, _, tabParamModel, ListeEQM, ListeDateI0 = fit([places, nbperiodes, decalage, UKF_filt, 0, 0])

		TAB_decalage.append(float(decalage))
		TAB_param_model.append(tabParamModel)
		TAB_ListeEQM.append(ListeEQM)

	#TAB_param_model[decalage][place][nbperiodes]

	# On enregistre le R0 moyen de la 3ieme période pour faire carte graphique
	rep  = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/', './figures/'+modeleString+'/')
	fileR0moyen = rep+'/R0Moyen_' + str(shift_mini)+ '_' + str(shift_maxi) +'.csv'
	if os.path.exists(fileR0moyen):
		os.remove(fileR0moyen)
	with open(fileR0moyen, 'a') as text_file:
		text_file.write('Place,R0MoyenP0,R0MoyenP1,R0MoyenP2,DateFirstCase\n')
	
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			placefull   = 'France-' + listnames[indexplace][0]
		else:
			placefull   = place

		# Repertoire des figures
		repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/'+placefull, './figures/'+modeleString+'/'+placefull)
		prefFig    = repertoire+'/'

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

			if plot==True:
				texte = list(map( lambda s: s.replace('$', '').replace('\\', '').replace('_', ''), labelsparam[:-1]))
				titre    = placefull + ' - ' + modeleString + ' parameters evolution for ' + labelsperiod[period]
				filename = prefFig   + 'EvolParam_Period' + str(period) + '_' + ''.join(texte) + '.png'
				plotData(TAB_decalage, Y1[:, :-1], titre, filename, labelsparam[:-1])
				titre    = placefull + ' - ' + modeleString + ' parameters evolution for ' + labelsperiod[period]
				filename = prefFig   + 'EvolParam_Period' + str(period) + '_R0.png'
				plotData(TAB_decalage, Y1[:, -1].reshape(shift_maxi-shift_mini, 1), titre, filename, [labelsparam[-1]])


		# plot pour les paramètres
		##########################################
		if os.path.exists(prefFig+'EvolParams.txt'):
			os.remove(prefFig+'EvolParams.txt')

		Y2 = np.zeros(shape=(shift_maxi-shift_mini, len(labelsperiod)))
		for param in range(len(labelsparam)):

			for decalage in range(shift_maxi-shift_mini):
				for period in range(len(labelsperiod)):
					try:
						Y2[decalage, period] = TAB_param_model[decalage][indexplace][period][param]
					except IndexError:
						Y2[decalage, period] = 0.
			if plot==True:
				titre    = placefull + ' - ' + modeleString + ' periods evolution for param ' + labelsparam[param]
				filename = prefFig   + 'EvolPeriod_Param' + labelsparam[param].replace('$', '') + '.png'
				plotData(TAB_decalage, Y2, titre, filename, labelsperiod)

			# Write parameters in a file
			with open(prefFig+'EvolParams.txt', 'a') as text_file:
				text_file.write('\n\nParam: %s' % labelsparam[param].replace('$', '').replace('\\', ''))
				for period in range(len(labelsperiod)):
					#text_file.write('\n  -->%s:\n' % labelsperiod[period])
					np.savetxt(text_file, Y2[:, period], delimiter=', ', newline=', ', fmt='%.4f', header='\n  -->'+labelsperiod[period]+': ')
			if param==len(labelsparam)-1: # c'est à dire R0
				with open(fileR0moyen, 'a') as text_file:
					Lieu = ''.join(filter(str.isdigit, listnames[indexplace][0]))
					if 'A' in listnames[indexplace][0]: # Dpt de Corse
						Lieu = Lieu+'A' 
					if 'B' in listnames[indexplace][0]: # Dpt de Corse
						Lieu = Lieu+'B' 
					chaine = Lieu+','
					for period in range(len(labelsperiod)):
						R0mean = round(np.mean(Y2[:, period]), 2)
						chaine += str(R0mean)+','
					if round(np.mean(Y2[:, 0]), 2) != -1.:
						chaine += ListeDateI0[indexplace]+'\n'
					else:
						chaine += 'Invalid\n'
					text_file.write(chaine)
					print('chaine=', chaine) # moyenne uniquement pour la 3ieme période

		# plot de l'EQM
		##########################################
		if plot==True:

			fig = plt.figure(facecolor='w', figsize=figsize)
			ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
			
			Y3 = np.zeros(shape=(len(X)))
			for k in range(len(Y3)):
				try:
					Y3[k] = TAB_ListeEQM[k][indexplace]
				except IndexError:
					Y3[k] = 0.
			ax.plot(TAB_decalage, Y3, alpha=1.0, lw=2, label='EQM')

			ax.set_xlabel('Time shift (days)')
			ax.yaxis.set_tick_params(length=0)
			ax.xaxis.set_tick_params(length=0)
			ax.grid(b=True, which='major', c='w', lw=1, ls='-')
			ax.xaxis.set_major_locator(MaxNLocator(integer=True))
			
			legend = ax.legend()
			legend.get_frame().set_alpha(0.5)
			for spine in ('top', 'right', 'bottom', 'left'):
				ax.spines[spine].set_visible(False)

			plt.xlim([TAB_decalage[0], TAB_decalage[-1]])
			#plt.ylim([0, 1.0])

			# ajout d'un text d'annotation
			plt.title(placefull + ' - ' + modeleString + ', EQM on the number of detected cases' )
			plt.savefig(prefFig + 'EvolPeriod_EQM_Deriv.png', dpi=dpi)
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


