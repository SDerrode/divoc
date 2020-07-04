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

dpi     = 300    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

# blackstart    = 2.0

def PlotTimeShift(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python PlotTimeShift.py 
		>> python PlotTimeShift.py France      0 SEIR1R2  0 13,14 1 1
		>> python PlotTimeShift.py Italy,Spain 1 SEIR1R2D 1 13,14 0 1 # Italy and Spain, with UKF filterin1

		For French Region (French database)
		>> python PlotTimeShift.py FRANCE,D69         0 SEIR1R2  0 13,14 0 1 # Code Insee Dpt 69 (Rhône)
		>> python PlotTimeShift.py FRANCE,R84         0 SEIR1R2  0 13,14 0 1 # Tous les dpts de la Région dont le code Insee est  regions)
		>> python PlotTimeShift.py FRANCE,R32+        0 SEIR1R2  0 13,14 0 1 # Somme de tous les dpts de la Région 32 (Hauts de F French regions)
		>> python PlotTimeShift.py FRANCE,MetropoleD  0 SEIR1R2  0 13,14 0 1 # Tous les départements de la France métropolitaine
		>> python PlotTimeShift.py FRANCE,MetropoleD+ 0 SEIR1R2  0 13,14 0 1 # Toute la France métropolitaine (en sommant les dpts)
		>> python PlotTimeShift.py FRANCE,MetropoleR+ 0 SEIR1R2  0 13,14 0 1 # Somme des dpts de toutes les régions françaises
		Toute combinaison de lieux est possible : exemple FRANCE,R32+,D05,R84
		
		argv[1] : List of countries (ex. France,Germany,Italy), or see above.  Default: France 
		argv[2] : Sex (male:1, female:2, male+female:0). Only for french database     Default: 0 
		argv[3] : EDO model (SEIR1R2 or SEIR1R2D).                             Default: SEIR2R2         
		argv[4] : UKF filtering of data (0/1).                                 Default: 0
		argv[5] : min shift, max shift, ex. 2,10.                              Default: 13,14
		argv[6] : Verbose level (debug: 3, ..., almost mute: 0).               Default: 1
		argv[7] : Plot graphique (0/1).                                        Default: 1
		argv[8] : stopDate.                                                    Default: None
	"""

	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	# Il y a 18 pays
	#Région Auvergne Rhone-Alpes
	# Ain (01), Allier (03), Ardèche (07), Cantal (15), Drôme (26), Isère (38), Loire (42), Haute-Loire (43), Puy-de-Dôme (63), Rhône (69), Savoie (73), Haute-Savoie (74)
	# 01,03,07,15,26,38,42,43,63,69,73,74
	
	# Interpretation of arguments - reparation
	######################################################@

	if len(sysargv) > 9:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	places                 = 'France'
	sexe, sexestr          = 0, 'male+female'
	listplaces             = list(places.split(','))
	modeleString           = 'SEIR1R2'
	UKF_filt, UKF_filt01   = False, 0  #True, 1
	shift_mini, shift_maxi = 13,14     #4,18
	verbose                = 1
	plot                   = True
	stopDate               = None
	
	# Parameters from argv
	if len(sysargv)>1: places, liste = sysargv[1], list(sysargv[1].split(','))
	if len(sysargv)>2: sexe = int(sysargv[2])
	if len(sysargv)>3: modeleString = sysargv[3]
	if len(sysargv)>4 and int(sysargv[4])==1: UKF_filt, UKF_filt01 = True, 1
	if len(sysargv)>5: shift_mini, shift_maxi = map(int, sysargv[5].split(','))
	if len(sysargv)>6: verbose = int(sysargv[6])
	if len(sysargv)>7 and int(sysargv[7])==0: plot = False
	if len(sysargv)>8: stopDate = sysargv[8]
	if shift_maxi-shift_mini==1:
		plot = False # Pas de plot possible s'il n'y a qu'une seule données
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
		places = [el[0] for el in listnames]
		places = 'FRANCE,'+','.join(places)
	else:
		listplaces = liste[:]
		FrDatabase = False

	# le modèle à traiter (SEIR1R2 or SEIR1R2D)
	if modeleString == 'SEIR1R2':
		fit           = fitProcessSEIR1R2
		modeleString2 =f'SEIR\N{SUPERSCRIPT ONE}R\N{SUPERSCRIPT TWO}'
	elif modeleString == 'SEIR1R2D':
		fit = fitProcessSEIR1R2D
		modeleString2 =f'SEIR\N{SUPERSCRIPT ONE}R\N{SUPERSCRIPT TWO}D'
	else:
		print('Wrong EDO model, only SEIR1R2 or SEIR1R2D available!')
		exit(1)

	if verbose>0:
		print('  Full command line : '+sysargv[0]+' '+places+' '+str(sexe)+' '+modeleString+' '+str(UKF_filt)+' '+str(shift_mini)+','+str(shift_maxi)+' '+str(verbose)+' '+str(plot), flush=True)
	

	# fit avec 3 périodes + décalage
	##################################
	nbperiodes = -1
	
	TAB_decalage    = []
	TAB_param_model = []
	TAB_ListeEQM    = []
	TAB_ListeDateI0 = []

	for decalage in range(shift_mini, shift_maxi):

		if verbose>0:
			print('TIME-SHIFT', str(decalage), 'OVER', str(shift_maxi))
	
		_, _, _, _, _, tabParamModel, ListeEQM, ListeDateI0 = fit([places, sexe, nbperiodes, decalage, UKF_filt, 0, 0, stopDate])

		TAB_decalage.append(float(decalage))
		TAB_param_model.append(tabParamModel)
		TAB_ListeEQM.append(ListeEQM)
		TAB_ListeDateI0.append(ListeDateI0)

	#TAB_param_model[decalage][place][nbperiodes]
	#input('apuse')

	# On enregistre le R0 moyen de la 3ieme période pour faire carte graphique
	rep  = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/TimeShift/', './figures/'+modeleString+'/TimeShift/')
	fileR0moyen = rep+'/R0Moyen_' + str(shift_mini)+ '_' + str(shift_maxi) +'.csv'
	if os.path.exists(fileR0moyen):
		os.remove(fileR0moyen)
	with open(fileR0moyen, 'a') as text_file:
		text_file.write('Place,R0MoyenP0,R0MoyenP1,R0MoyenP2,DateFirstCase\n')
	
	ListeChaines = []
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			placefull   = 'France-' + listnames[indexplace][0]
			placefull1   = listnames[indexplace][0]
		else:
			placefull   = place

		# Repertoire des figures
		ch1 = placefull+'/sexe_'+str(sexe)+'_shift_'+str(shift_mini)+'_'+str(shift_maxi)
		repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/'+ch1, './figures/'+modeleString+'/'+ch1)
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
				titre    = placefull + ', Sex:' + sexestr + ' - ' + modeleString2 + ' parameters evolution for ' + labelsperiod[period]
				filename = prefFig   + 'Plot_TS' + str(period) + '_' + ''.join(texte) + '.png'
				plotData(TAB_decalage, Y1[:, :-1], titre, filename, labelsparam[:-1])
				titre    = placefull + ', Sex:' + sexestr + ' - ' + modeleString2 + ' parameters evolution for ' + labelsperiod[period]
				filename = prefFig   + 'Plot_TS' + str(period) + '_R0.png'
				plotData(TAB_decalage, Y1[:, -1].reshape(shift_maxi-shift_mini, 1), titre, filename, [labelsparam[-1]])


		# plot pour les paramètres
		##########################################
		if os.path.exists(prefFig+'Plot_TS.txt'):
			os.remove(prefFig+'Plot_TS.txt')

		Y2 = np.zeros(shape=(shift_maxi-shift_mini, len(labelsperiod)))
		for param in range(len(labelsparam)):

			for decalage in range(shift_maxi-shift_mini):
				for period in range(len(labelsperiod)):
					try:
						Y2[decalage, period] = np.round(TAB_param_model[decalage][indexplace][period][param],3)
					except IndexError:
						Y2[decalage, period] = 0.
			if plot==True:
				titre    = placefull + ', Sex:' + sexestr + ' - ' + modeleString2 + ' periods evolution for param ' + labelsparam[param]
				filename = prefFig   + 'Plot_TS_Param' + labelsparam[param].replace('$', '') + '.png'
				plotData(TAB_decalage, Y2, titre, filename, labelsperiod)

			# Write parameters in a file
			with open(prefFig+'Plot_TS.txt', 'a') as text_file:
				text_file.write('\n\nParam: %s' % labelsparam[param].replace('$', '').replace('\\', ''))
				for period in range(len(labelsperiod)):
					#text_file.write('\n  -->%s:\n' % labelsperiod[period])
					np.savetxt(text_file, Y2[:, period], delimiter=', ', newline=', ', fmt='%.4f', header='\n  -->'+labelsperiod[period]+': ')
			
			if param==len(labelsparam)-1: # c'est à dire R0
				with open(fileR0moyen, 'a') as text_file:
					# print('placefull1=', placefull1)
					# print('place=', place)
					# print('listnames=', listnames)
					#input('apuse')
					#Lieu = ''.join(filter(str.isdigit, placefull1))
					Lieu = placefull1
					if placefull1[0]=='D' or placefull1[0]=='R':
						Lieu = Lieu[1:]
					if placefull1[-1]=='+':
						Lieu = Lieu[:-1]
					if Lieu[-1]=='D':
						Lieu = Lieu[:-1]
					chaine = Lieu+','
					for period in range(len(labelsperiod)):
						if -1. in Y2[:, period]:
							R0Est = -1.
						else:
							#R0Est = np.mean(Y2[:, period])
							# Valeur médiane
							R0Est = sorted(Y2[:, period])[int((len(Y2[:, period])-1)/2)]
							if period==0:
								#print(Y2[:, period], R0Est)
								Indice = np.where(Y2[:, period]==R0Est)[0][0]
								# print('Indice=', Indice)
								# input('apuse')
						chaine += str(R0Est)+','
					
					#Si -1. pour la periode 0, alors pas de date
					if -1. in Y2[:, 0]:
						chaine += 'Invalid'
					else:
						chaine += TAB_ListeDateI0[Indice][indexplace]
					text_file.write(chaine+'\n')
					if verbose>0:
						print('chaine=', chaine)
					ListeChaines.append(chaine)

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
			plt.title(placefull + ', Sex:' + sexestr + ' - ' + modeleString2 + ', EQM on ' + f'R\N{SUPERSCRIPT ONE}' )
			plt.savefig(prefFig + 'Plot_TS_EQM_Deriv.png', dpi=dpi)
			plt.close()

	return ListeChaines


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
	ListeChaines = PlotTimeShift(sys.argv)


