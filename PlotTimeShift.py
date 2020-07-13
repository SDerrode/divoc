#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy             as np
import matplotlib.pyplot as plt
import seaborn           as sns

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

def PlotTimeShift(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python PlotTimeShift.py 
		>> python PlotTimeShift.py France      0 SEIR1R2  0 18,19 1 1
		>> python PlotTimeShift.py Italy,Spain 1 SEIR1R2D 1 18,19 0 1 # Italy and Spain, with UKF filtering

		For French Region (French database)
		>> python PlotTimeShift.py FRANCE,D69         0 SEIR1R2  0 18,19 0 1 # Code Insee Dpt 69 (Rhône)
		>> python PlotTimeShift.py FRANCE,R84         0 SEIR1R2  0 18,19 0 1 # Tous les dpts de la Région dont le code Insee est  regions)
		>> python PlotTimeShift.py FRANCE,R32+        0 SEIR1R2  0 18,19 0 1 # Somme de tous les dpts de la Région 32 (Hauts de F French regions)
		>> python PlotTimeShift.py FRANCE,MetropoleD  0 SEIR1R2  0 18,19 0 1 # Tous les départements de la France métropolitaine
		>> python PlotTimeShift.py FRANCE,MetropoleD+ 0 SEIR1R2  0 18,19 0 1 # Toute la France métropolitaine (en sommant les dpts)
		>> python PlotTimeShift.py FRANCE,MetropoleR+ 0 SEIR1R2  0 18,19 0 1 # Somme des dpts de toutes les régions françaises
		Toute combinaison de lieux est possible : exemple FRANCE,R32+,D05,R84
		
		argv[1] : List of countries (ex. France,Germany,Italy), or see above.  Default: France 
		argv[2] : Sex (male:1, female:2, male+female:0). Only for french database     Default: 0 
		argv[3] : EDO model (SEIR1R2 or SEIR1R2D).                             Default: SEIR2R2         
		argv[4] : UKF filtering of data (0/1).                                 Default: 0
		argv[5] : min shift, max shift, ex. 2,10.                              Default: 18,19
		argv[6] : Verbose level (debug: 3, ..., almost mute: 0).               Default: 1
		argv[7] : Plot graphique (0/1).                                        Default: 1
		argv[8] : stopDate.                                                    Default: None
	"""

	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	# Il y a 18 pays
	
	# Interpretation of arguments - reparation
	######################################################@


	SMALL_SIZE = 16
	MEDIUM_SIZE = 20
	BIGGER_SIZE = 24

	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

	if len(sysargv) > 9:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	places                 = 'France'
	sexe, sexestr          = 0, 'male+female'
	listplaces             = list(places.split(','))
	modeleString           = 'SEIR1R2'
	UKF_filt, UKF_filt01   = False, 0  #True, 1
	shift_mini, shift_maxi = 18,19     #4,18
	verbose                = 1
	plot                   = True
	stopDate               = "2020-07-01"
	
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
	TAB_IEnd = []
	TAB_ListeEQM    = []
	TAB_ListeDateI0 = []

	for decalage in range(shift_mini, shift_maxi):

		if verbose>0:
			print('TIME-SHIFT', str(decalage), 'OVER', str(shift_maxi))
	
		_, _, _, _, _, tabParamModel, tabIEnd, ListeEQM, ListeDateI0 = fit([places, sexe, nbperiodes, decalage, UKF_filt, 0, 0, stopDate])

		TAB_decalage.append(float(decalage))
		TAB_param_model.append(tabParamModel)
		TAB_IEnd.append(tabIEnd)
		TAB_ListeEQM.append(ListeEQM)
		TAB_ListeDateI0.append(ListeDateI0)

	#TAB_param_model[decalage][place][nbperiodes]
	#input('apuse')

	X = np.linspace(shift_mini, shift_maxi-1, shift_maxi-shift_mini)
	if modeleString == 'SEIR1R2':
		labelsparam  = [r'$a$', r'$b$', r'$c$', r'$f$', r'$R_0$']
	else:
		labelsparam  = [r'$a$', r'$b$', r'$c$', r'$f$', r'$\mu$', r'$\xi$', r'$R_0$']
	labelsperiod = ['Period 1', 'Period 2', 'Period 3']

	# On enregistre le R0 moyen de la 3ieme période pour faire carte graphique
	rep  = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/TimeShift/', './figures/'+modeleString+'/TimeShift/')
	fileR0moyen = rep+'/R0Moyen_' + str(shift_mini)+ '_' + str(shift_maxi) +'.csv'
	if os.path.exists(fileR0moyen):
		os.remove(fileR0moyen)
	with open(fileR0moyen, 'a') as text_file:
		text_file.write('Place,R0MoyenP0,R0MoyenP1,R0MoyenP2,DateFirstCase,IEndP0,IEndP1,IEndP2\n')


	ListeChaines = []
	for indexplace, place in enumerate(listplaces):

		# Get the full name of the place to process, and the special dates corresponding to the place
		if FrDatabase == True: 
			if 'MetropoleD+' in listnames[indexplace][0]:
				placefull = 'France'
			else:
				placefull  = 'France-' + listnames[indexplace][0]
		else:
			placefull  = place

		# Repertoire des figures
		if plot==True:
			ch1 = placefull+'/sexe_'+str(sexe)+'_delay_'+str(shift_mini)+'_'+str(shift_maxi)
			repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/'+ch1, './figures/'+modeleString+'/'+ch1)
			prefFig    = repertoire+'/'

		nbperiodes = len(TAB_param_model[0][indexplace][:])

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
				if sexe==0:
					titre = placefull + ' - ' + modeleString2 + ' parameters evolution for ' + labelsperiod[period]
				else:
					titre = placefull + ' - Sex=' + sexestr + ', ' + modeleString2 + ' parameters evolution for ' + labelsperiod[period]
				texte = list(map( lambda s: s.replace('$', '').replace('\\', '').replace('_', ''), labelsparam[:-1]))
				filename = prefFig   + 'Plot_TS_Period' + str(period) + '_' + ''.join(texte) + '.png'
				plotData(TAB_decalage, Y1[:, :-1], titre, filename, labelsparam[:-1])
				filename = prefFig   + 'Plot_TS_Period' + str(period) + '_R0.png'
				plotData(TAB_decalage, Y1[:, -1].reshape(shift_maxi-shift_mini, 1), titre, filename, [labelsparam[-1]])


		# plot pour les paramètres
		##########################################
		if plot==True:
			if os.path.exists(prefFig+'Plot_TS.txt'):
				os.remove(prefFig+'Plot_TS.txt')

		Y2 = np.zeros(shape=(shift_maxi-shift_mini, len(labelsperiod)))
		for param in range(len(labelsparam)):

			for decalage in range(shift_maxi-shift_mini):
				for period in range(len(labelsperiod)):
					try:
						Y2[decalage, period] = np.round(TAB_param_model[decalage][indexplace][period][param], 3)
					except IndexError:
						Y2[decalage, period] = 0.
			if plot==True:
				if sexe==0:
					titre = placefull + ' - ' + modeleString2 + ' periods evolution for param ' + labelsparam[param]
				else:
					titre = placefull + ' - Sex=' + sexestr + ', ' + modeleString2 + ' periods evolution for param ' + labelsparam[param]
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
					# print('placefull=', placefull)
					# print('place=', place)
					# print('listnames=', listnames)
					#input('apuse')
					#Lieu = ''.join(filter(str.isdigit, placefull))
					Lieu = placefull
					if placefull[0]=='D' or placefull[0]=='R':
						Lieu = Lieu[1:]
					if placefull[-1]=='+':
						Lieu = Lieu[:-1]
					if Lieu[-1]=='D':
						Lieu = Lieu[:-1]
					chaine = Lieu+','
					# Les R0 moyens pour les 3 périodes
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
						chaine += str(R0Est) + ','
					
					# La date du premier infecté
					#Si -1. pour la periode 0, alors pas de date
					if -1. in Y2[:, 0]:
						chaine += 'Invalid'
					else:
						chaine += TAB_ListeDateI0[Indice][indexplace]

					# Les infectés en fin de période pour les 3 périodes
					for period in range(len(labelsperiod)):
						chaine += ',' + str(TAB_IEnd[decalage][indexplace][period])
					
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
			ax.plot(TAB_decalage, Y3, alpha=1.0, lw=2, label='MSE for ' + f'R\N{SUPERSCRIPT ONE}')

			ax.set_xlabel('Delay (delta) in days')
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
			if sexe==0:
				titre = placefull + ' - ' + modeleString2# + ', EQM on ' + f'R\N{SUPERSCRIPT ONE}'
			else:
				titre = placefull + ' - Sex=' + sexestr + ', ' + modeleString2# + ', EQM on ' + f'R\N{SUPERSCRIPT ONE}'
			plt.tight_layout(rect=(0, 0.03, 1., 0.95))
			plt.title(titre)
			plt.savefig(prefFig + 'Plot_TS_EQM_cumul.png', dpi=dpi)
			plt.close()

	# Plot des distributions des paramètres
	##########################################

	if plot==True and len(listplaces)>1:
		rep = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/TimeShift/', './figures/'+modeleString+'/TimeShift/')
		file = rep+'/Distrib_'

		for i, param in enumerate(labelsparam):
			# print('Param:', param, ', i=', i)
			xlim = True
			if i==len(labelsparam)-1:
				xlim = False
			for decalage in range(shift_maxi-shift_mini):
				Y3 = np.zeros(shape=(len(labelsperiod), len(listplaces)))
				for period in range(len(labelsperiod)):
					for indexplace, place in enumerate(listplaces):
						Y3[period, indexplace] = TAB_param_model[decalage][indexplace][period][i]

				titre = 'Distribution of parameter ' +  param + ' for France departements'
				filename = file + str(decalage+shift_mini)+'_param'+param
				PlotDistribParam(Y3, titre, filename, labelsperiod, xlim)

	return ListeChaines


def PlotDistribParam(Y3, titre, filename, label, xlim=False):

	fig = plt.figure(facecolor='w', figsize=figsize)
	ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
	# pour toutes les périodes sur le même dessin
	for i,l in enumerate(label):
		sns.distplot(Y3[i, :], kde=True, rug=True, hist=False, label=l, ax=ax, kde_kws={'linewidth': 1.5})

	if xlim==True:
		ax.set_xlim([-0.05, 1.05])
	ax.set_ylabel('')
	plt.title(titre)
	plt.tight_layout()
	plt.legend()
	plt.savefig(filename, dpi=dpi)
	plt.close()

	return 0

def plotData(X, Y, titre, filename, labels):
	
	fig = plt.figure(facecolor='w', figsize=figsize)
	ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
	
	for i in range(len(labels)):
		ax.plot(X, Y[:, i], alpha=1.0, lw=2, label=labels[i])

	ax.set_xlabel('Delay (delta) in days')
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
	#plt.gcf().subplots_adjust(bottom=0.15)
	plt.tight_layout(rect=(0, 0.03, 1., 0.95))
	plt.savefig(filename, dpi=dpi)
	plt.close()


if __name__ == '__main__':
	ListeChaines = PlotTimeShift(sys.argv)


