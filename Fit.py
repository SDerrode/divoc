#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from datetime          import datetime, timedelta

from common            import readDates, addDaystoStrDate, get_WE_indice, drawAnnotation
from common            import getLowerDateFromString, getNbDaysBetweenDateFromString, getRepertoire
from France            import getPlace
from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2
from SolveEDO_SEIR1R2D import SolveEDO_SEIR1R2D
from ProcessSEIR1R2D   import fit as fitProcessSEIR1R2D
from ProcessSEIR1R2    import fit as fitProcessSEIR1R2

dpi     = 300    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python Fit.py 
		>> python Fit.py France      0 SEIR1R2 18 0 1 1
		>> python Fit.py Italy,Spain 1 SEIR1R2D 18 1 1 1 # Italy and Spain, with UKF filtering

		For French Region (French database)
		>> python Fit.py FRANCE,D69         0 SEIR1R2  18 0 1 1 # Code Insee Dpt 69 (Rhône)
		>> python Fit.py FRANCE,R84         0 SEIR1R2  18 0 1 1 # Tous les dpts de la Région dont le code Insee est 
		>> python Fit.py FRANCE,R32+        0 SEIR1R2  18 0 1 1 # Somme de tous les dpts de la Région 32 (Hauts de F
		>> python Fit.py FRANCE,MetropoleD  0 SEIR1R2  18 0 1 1 # Tous les départements de la France métropolitaine
		>> python Fit.py FRANCE,MetropoleD+ 0 SEIR1R2  18 0 1 1 # Toute la France métropolitaine (en sommant les dpts)
		>> python Fit.py FRANCE,MetropoleR+ 0 SEIR1R2  18 0 1 1 # Somme des dpts de toutes les régions françaises
		Toute combinaison de lieux est possible : exemple FRANCE,R32+,D05,R84
		
		argv[1] : List of countries (ex. France,Germany,Italy), or see above.  Default: France 
		argv[2] : Sex (male:1, female:2, male+female:0). Only for french database     Default: 0 
		argv[3] : EDO model (SEIR1R2 or SEIR1R2D).                             Default: SEIR2R2         
		argv[4] : Delay (in days).                                             Default: 18
		argv[5] : UKF filtering of data (0/1).                                 Default: 0
		argv[6] : Verbose level (debug: 3, ..., almost mute: 0).               Default: 1
		argv[7] : Plot graphique (0/1).                                        Default: 1
	"""

	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	#Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
	# Il y a 18 pays

	# Interpetation of arguments - reparation
	######################################################@

	if len(sysargv) > 8:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	places        = 'France'
	sexe, sexestr = 0, 'male+female'
	listplaces    = list(places.split(','))
	modeleString  = 'SEIR1R2'
	decalage3P    = 18
	UKF_filt, UKF_filt01 = False, 0
	verbose       = 1
	plot          = True
	France        = 'France'

	# Parameters from argv
	if len(sysargv)>1: places, liste = sysargv[1], list(sysargv[1].split(','))
	if len(sysargv)>2: sexe = int(sysargv[2])
	if len(sysargv)>3: modeleString = sysargv[3]
	if len(sysargv)>4: decalage3P   = int(sysargv[4])
	if len(sysargv)>5 and int(sysargv[5])==1: UKF_filt, UKF_filt01 = True, 1
	if len(sysargv)>6: verbose      = int(sysargv[6])
	if len(sysargv)>7 and int(sysargv[7])==0: plot = False
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
		fit = fitProcessSEIR1R2
	elif modeleString == 'SEIR1R2D':
		fit = fitProcessSEIR1R2D
	else:
		print('Wrong EDO model, only SEIR1R2 or SEIR1R2D available!')
		exit(1)

	if verbose>0:
		print('  Full command line : '+sysargv[0]+' '+places+' '+str(sexe)+' '+modeleString+' '+str(decalage3P)+' '+str(UKF_filt)+' '+str(verbose)+' '+str(plot), flush=True)

	# fit avec 3 périodes + décalage
	##################################
	nbperiodes = -1

	model_piecewise, ListeTextParamPlace_piecewise, liste_pd_piecewise, data_deriv_piecewise, model_deriv_piecewise, _, _, _, ListeDateI0 = \
			fit([places, sexe, nbperiodes, decalage3P, UKF_filt01, 0, 0])

	ListeTestPlace = []
	for indexplace in range(len(listplaces)):
		texteplace = ''
		for texte in ListeTextParamPlace_piecewise[indexplace]:
			texteplace += '\n' + texte
		ListeTestPlace.append(texteplace)


	# Plot the multi-period strategy
	##################################
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
		repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/'+placefull+'/sexe_'+str(sexe)+'_delay_'+str(decalage3P), './figures/'+modeleString+'/'+placefull+'/sexe_'+str(sexe)+'_delay_'+str(decalage3P))
		prefFig    = repertoire + '/Fit_'

		# Preparation plot pandas
		listheader = list(liste_pd_piecewise[indexplace])

		if FrDatabase==True:
			DatesString = readDates(France, verbose)
		else:
			DatesString = readDates(place, verbose)

		
		#####################################################@
		# DERIVEES
		# on ajoute les dérivées numériques des cas et des morts
		liste_pd_piecewise[indexplace]['dcases']           = liste_pd_piecewise[indexplace][listheader[0]].diff()
		liste_pd_piecewise[indexplace]['ddeaths']          = liste_pd_piecewise[indexplace][listheader[1]].diff()
		liste_pd_piecewise[indexplace]['dcasesplusdeaths'] = liste_pd_piecewise[indexplace][listheader[2]].diff()
		longueur = len(liste_pd_piecewise[indexplace].loc[:, ('dcases')])

		# on ajoute les dérivées numériques des cas et des morts
		liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise')]          = model_deriv_piecewise[indexplace, 0:longueur, 0]
		liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise_residual')] = liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise')]-liste_pd_piecewise[indexplace].loc[:, ('dcasesplusdeaths')]
		if modeleString=='SEIR1R2D':
			liste_pd_piecewise[indexplace].loc[:, ('md_piecewise')]          = model_deriv_piecewise[indexplace, 0:longueur, 1]
			liste_pd_piecewise[indexplace].loc[:, ('md_piecewise_residual')] = liste_pd_piecewise[indexplace].loc[:, ('md_piecewise')]-liste_pd_piecewise[indexplace].loc[:, ('ddeaths')]

		# Dessin des dérivées
		filename = prefFig + str(decalage3P) + '_Diff_Piecewise.png'
		if sexe==0:
			title = placefull + ' - Delay (delta)=' + str(decalage3P) + ' day(s)'
		else:
			title = placefull + ' - Sex=' + sexestr + ', Delay (delta)=' + str(decalage3P) + ' day(s)'
		if modeleString=='SEIR1R2':
			listPlots = ['dcasesplusdeaths', 'mc_piecewise']
		if modeleString=='SEIR1R2D':
			listPlots = ['dcases', 'mc_piecewise', 'ddeaths', 'md_piecewise']
		PlotFitPiecewise(modeleString, liste_pd_piecewise[indexplace], title, filename, y=listPlots, Dates=DatesString, textannotation=ListeTestPlace[indexplace])

		# Dessin des résidus des dérivées
		filename = prefFig + str(decalage3P) + '_Diff_PiecewiseResiduals.png'
		if sexe==0:
			title = placefull + ' - Delay (delta)=' + str(decalage3P) + ' day(s)'
		else:
			title = placefull + ' - Sex=' + sexestr + ', Delay (delta)=' + str(decalage3P) + ' day(s)'

		if modeleString=='SEIR1R2':
			listPlots = ['mc_piecewise_residual']
		if modeleString=='SEIR1R2D':
			listPlots = ['mc_piecewise_residual', 'md_piecewise_residual']
		PlotFitPiecewiseResidual(modeleString, liste_pd_piecewise[indexplace], title, filename, y=listPlots, Dates=DatesString, textannotation=ListeTestPlace[indexplace])


	# Plot the three SEIR1R2
	##################################

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
		repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/'+placefull+'/sexe_'+str(sexe)+'_delay_'+str(decalage3P), './figures/'+modeleString+'/'+placefull+'/sexe_'+str(sexe)+'_delay_'+str(decalage3P))
		prefFig    = repertoire + '/Fit_'

		# Preparation plot pandas
		listheader = list(liste_pd_piecewise[indexplace])
		longueur = len(liste_pd_piecewise[indexplace].loc[:, (listheader[0])])
		
		if FrDatabase==True:
			DatesString = readDates(France, verbose)
		else:
			DatesString = readDates(place, verbose)

		liste_pd_piecewise[indexplace].loc[:, ('S(t)')]  = model_piecewise [indexplace, :, 0]
		liste_pd_piecewise[indexplace].loc[:, ('E(t)')]  = model_piecewise [indexplace, :, 1]
		liste_pd_piecewise[indexplace].loc[:, ('I(t)')]  = model_piecewise [indexplace, :, 2]
		liste_pd_piecewise[indexplace].loc[:, ('R1(t)')] = model_piecewise [indexplace, :, 3]
		liste_pd_piecewise[indexplace].loc[:, ('R2(t)')] = model_piecewise [indexplace, :, 4]
		if modeleString == 'SEIR1R2D':
			liste_pd_piecewise[indexplace].loc[:, ('D(t)')]  = model_piecewise [indexplace, :, 5]

		if sexe==0:
			titre = placefull + ' - Delay (delta)=' + str(decalage3P) + ' day(s)'
		else:
			titre = placefull + ' - Sex=' + sexestr + ', Delay (delta)=' + str(decalage3P) + ' day(s)'

		listePlot=['E(t)', 'I(t)', 'R1(t)']
		if modeleString=='SEIR1R2D':
			listePlot.append('D(t)')
		filename = prefFig + str(decalage3P) + '_' + ''.join(map(str, listePlot)) + '_piecewise.png'
		PlotModel(modeleString, liste_pd_piecewise[indexplace], titre, filename, y=listePlot, Dates=DatesString, textannotation=ListeTestPlace[indexplace])

		listePlot=['R1(t)']
		if modeleString=='SEIR1R2D':
			listePlot.append('D(t)')
		filename = prefFig + str(decalage3P) + '_' + ''.join(map(str, listePlot)) + '_piecewise.png'
		
		PlotModel(modeleString, liste_pd_piecewise[indexplace], titre, filename, y=listePlot, Dates=DatesString, textannotation=ListeTestPlace[indexplace])


def PlotModel(modeleString, pd, titre, filenameFig, y, Dates=None, textannotation=''):
	
	if len(y)==0 or y is None: pass

	# on enlev la premiere et la dernière valeurs
	pd = pd.iloc[1:-1]

	# juste pour récupérer des couleurs sohérentes
	if modeleString=='SEIR1R2':
		solveur = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
	else:
		solveur = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)

	mod = solveur.modele
	listeColor = []
	for p in y:
		if p=='S(t)' : listeColor.append(mod.getColor(0))
		if p=='E(t)' : listeColor.append(mod.getColor(1))
		if p=='I(t)' : listeColor.append(mod.getColor(2))
		if p=='R1(t)': listeColor.append(mod.getColor(3))
		if p=='R2(t)': listeColor.append(mod.getColor(4))
		if p=='D(t)' : listeColor.append(mod.getColor(5))

	fig = plt.figure(facecolor='w', figsize=figsize)
	ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

	# Dessin des courbes théoriques
	pd.plot(ax=ax, y=y, title=titre, color=listeColor)
	
	# ajout des dates spéciales
	if Dates!=None:
		# for d in Dates.listFirstCaseDates:
		#     drawAnnotation(ax, 'First case date\n', d, color='blue')
		for d in Dates.listConfDates:
			drawAnnotation(ax, 'Conf. date\n', d, color='red')
		for d in Dates.listDeconfDates:
			drawAnnotation(ax, 'Deconf. date\n', d, color='green')
		for d in Dates.listOtherDates:
			drawAnnotation(ax, 'Other date\n', d)

	# surlignage des jours de WE
	WE_indices = get_WE_indice(pd)
	i = 0
	while i < len(WE_indices)-1:
		ax.axvspan(pd.index[WE_indices[i]], pd.index[WE_indices[i+1]], facecolor='gray', edgecolor='none', alpha=.15, zorder=-100)
		i += 2

	# ajout d'un text d'annotation
	if textannotation != '':
		bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
		ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/3.), fontsize=6, bbox=bbox, ha="left", va="center") 

	# axes
	ax.grid(True, which='major', axis='both')
	ax.grid(True, which='minor', axis='both')
	ax.grid(True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
	ax.grid(True, which='minor', c='w', lw=0.5, ls='-')
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)

	plt.ticklabel_format(style='plain', axis='y')

	# On enlève le label sur l'axe x
	x_label = ax.axes.get_xaxis().get_label().set_visible(False)

	# legende
	legend = ax.legend().get_frame().set_alpha(0.8)
	plt.tight_layout(rect=(0, 0, 1., 0.95))
	plt.legend(fontsize=8)
	plt.tight_layout()
	plt.savefig(filenameFig, dpi=dpi)
	plt.close()

def PlotFitPiecewise(modeleString, pd, titre, filenameFig, y, Dates=None, textannotation=''):

	if len(y)==0 or y is None: pass

	# on enlev la premiere et la dernière valeurs
	pd = pd.iloc[1:-1]
	
	fig = plt.figure(facecolor='w', figsize=figsize)
	ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

	# Dessin des courbes
	if modeleString=='SEIR1R2':
		solveur    = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
		colori     = solveur.modele.getColor(3)
		linestyles = ['-', '-']#, '-']
		markers    = ['x', '']#, '']
		colors     = [colori, colori]#, 'blue']
		linewidths = [0.5, 1.5]#, 1.5]
		alphas     = [1.0, 0.7]#, 0.55]
		labels     = [r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise']
	else:
		solveur    = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)
		colori     = [solveur.modele.getColor(3), solveur.modele.getColor(5)]
		linestyles = ['-', '-', '-', '-']#, '-']
		markers    = ['x', '', '+', '']#, '']
		colors     = [colori[0], colori[0], colori[1], colori[1]]#, 'blue']
		linewidths = [0.5, 1.5, 0.5, 1.5]#, 1.5]
		alphas     = [1.0, 0.7, 1.0, 0.7]#, 0.55]
		labels     = [r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise', \
					  r'$\frac{\partial F(n)}{\partial n}$',   r'$\frac{\partial D(t)}{\partial t}$ - Piecewise']
	for col, ls, lw, l, a, c, m in zip(y, linestyles, linewidths, labels, alphas, colors, markers):
		pd[col].plot(title=titre, ax=ax, ls=ls, lw=lw, label=l, alpha=a, color=c, marker=m)
 
	# ajout des dates spéciales
	if Dates!=None:
		# for d in Dates.listFirstCaseDates:
		#     drawAnnotation(ax, 'First case date\n', d, color='blue')
		for d in Dates.listConfDates:
			drawAnnotation(ax, 'Conf. date\n', d, color='red')
		for d in Dates.listDeconfDates:
			drawAnnotation(ax, 'Deconf. date\n', d, color='green')
		for d in Dates.listOtherDates:
			drawAnnotation(ax, 'Other date\n', d)

	# surlignage des jours de WE
	WE_indices = get_WE_indice(pd)
	i = 0
	while i < len(WE_indices)-1:
		ax.axvspan(pd.index[WE_indices[i]], pd.index[WE_indices[i+1]], facecolor='gray', edgecolor='none', alpha=.15, zorder=-100)
		i += 2

	# ajout d'un text d'annotation
	if textannotation != '':
		bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
		ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/3.), fontsize=6, bbox=bbox, ha="left", va="center") 

	# axes
	ax.grid(True, which='major', axis='both')
	ax.grid(True, which='minor', axis='both')
	ax.grid(True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
	ax.grid(True, which='minor', c='w', lw=0.5, ls='-')
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)

	plt.ticklabel_format(style='plain', axis='y')
	
	# On enlève le label sur l'axe x
	x_label = ax.axes.get_xaxis().get_label().set_visible(False)

	# legende
	legend = ax.legend().get_frame().set_alpha(0.5)
	plt.legend(fontsize=8)
	plt.tight_layout(rect=(0, 0, 1., 0.95))
	plt.savefig(filenameFig, dpi=dpi)
	plt.close()


def PlotFitPiecewiseResidual(modeleString, pd, titre, filenameFig, y, Dates=None, textannotation=''):

	if len(y)==0 or y is None: pass

	# on enlev la premiere et la dernière valeurs
	pd = pd.iloc[1:-1]
	
	fig = plt.figure(facecolor='w', figsize=figsize)
	ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

	# Dessin des courbes
	if modeleString=='SEIR1R2':
		solveur    = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
		colori     = solveur.modele.getColor(3)
		linestyles = ['-']
		markers    = ['x']
		colors     = [colori]
		linewidths = [0.5]
		alphas     = [1.0]
		labels     = [r'$\frac{\partial R^1(n)}{\partial n}-\frac{\partial R^1(t)}{\partial t}$ - Piecewise']
	else:
		solveur    = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)
		colori     = [solveur.modele.getColor(3), solveur.modele.getColor(5)]
		linestyles = ['-', '-']
		markers    = ['x', '+']
		colors     = [colori[0], colori[1]]
		linewidths = [0.5, 0.5]
		alphas     = [1.0, 1.0]
		labels     = [r'$\frac{\partial R^1(n)}{\partial n}-\frac{\partial R^1(t)}{\partial t}$ - Piecewise', \
					  r'$\frac{\partial F(n)}{\partial n}-\frac{\partial D(t)}{\partial t}$ - Piecewise']
	for col, ls, lw, l, a, c, m in zip(y, linestyles, linewidths, labels, alphas, colors, markers):
		pd[col].plot(title=titre, ax=ax, ls=ls, lw=lw, label=l, alpha=a, color=c, marker=m)
 
	# ajout des dates spéciales
	if Dates!=None:
		# for d in Dates.listFirstCaseDates:
		#     drawAnnotation(ax, 'First case date\n', d, color='blue')
		for d in Dates.listConfDates:
			drawAnnotation(ax, 'Conf. date\n', d, color='red')
		for d in Dates.listDeconfDates:
			drawAnnotation(ax, 'Deconf. date\n', d, color='green')
		for d in Dates.listOtherDates:
			drawAnnotation(ax, 'Other date\n', d)

	# surlignage des jours de WE
	WE_indices = get_WE_indice(pd)
	i = 0
	while i < len(WE_indices)-1:
		ax.axvspan(pd.index[WE_indices[i]], pd.index[WE_indices[i+1]], facecolor='gray', edgecolor='none', alpha=.15, zorder=-100)
		i += 2

	# ajout d'un text d'annotation
	if textannotation != '':
		bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
		ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/3.), fontsize=6, bbox=bbox, ha="left", va="center") 

	# axes
	ax.grid(True, which='major', axis='both')
	ax.grid(True, which='minor', axis='both')
	ax.grid(True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
	ax.grid(True, which='minor', c='w', lw=0.5, ls='-')
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)

	plt.ticklabel_format(style='plain', axis='y')

	# On enlève le label sur l'axe x
	x_label = ax.axes.get_xaxis().get_label().set_visible(False)

	# legende
	legend = ax.legend().get_frame().set_alpha(0.5)
	plt.legend(fontsize=8)
	plt.tight_layout(rect=(0, 0, 1., 0.95))
	plt.savefig(filenameFig, dpi=dpi)
	plt.close()

# def PlotAll(data_deriv_allinone, model_deriv_allinone, model_deriv_piecewise, title, filename):

#     fig = plt.figure(facecolor='w', figsize=figsize)
#     ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#     # juste pour récupérer la bonne coukeur
#     solveur = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
#     mod   = solveur.modele
#     color = mod.getColor(3)

#     time = np.linspace(1, len(data_deriv_allinone[1:, 0]), len(data_deriv_allinone[1:, 0]))
#     ax.plot(time, data_deriv_allinone     [1:, 0], label=r'$\frac{\partial R^1(n)}{\partial n}$',              color=color,  alpha=1.0,  lw=0.5, marker='x')
#     ax.plot(time, model_deriv_allinone [1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - All-in-one', color='blue', alpha=0.55, lw=1.5)
#     ax.plot(time, model_deriv_piecewise[1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise',  color=color,  alpha=0.7,  lw=1.5)

#     ax.set_xlabel('Time (days)')
#     ax.yaxis.set_tick_params(length=0)
#     ax.xaxis.set_tick_params(length=0)
#     ax.grid(b=True, which='major', c='w', lw=1, ls='-')
	
#     legend = ax.legend()
#     legend.get_frame().set_alpha(0.5)
#     for spine in ('top', 'right', 'bottom', 'left'):
#         ax.spines[spine].set_visible(False)

#     plt.title(title)
#     plt.savefig(filename, dpi=dpi)
#     plt.close()


# def PlotFit(modeleString, data_deriv, model_deriv, title, filename, ch, textannotation=''):

#     fig = plt.figure(facecolor='w', figsize=figsize)
#     ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#     time = np.linspace(1, len(data_deriv[1:]), len(data_deriv[1:]))

#     if modeleString=='SEIR1R2':
#         solveur = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
#         colori  = [solveur.modele.getColor(3)]
#     elif modeleString=='SEIR1R2D':
#         solveur = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)
#         colori  = [solveur.modele.getColor(3), solveur.modele.getColor(5)]
#     else:
#         print('PlotFit Impossible')
#         exit(1)

#     # les plots
#     ax.plot(time, data_deriv [1:, 0], label=r'$\frac{\partial R^1(n)}{\partial n}$ - ' + ch, color=colori[0], alpha=1.0, lw=0.5, marker='x')
#     ax.plot(time, model_deriv[1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - ' + ch, color=colori[0], alpha=1.0, lw=1.5)
#     if modeleString=='SEIR1R2D':
#         ax.plot(time, data_deriv [1:, 1], label=r'$\frac{\partial D(n)}{\partial n}$ - '   + ch, color=colori[1], alpha=1.0, lw=0.5, marker='+')
#         ax.plot(time, model_deriv[1:, 1], label=r'$\frac{\partial D(t)}{\partial t}$ - '   + ch, color=colori[1], alpha=1.0, lw=1.5)

#     # ajout d'un text d'annotation
#     if textannotation != '':
#         bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
#         ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=6, bbox=bbox, ha="left", va="center") 

#     ax.set_xlabel('Time (days)')
#     ax.yaxis.set_tick_params(length=0)
#     ax.xaxis.set_tick_params(length=0)
#     ax.grid(b=True, which='major', c='w', lw=1, ls='-')
	
#     legend = ax.legend()
#     legend.get_frame().set_alpha(0.5)
#     for spine in ('top', 'right', 'bottom', 'left'):
#         ax.spines[spine].set_visible(False)

#     plt.tight_layout(rect=(0, 0, 1., 0.95))
#     plt.title(title)
#     plt.savefig(filename, dpi=dpi)
#     plt.close()



if __name__ == '__main__':
	main(sys.argv)
