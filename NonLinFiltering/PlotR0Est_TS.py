#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd

#from matplotlib.ticker import MaxNLocator
from datetime          import datetime, timedelta
  
from common            import addDaystoStrDate, getRepertoire, get_WE_indice
# from France            import getPlace
# from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2
# from SolveEDO_SEIR1R2D import SolveEDO_SEIR1R2D
# from ProcessSEIR1R2D   import fit as fitProcessSEIR1R2D
# from ProcessSEIR1R2    import fit as 
from PlotTimeShift       import PlotTimeShift


dpi     = 300    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)
strDate = "%Y-%m-%d"

def main(sysargv):

	print('sysargv=', sysargv)

	startAnalyseDate = "2020-06-01"
	TS               = 13
	lenAnalysis      = 20
	sexe, sexestr    = 0, 'male+female'
	modelString      = 'SEIR1R2'
	UKF_filt         = 0
	verbose          = 1
	plot             = 1

	#ListPlaces = 'FRANCE,MetropoleD+,MetropoleR+'
	ListPlaces = 'FRANCE,MetropoleD+'
	#ListPlaces = 'FRANCE,MetropoleD+,R84+,R52+'
	#ListPlaces = 'FRANCE,D28,D84'

	repertoire = getRepertoire(UKF_filt, './figures/'+modelString+'_UKFilt/R0Est_TS_sexe_'+str(sexe)+'_shift_'+str(TS), './figures/'+modelString+'/R0Est_TS_sexe_'+str(sexe)+'_shift_'+str(TS))

	if len(sysargv) == 1:

		# Recupération des données et creation d'un dataframe pandas
		#################################################################

		df = pd.DataFrame(columns = ['EndEstDate','Place','R0MoyenP2'])
		for delta in range(lenAnalysis): #range(0, 13, 5):
			stopDate = addDaystoStrDate(startAnalyseDate, delta)
			print('-->STOP DATE=', stopDate)

			# Preparation de la liste des arguments por appel aux fonctions
			ListeArg=['PlotTimeShift.py', ListPlaces, str(sexe), \
					modelString, str(UKF_filt), str(TS)+','+str(TS+1), str(0), str(plot), stopDate]
			
			# Appel à la fonction
			ListeChaines = PlotTimeShift(ListeArg)
			if verbose>1:
				print('ListeChaines=', ListeChaines)
			
			# recupération en enregistrement des infos pour tous les sites demandés
			for elt in ListeChaines:
				ListeElt = elt.split(',')
				new_row = {'EndEstDate':stopDate, 'Place':ListeElt[0], 'R0MoyenP2':np.float64(ListeElt[3])}
				df = df.append(new_row, ignore_index=True)

		# Sauvegarde des résultats dans un fichier
		filename = repertoire+'/ROEst_TS.csv'
		print('filename=', filename)
		df.to_csv(filename, index=False)

	else:
		df = pd.read_csv(sysargv[1], sep=',', dtype={'R0MoyenP2':np.float64})

	# listeHeader = list(df)
	df['EndEstDate']    = pd.to_datetime(df['EndEstDate'])
	df.set_index('EndEstDate', inplace=True)
	# print('df.head()=', df.head())
	# input('attente')

	# Filtage de quelques régions
	# filterArea=['53', '76', '94']
	# print(df.info())
	# df.query('Place not in @filterArea', inplace=True)
	# print(df.head(15))

	# Affichages graphiques
	#################################################################
	# listplaces = ListPlaces.split(',')[1:]

	# Répertoire des figures
	prefFig = repertoire + '/TS_'

	# On s'occupe de l'évolution du R0 sur la période 2
	filename = prefFig + str(TS) + '_R0P2.png'
	if sexe==0:	
		title    = f'R\N{SUBSCRIPT ZERO}' + ' estimation w.r.t. to the end date for estimation - Delay (delta)=' + str(TS) + ' day(s)'
	else:
		title    = f'R\N{SUBSCRIPT ZERO}' + ' estimation w.r.t. to the end date for estimation - Sex=' + sexestr + ', Delay (delta)=' + str(TS) + ' day(s)'

	y='R0MoyenP2'
	PlotR0Est_TS(modelString, df, y, title, filename)

	df1 = df.groupby('Place').get_group('Metropole')
	filename = prefFig + str(TS) + '_R0P2Metropole.png'
	PlotR0Est_TS(modelString, df1, y, title, filename)

	print('Mean of estimated R0=', df.groupby('Place').mean())
	print('Std  of estimated R0=', df.groupby('Place').std())


def PlotR0Est_TS(modelString, df, y, title, filename):

	fig = plt.figure(facecolor='w', figsize=figsize)
	ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

	import seaborn as sns
	sns.set_palette(sns.color_palette("hls", 20))

	# Dessin
	df.groupby('Place')[y].plot(title=title, legend=True)

	# surlignage des jours de WE
	WE_indices = get_WE_indice(df)
	i = 0
	while i < len(WE_indices)-1:
		ax.axvspan(df.index[WE_indices[i]], df.index[WE_indices[i+1]], facecolor='gray', edgecolor='none', alpha=.15, zorder=-100)
		i += 2

	# axes
	ax.grid(True, which='major', axis='both')
	ax.grid(True, which='minor', axis='both')
	ax.grid(True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
	ax.grid(True, which='minor', c='w', lw=0.5, ls='-')
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False, useLocale=False)

	# On enlève le label sur l'axe x
	x_label = ax.axes.get_xaxis().get_label().set_visible(False)
	# On enleve les exponeent sur l'axe Y
	#ax.ticklabel_format(useOffset=False)

	# legende
	legend = ax.legend().get_frame().set_alpha(0.5)
	ax.legend(fontsize=8)
	ax.set_ylim([0.5, 1.1])
	fig.tight_layout(rect=(0, 0, 1., 0.95))
	fig.savefig(filename, dpi=dpi)
	plt.close()


if __name__ == '__main__':
	main(sys.argv)


