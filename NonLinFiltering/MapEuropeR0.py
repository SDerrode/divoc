#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy                as np
import pandas               as pd
import plotly.graph_objects as go

from common import getRepertoire, getColorMap


# constantes
# figsize       = (15, 15)
alpha         = 0.90
blackstart    = 2.0
dpi           = 120
indexmaxcolor = 2000

def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python3 MapEuropeR0.py R0Moyen_11_12.csv            # default value : SEIR1R2, 0
		>> python3 MapEuropeR0.py R0Moyen_11_12.csv SEIR1R2D 1

		argv[1] : Name of the file to process
		argv[2] : EDO model (SEIR1R2 or SEIR1R2D)                    Default: SEIR2R2
		argv[3] : UKF filtering of data (0/1).                       Default: 0
		argv[4] : Verbose level (debug: 2, ..., almost mute: 0).     Default: 1
	"""

	# constante
	datesestimation='[2020-05-18\u21922020-06-23]'

	##################################################################@
	# Gestion des arguments
	if len(sysargv)>4:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	UKF_filt, UKF_filt01 = False, 0
	modeleString         = 'SEIR1R2'
	verbose              = 1

	if len(sysargv)>1: filename = sysargv[1]
	if len(sysargv)>2: modeleString = sysargv[2]
	if len(sysargv)>3 and int(sysargv[3])==1: UKF_filt, UKF_filt01 = True, 1
	if len(sysargv)>4: verbose = int(sysargv[4])
	filenamewithoutext = os.path.splitext(filename)[0]

	# Lecture des données
	repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/', './figures/'+modeleString+'/')
	df1 = pd.read_csv(repertoire+filename, dtype={'Place': 'str'}) #, parse_dates=[[2]]
	# on filtre les départements dont le R0 n'est pas significatif
	df = df1.query(expr='R0MoyenP2!=-1')
	# on filtre les valeurs supérieures à 2
	# a = np.array(df['R0MoyenP2'].values.tolist())
	# df['R0MoyenP2'] = np.where(a > blackstart, blackstart+0.005, a).tolist()

	minRO = df['R0MoyenP2'].min()
	maxRO = df['R0MoyenP2'].max()
	if verbose>0:
		print('minRO=', minRO)
		print('maxRO=', maxRO)

	mycolormap, newcmp = getColorMap(indexmaxcolor, minRO, maxRO, blackstart, alpha)
	pl_colorscale = matplotlib_to_plotly(mycolormap, indexmaxcolor)

	# Dessin de la carte
	fig = go.Figure(
		data=go.Choropleth(
			locations     = df['Place'], # Spatial coordinates
			z             = df['R0MoyenP2'].astype(float), # Data to be color-coded
			locationmode  = "country names",  # set of locations match entries in `locations`
			colorscale    = 'agsunset', #pl_colorscale,
			#colorbar_title = "R0 value",
		)
	)

	fig.update_layout(
		title_text = 'R0 in Europe ' + datesestimation + ' - model: '+modeleString,
		geo_scope  = 'europe',
	)

	fig.show()
	
# carte de couleurs
def matplotlib_to_plotly(cmap, pl_entries):
	h = 1.0/(pl_entries-1)
	pl_colorscale = []

	for k in range(pl_entries):
		# print('k=', k, ', k*h=', k*h)
		# print('cmap[k,:3]=', cmap[k,:3])
		C = list(map(np.uint8, np.array(cmap[k,:3])*255))
		# print('C=', C)
		pl_colorscale.append([k*h, 'rgba'+str((C[0], C[1], C[2], alpha))])
		#print('pl_colorscale=', pl_colorscale)

	return pl_colorscale

if __name__ == '__main__':
	main(sys.argv)
