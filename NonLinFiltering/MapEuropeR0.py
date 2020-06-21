#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import plotly.graph_objects as go
import pandas as pd

from common import getRepertoire


def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python3 MapEurope.py             # default value : SEIR1R2, 0
		>> python3 MapEurope.py SEIR1R2D 1
	"""

	# constante
	filename   = 'R0moyen_10_16.csv'
	filename   = 'R0moyen_7_12.csv'

	if len(sysargv) > 3:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	UKF_filt, UKF_filt01 = False, 0
	modeleString = 'SEIR1R2'

	if len(sysargv)>1: modeleString = sysargv[1]
	if len(sysargv)>2 and int(sysargv[2])==1: UKF_filt, UKF_filt01 = True, 1

	repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/', './figures/'+modeleString+'/')
	
	# Lecture des données
	df = pd.read_csv(repertoire+filename)

	minRO = df['R0moyen'].min()
	maxRO = df['R0moyen'].max()

	
	# Dessin de la carte
	fig = go.Figure(data=go.Choropleth(
	    locations=df['Country'], # Spatial coordinates
	    z = df['R0moyen'].astype(float), # Data to be color-coded
	    locationmode = "country names", # set of locations match entries in `locations`
	    colorscale = [[0, 'green'], [(1.-minRO)/(maxRO-minRO), 'white'], [1, 'red']],
	    colorbar_title = "R0 value",
	))

	fig.update_layout(
	    title_text = 'Europe ' + '$R_0$' + '(date to define) - molde: '+modeleString,
	    geo_scope  = 'europe',
	)

	fig.show()
	
if __name__ == '__main__':
	main(sys.argv)
