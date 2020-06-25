#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Largement inspiré de 
# https://www.data.gouv.fr/fr/reuses/open-source-script-python-pour-visualiser-levolution-des-donnees-covid-19-par-departement-courbes-levolution-sur-la-france-jour-apres-jour-carte/
# script github
# https://github.com/thomasdubdub/covid-france/blob/master/demo-covid.ipynb

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import os
import requests
import zipfile
import io
import pathlib
import numpy                 as np
import matplotlib            as mpl
import matplotlib.pyplot     as plt
import pandas                as pd
import geopandas             as gpd
from datetime                import timedelta, date

from common                  import getRepertoire, getColorMap, save_mapFrance

# constantes
#url_region    = 'http://osm13.openstreetmap.fr/~cquest/openfla/export/regions-20190101-shp.zip'
local_path    = 'tmp/'
name_shp_reg  = 'regions-20190101'
figsize       = (15, 15)
tile_zoom     = 6
indexmaxcolor = 2000
alpha         = 0.70
blackstartP   = 5

def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python MapRegionR0.py R0Moyen_11_12.csv

		argv[1] : Name of the file to process
		argv[2] : EDO model (SEIR1R2 or SEIR1R2D)                Default: SEIR2R2
		argv[3] : UKF filtering of data (0/1).                   Default: 0
		argv[4] : Verbose level (debug: 2, ..., almost mute: 0). Default: 1
	"""

	# constante
	filter_reg      = ['01', '02', '03', '04', '06'] # Guadeloupe, Martinique, Guyane, La Réunion, Mayotte
	datesestimation = '[2020-05-18\u21922020-06-20]'

	##################################################################@
	# Gestion des arguments
	if len(sysargv)>5:
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

	
	##################################################################@	
	# Preparation de la carte de france des régions
	
	# Load French regions data into a GeoPandas GeoSeries
	# lecture à distance
	# r = requests.get(url_region)
	# z = zipfile.ZipFile(io.BytesIO(r.content))
	# z.extractall(path=local_path)
	
	# Lecture en local
	p = pathlib.Path(local_path)
	filen = [j.name for j in p.glob(name_shp_reg+'.*')]
	filenames = [
		y
		for y in sorted(filen)
		for ending in ['dbf', 'prj', 'shp', 'shx']
		if y.endswith(ending)
	]
	
	dbf, prj, shp, shx = [fname for fname in filenames]
	fr = gpd.read_file(local_path + shp)  #  + encoding='utf-8' if needed
	# geométrie
	fr.crs = 'epsg:4326'  # {'init': 'epsg:4326'}
	met = fr.query('code_insee not in @filter_reg')
	met.set_index('code_insee', inplace=True)
	met = met['geometry']

	# Load labelRO data into a pandas DataFrame
	repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/', './figures/'+modeleString+'/')
	df1        = pd.read_csv(repertoire+filename, dtype={'Place': 'str'}) #, parse_dates=[[2]]

	# Le min et le max des 3 colonnes
	a = df1[list(df1)[1:-1]].values
	minRO=np.amin(np.amin(a))
	if minRO==-1.:
		#find th second minimum value (to avoid -1 value whose meaning is to say that the place has a non-meaning RO)
		minRO = np.amin(np.array(a)[a != np.amin(a)])
	maxRO = df1[list(df1)[1:-1]].max().max()
	if verbose>0:
		print('minRO=', minRO)
		print('maxRO=', maxRO)

	# replace les -1 par des nan poiur etre traité comme des données manquantes
	df1.replace(-1, np.nan, inplace=True)
	print(df1.head(15))

	# on rajoute la géométrie
	df1.loc[:, ('geometry')] = df1.loc[:, ('Place')].map(met)
	if verbose>1:
		print(df1.head())
		#input('attente')

	# on rajoute les centroids
	df1['coords'] = df1['geometry'].apply(lambda x: x.centroid.coords[:])
	df1['coords'] = [coords[0] for coords in df1['coords']]
	
	# carte de couleurs (commune aux 3 périodes)
	mycolormap, newcmp = getColorMap(indexmaxcolor, minRO, maxRO, blackstartP, alpha)

	# On dessone les cartes pour les 3 R0
	for p in range(3):
		labelRO = 'R0MoyenP'+str(p)
		print('PROCESSING of', labelRO)
		
		# display the map with the RO data
		img_name = repertoire + filenamewithoutext + '_P' + str(p) + '.png'
		title    = 'R0 par Dpt (vert' + '\u2192' + 'R0 non significatif) \n' + datesestimation + ' - model: '+modeleString
		save_mapFrance(df1, met, newcmp, title, img_name, labelRO, minRO, maxRO, tile_zoom, alpha, figsize)


if __name__ == '__main__':
	main(sys.argv)
