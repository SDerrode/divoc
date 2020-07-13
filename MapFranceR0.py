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
import numpy              as np
import matplotlib         as mpl
import matplotlib.pyplot  as plt
import pandas             as pd
import geopandas          as gpd

from datetime             import date, datetime, timedelta
from common               import getRepertoire, getColorMap
from France               import save_mapFranceR0, save_mapFranceI0

strDate = "%Y-%m-%d"

def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python MapFranceR0.py DPT R0Moyen_18_19_DPT.csv
		>> python MapFranceR0.py REG R0Moyen_18_19_REG.csv

		argv[1] : France departments ('DPT') or France regions ('REG')
		argv[2] : Name of the file to process
		argv[3] : EDO model (SEIR1R2 or SEIR1R2D)                Default: SEIR2R2
		argv[4] : UKF filtering of data (0/1).                   Default: 0
		argv[5] : Verbose level (debug: 2, ..., almost mute: 0). Default: 1
	"""

	# constantes
	local_path    = 'shapefileFrance/'
	figsize       = (15, 15)
	tile_zoom     = 6
	indexmaxcolor = 4000
	alpha         = 0.70
	blackstartP   = 100

	##################################################################@
	# Gestion des arguments
	if len(sysargv)>6:
		print('  CAUTION : bad number of arguments - see help')
		exit(1)

	# Default value for parameters
	UKF_filt, UKF_filt01 = False, 0
	modeleString         = 'SEIR1R2'
	verbose              = 1

	if len(sysargv)>1: mapType  = sysargv[1]
	if len(sysargv)>2: filename = sysargv[2]
	if len(sysargv)>3: modeleString = sysargv[3]
	if len(sysargv)>4 and int(sysargv[4])==1: UKF_filt, UKF_filt01 = True, 1
	if len(sysargv)>5: verbose = int(sysargv[5])
	filenamewithoutext = os.path.splitext(filename)[0]
	modeleString2 =f'SEIR\N{SUPERSCRIPT ONE}R\N{SUPERSCRIPT TWO}'

	if mapType=='DPT':
		name_shp    = 'departements-20140306-5m'
		filterArea  = ['971', '972', '973', '974', '976']
		strTile     = '[2020-06-24]'
		textMapType = "Departments"
		substring   = "France-D"
		#textOnMap   = False
		#url_dep    = 'http://osm13.openstreetmap.fr/~cquest/openfla/export/departements-20140306-5m-shp.zip'
	elif mapType=='REG':
		name_shp    = 'regions-20190101'
		filterArea  = ['01', '02', '03', '04', '06'] # Guadeloupe, Martinique, Guyane, La Réunion, Mayotte
		strTile     = '[2020-06-24]'
		textMapType = "Regions"
		substring   = "France-R"
		#textOnMap   = True
		#url_region = 'http://osm13.openstreetmap.fr/~cquest/openfla/export/regions-20190101-shp.zip'
	else:
		print('Only DPT or REG accepted! --> exit!')
		exit(1)

	if verbose>0:
		print('  Full command line : '+sysargv[0]+' '+mapType+' '+filename+' '+modeleString+' '+str(UKF_filt)+' '+str(verbose), flush=True)

	##################################################################@	
	# Preparation de la carte de france (par dpts or par régions)

	# Load French departements data into a GeoPandas GeoSeries
	# lecture à distance
	# r = requests.get(url_dep)
	# z = zipfile.ZipFile(io.BytesIO(r.content))
	# z.extractall(path=local_path)

	# Lecture en local
	p = pathlib.Path(local_path)
	filen = [j.name for j in p.glob(name_shp+'.*')]
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
	met = fr.query('code_insee not in @filterArea')
	met.set_index('code_insee', inplace=True)
	met = met['geometry']

	# Load labelRO data into a pandas DataFrame
	repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/TimeShift/', './figures/'+modeleString+'/TimeShift/')
	df1        = pd.read_csv(repertoire+filename, dtype={'Place': 'str'}) #, parse_dates=[[2]]
	# On simplifie Place pour ne garder que le numéro du dpt ou de la région
	df1['Place'] = df1['Place'].str.replace(substring, '')

	# on rajoute la géométrie
	df1.loc[:, ('geometry')] = df1.loc[:, ('Place')].map(met)
	if verbose>1:
		print(df1.head())
		#input('attente')

	# on rajoute les centroids
	# df1['coords'] = df1['geometry'].apply(lambda x: x.centroid.coords[:])
	# df1['coords'] = [coords[0] for coords in df1['coords']]

	# Le min et le max des 3 colonnes
	a = df1[list(df1)[1:4]].values
	minRO=np.amin(np.amin(a))
	if minRO == -1:
		#find th second minimum value (to avoid -1 value whose meaning is to say that the place has a non-meaning RO)
		minRO = np.amin(np.array(a)[a != np.amin(a)])
		print('-->minRO=', minRO)
	maxRO = df1[list(df1)[1:4]].max().max()
	if verbose>0:
		print('minRO=', minRO)
		print('maxRO=', maxRO)

	# carte de couleurs (commune aux 3 périodes)
	mycolormapR0, newcmpR0 = getColorMap(indexmaxcolor, minRO, maxRO, blackstartP, alpha)
	
	# PARTIE SUR R0
	##################################################################
	# On dessine les cartes pour les 3 R0
	for p in range(3):
		labelRO = 'R0MoyenP'+str(p)
		labelI  = 'IEndP'+str(p)
		print('PROCESSING of', labelRO)

		# # Le min et le max de la colonne
		# a = df1[list(df1)[p+1]].values
		# minRO=np.amin(np.amin(a))
		# if minRO == -1:
		# 	#find th second minimum value (to avoid -1 value whose meaning is to say that the place has a non-meaning RO)
		# 	minRO = np.amin(np.array(a)[a != np.amin(a)])
		# 	print('-->minRO=', minRO)
		# maxRO = df1[list(df1)[p+1]].max().max()
		# if verbose>0:
		# 	print('minRO=', minRO)
		# 	print('maxRO=', maxRO)

		# replace les -1 par des nan pour être traités comme des données manquantes
		df1[list(df1)[p+1]].replace(-1, np.nan, inplace=True)
		
		# display the map with the RO data
		img_name = repertoire + filenamewithoutext + '_P' + str(p) + '.png'
		title    = 'Estimated ' + f'R\N{SUBSCRIPT ZERO}'+ ', scale: ' + textMapType + ' - ' + modeleString2 + ' model'
		save_mapFranceR0(df1, met, newcmpR0, title, img_name, labelRO, labelI, minRO, maxRO, tile_zoom, alpha, figsize, mapType)

	# PARTIE SUR I0
	##################################################################
	
	# on enleve du min, max les dates des territoires exclus durant P0
	for index, row in df1.iterrows():
	     # access data using column names
	     #print(index, row['delay'], row['distance'], row['origin'])
	     if row['R0MoyenP0']==-1.: 
	     	df1.at[index, 'DateFirstCase'] = 'Invalid'
	     	# row['DateFirstCase']='Invalid'
	     	# print('row =', row)
	     	# print('index =', index)
	     	# input('pause')
	dateFirstCase = df1['DateFirstCase'].values.tolist()

	while 'Invalid' in dateFirstCase: 
		dateFirstCase.remove('Invalid') 
	
	dateFirstCase.sort()

	minDateIO = datetime.strptime(dateFirstCase[0],  strDate)
	maxDateIO = datetime.strptime(dateFirstCase[-1], strDate)
	if verbose>0:
		print('minDateIO=', minDateIO)
		print('maxDateIO=', maxDateIO)

	# On dessine la carte des dates du 1er infecté
	img_name = repertoire + filenamewithoutext + '_I0.png'
	title    = 'Date of first infection, scale: ' + textMapType + ' - ' + modeleString2 + ' model'

	dateFirstCase = df1['DateFirstCase'].values.tolist()
	deltaFirstCase = []
	for index in range(len(dateFirstCase)):
		if dateFirstCase[index] != 'Invalid':
			deltaFirstCase.append(float((datetime.strptime(dateFirstCase[index], strDate)-minDateIO).days))
		else:
			deltaFirstCase.append(np.nan)
	df1['deltaI0'] = np.resize(deltaFirstCase,len(df1))
	
	save_mapFranceI0(df1, met, title, img_name, 'deltaI0', 0., float((maxDateIO-minDateIO).days), tile_zoom, alpha, figsize, mapType)


	# # Parse recorded days and save one image for each day
	# vmax = cov1.hosp.max()
	# for i, dt in enumerate(daterange(cov1.index.min(), cov1.index.max())):
	# 	title = dt.strftime('%d-%b-%Y')
	# 	df = cov1.query('jour == @dt')
	# 	df = df.drop_duplicates(subset=['dep'], keep='first')
	# 	img_name = 'figures/' + str(i) + '.png'
	# 	save_img(df, met, title, img_name, 0, vmax)


# def daterange(date1, date2):
# 	for n in range(int((date2 - date1).days) + 1):
# 		yield date1 + timedelta(n)


if __name__ == '__main__':
	main(sys.argv)
