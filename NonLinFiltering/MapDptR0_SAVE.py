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
import contextily            as ctx
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime                import timedelta, date

from common                  import getRepertoire, getColorMap

# constantes
#url_dep       = 'http://osm13.openstreetmap.fr/~cquest/openfla/export/departements-20140306-5m-shp.zip'
local_path    = 'tmp/'
name_shp_dpt  = 'departements-20140306-5m'
figsize       = (15, 15)
tile_zoom     = 6
indexmaxcolor = 2000
alpha         = 0.70
blackstartP   = 5
dpi           = 120

def main(sysargv):

	"""
		:Example:

		For countries (European database)
		>> python MapDptR0.py R0Moyen_11_12.csv

		argv[1] : Name of the file to process
		argv[2] : EDO model (SEIR1R2 or SEIR1R2D)                    Default: SEIR2R2
		argv[3] : UKF filtering of data (0/1).                       Default: 0
		argv[4] : Verbose level (debug: 2, ..., almost mute: 0).     Default: 1
	"""

	# constante
	filter_dep      = ['971', '972', '973', '974', '976']
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
	# Preparation de la carte de france des dpts

	# Load French departements data into a GeoPandas GeoSeries
	# lecture à distance
	# r = requests.get(url_dep)
	# z = zipfile.ZipFile(io.BytesIO(r.content))
	# z.extractall(path=local_path)

	# Lecture en local
	p = pathlib.Path(local_path)
	filen = [j.name for j in p.glob(name_region+'.*')]
	filenames = [
		y
		for y in sorted(filen)
		for ending in ['dbf', 'prj', 'shp', 'shx']
		if y.endswith(ending)
	]

	dbf, prj, shp, shx = [fname for fname in filenames]
	fr = gpd.read_file(local_path + shp)  #  + encoding='utf-8' if needed
	fr.crs = 'epsg:4326'  # {'init': 'epsg:4326'}
	met = fr.query('code_insee not in @filter_dep')
	met.set_index('code_insee', inplace=True)
	met = met['geometry']

	# Load labelRO data into a pandas DataFrame
	repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/', './figures/'+modeleString+'/')
	df1        = pd.read_csv(repertoire+filename, dtype={'Place': 'str'}) #, parse_dates=[[2]]

	# Le min et le max des 3 colonnes
	minRO = df1[list(df1)[1:-1]].min().min()
	maxRO = df1[list(df1)[1:-1]].max().max()
	if verbose>0:
		print('minRO=', minRO)
		print('maxRO=', maxRO)

	# On dessone les cartes pour les 3 R0
	for p in range(3):
		labelRO = 'R0MoyenP'+str(p)
		print('PROCESSING of', labelRO)

		# on filtre les départements dont le R0 n'est pas significatif
		if p == 0:
			df = df1.query(expr='R0MoyenP0!=-1.')
		elif p == 1:
			df = df1.query(expr='R0MoyenP1!=-1.')
		else:
			df = df1.query(expr='R0MoyenP2!=-1.')

		# on filtre les valeurs supérieures à 2 pour P3 uniquement
		if p==2:
			blackstart = 1000
			a = np.array(df[labelRO].values.tolist())
			df[labelRO] = np.where(a > blackstart, blackstart+0.005, a).tolist()
		else:
			blackstart = 1000

		# Add geometry data to COVID DataFrame
		df['geometry'] = df['Place'].map(met)
		if verbose>1:
			print(df.head())
			#input('attente')

		# carte de couleurs
		mycolormap, newcmp = getColorMap(indexmaxcolor, minRO, maxRO, blackstart, alpha)

		# display the map with the RO data
		img_name = repertoire + filenamewithoutext + '_P' + str(p) + '.png'
		title    = 'R0 par Dpt (sans couleur' + '\u2192' + 'R0 non significatif) ' + datesestimation + ' - model: '+modeleString
		save_img(df, met, newcmp, title, img_name, labelRO, minRO, maxRO)


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


def save_img(df, met, newcmp, title, img_name, labelRO, vmin, vmax):
	
	# Load the map tile with contextily
	w, s, e, n = met.total_bounds
	bck, ext = ctx.bounds2img(w, s, e, n, zoom=tile_zoom, ll=True)

	gdf = gpd.GeoDataFrame(df, crs={'init': 'epsg:4326'})
	gdf_3857 = gdf.to_crs(epsg=3857)  # web mercator
	f, ax = plt.subplots(figsize=figsize)
	ax.imshow(
		bck, extent=ext, interpolation='sinc', aspect='equal'
	)  # load background map
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.15)  # GeoPandas trick to adjust the legend bar
	gdf_3857.plot(
		column=labelRO,  # R0 for the county
		ax        = ax,
		cax       = cax,
		alpha     = alpha,
		edgecolor = 'k',
		legend    = True,
		cmap      = newcmp,
		vmin      = vmin,
		vmax      = vmax,
	)

	ax.set_axis_off()
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	ax.set_title(title, fontsize=20)
	plt.savefig(img_name, bbox_inches='tight') # to remove border
	plt.close(f)

if __name__ == '__main__':
	main(sys.argv)
