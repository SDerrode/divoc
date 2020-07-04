#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy             as np
import matplotlib.pyplot as plt
import contextily        as ctx
import geopandas         as gpd
import seaborn           as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

dpi=300

def save_mapFranceR0(df, met, newcmp, title, img_name, labelRO, vmin, vmax, tile_zoom, alpha, figsize, mapType):
	
	# Load the map tile with contextily
	w, s, e, n = met.total_bounds
	bck, ext   = ctx.bounds2img(w, s, e, n, zoom=tile_zoom, ll=True)

	gdf      = gpd.GeoDataFrame(df, crs={'init': 'epsg:4326'})	
	gdf_3857 = gdf.to_crs(epsg=3857)  # web mercator

	# on rajoute les centroids
	gdf_3857['coords'] = gdf_3857['geometry'].apply(lambda x: x.centroid.coords[:])
	gdf_3857['coords'] = [coords[0] for coords in gdf_3857['coords']]

	# Prepare the figure
	fig, (ax, hax) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [10, 1]}, figsize=(12, 6))

	#######################@ Pour ax
	ax.imshow(bck, extent=ext, interpolation='sinc', aspect='equal', alpha=0.4)  # load background map
	
	# Plot the map with transparency
	divider = make_axes_locatable(ax)
	cax     = divider.append_axes('right', size='5%', pad=0.15)  # GeoPandas trick to adjust the legend bar
	gdf_3857.plot(
		column    = labelRO,  # R0 for the places
		ax        = ax,
		cax       = cax,
		alpha     = alpha,
		edgecolor = 'k',
		linewidth = 0.7,
		legend    = True,
		cmap      = newcmp,
		vmin      = vmin,
		vmax      = vmax,
		#legend_kwds={'label': "Mean R0"},
		#missing_kwds={'color': 'lightgrey'},
		missing_kwds={
				"color": "lightgreen",
				"edgecolor": "k",
				"alpha": 0.5,
				"linewidth": 0.7,
				#"hatch": "///",
				"label": "Missing values",
			},
	)

	if mapType=='REG':
		fontsize = 10
		weight   = 'bold'
		for _, row in gdf_3857.iterrows():
			if np.isnan(row[labelRO])==False:
				s='#'+row['Place']+'\n'+f'R\N{SUBSCRIPT ZERO}:' + str(round(row[labelRO], 2))
			else:
				s='#'+row['Place']
			ax.text(s=s, x=row['coords'][0], y = row['coords'][1],
				   horizontalalignment='center', verticalalignment='center', fontdict = {'size': fontsize, 'weight': weight})
	if mapType=='DPT':
		fontsize = 7
		weight   = 'normal' #bold, normal
		for _, row in gdf_3857.iterrows():
			if row['Place'] not in ['75', '92', '93', '94']:
				s='#'+row['Place']
				ax.text(s=s, x=row['coords'][0], y = row['coords'][1],
					   horizontalalignment='center', verticalalignment='center', fontdict = {'size': fontsize, 'weight': weight})


	# Affine la visu du plot
	ax.set_axis_off()
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	ax.set_title(title, fontsize=14)

	#######################@ Pour hax
	sns.distplot(df[labelRO], kde=True, 
			 vertical=True,
			 bins=int(20), color = 'darkblue',
			 hist_kws={'edgecolor': 'black'},
			 kde_kws={'linewidth': 2}, ax=hax,
			 rug=True, hist=False)
	hax.set_ylim([vmin, vmax])
	hax.set_ylabel('')
	# Affine la visu du plot
	hax.set_axis_off()
	hax.get_xaxis().set_visible(False)
	hax.get_yaxis().set_visible(False)

	plt.subplots_adjust(wspace = -0.45, hspace = -0.45)
	#plt.subplots_adjust(wspace = 0., hspace = -0.)

	plt.savefig(img_name, bbox_inches='tight', dpi=dpi) # to remove border
	plt.close(fig)


def save_mapFranceI0(df, met, title, img_name, deltaIO, vmin, vmax, tile_zoom, alpha, figsize, mapType):

	# Load the map tile with contextily
	w, s, e, n = met.total_bounds
	bck, ext   = ctx.bounds2img(w, s, e, n, zoom=tile_zoom, ll=True)

	gdf      = gpd.GeoDataFrame(df, crs={'init': 'epsg:4326'})	
	gdf_3857 = gdf.to_crs(epsg=3857)  # web mercator

	# on rajoute les centroids
	gdf_3857['coords'] = gdf_3857['geometry'].apply(lambda x: x.centroid.coords[:])
	gdf_3857['coords'] = [coords[0] for coords in gdf_3857['coords']]

	# Prepare the figure and plot the background
	fig, (ax, hax) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [10, 1]}, figsize=(12, 6))

	#######################@ Pour ax
	ax.imshow(bck, extent=ext, interpolation='sinc', aspect='equal', alpha=0.4)  # load background map
	
	# Plot the map with transparency
	divider = make_axes_locatable(ax)
	cax     = divider.append_axes('right', size='5%', pad=0.15)  # GeoPandas trick to adjust the legend bar
	gdf_3857.plot(
		column    = deltaIO,  # delta I0 for the places
		ax        = ax,
		cax       = cax,
		alpha     = alpha,
		edgecolor = 'k',
		linewidth = 0.7,
		legend    = True,
		cmap      = 'YlOrRd_r',
		vmin      = vmin,
		vmax      = vmax,
		#legend_kwds={'label': "Dates I0 from " + minDateIO.strftime("%Y-%m-%d")},
		missing_kwds={
				"color": "lightgreen",
				"edgecolor": "k",
				"alpha": 0.5,
				"linewidth": 0.7,
				#"hatch": "///",
				"label": "Missing values",
			},
	)

	if mapType=='REG':
		fontsize = 10
		weight   = 'bold'
		for _, row in gdf_3857.iterrows():
			if row['DateFirstCase'] != 'Invalid':
				s='#'+row['Place']+'\n'+str(row['DateFirstCase']).replace('2020-', '')
			else:
				s='#'+row['Place']
			ax.text(s=s, x=row['coords'][0], y = row['coords'][1],
				   horizontalalignment='center', verticalalignment='center', fontdict = {'size': fontsize, 'weight': weight})
	if mapType=='DPT':
		fontsize = 7
		weight   = 'normal' #bold, normal
		for _, row in gdf_3857.iterrows():
			if row['Place'] not in ['75', '92', '93', '94']:
				s='#'+row['Place']
				ax.text(s=s, x=row['coords'][0], y = row['coords'][1],
					   horizontalalignment='center', verticalalignment='center', fontdict = {'size': fontsize, 'weight': weight})

	
	# Zoom sur la région parisienne
	# create some data to use for the plot
	# dt = 0.001
	# t = np.arange(0.0, 10.0, dt)
	# r = np.exp(-t[:1000]/0.05)               # impulse response
	# x = np.random.randn(len(t))
	# s = np.convolve(x, r)[:len(x)]*dt  # colored noise

	# a = plt.axes([.13, .68, .17, .17], facecolor='y')
	# #n, bins, patches = plt.hist(s, 400, normed=1)
	# #plt.title('Probability')
	# plt.xticks([])
	# plt.yticks([])

	# Affine la visu du plot
	ax.set_axis_off()
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	ax.set_title(title, fontsize=14)

	#######################@ Pour hax
	sns.distplot(df[deltaIO], kde=True, 
			 vertical=True,
			 bins=int(20), color = 'darkblue',
			 hist_kws={'edgecolor': 'black'},
			 kde_kws={'linewidth': 2}, ax=hax,
			 rug=True, hist=False)
	hax.set_ylim([vmin, vmax])
	hax.set_ylabel('')
	# Affine la visu du plot
	hax.set_axis_off()
	hax.get_xaxis().set_visible(False)
	hax.get_yaxis().set_visible(False)

	# Enregistrement de la figure
	plt.subplots_adjust(wspace = -0.45, hspace = -0.45)
	plt.savefig(img_name, bbox_inches='tight', dpi=dpi) # to remove border
	plt.close(fig)


def getRegionFromCodeInsee(codeR):
	
	if '+' in codeR:
		code=int(codeR[1:-1])
	else:
		code=int(codeR[1:])

	if code==84:
		listDpt=[['D01'],['D03'],['D07'],['D15'],['D26'],['D38'],['D42'],['D43'],['D63'],['D69'],['D73'],['D74']]
		codeISO='FR-ARA'
	if code==27:
		listDpt=[['D21'],['D25'],['D39'],['D58'],['D70'],['D71'],['D89'],['D90']]
		codeISO='FR-BFC'
	if code==53:
		listDpt=[['D22'],['D29'],['D35'],['D56']]
		codeISO='FR-BRE'
	if code==24:
		listDpt=[['D18'],['D28'],['D36'],['D37'],['D41'],['D45']]
		codeISO='FR-CVL'
	if code==94:
		listDpt=[['D2A'],['D2B']]
		codeISO='FR-COR'
	if code==44:
		listDpt=[['D08'],['D10'],['D51'],['D52'],['D54'],['D55'],['D57'],['D67'],['D68'],['D88']]
		codeISO='FR-GES'
	if code==32:
		listDpt=[['D02'],['D59'],['D60'],['D62'],['D80']]
		codeISO='FR-HDF'
	if code==11:
		listDpt=[['D75'],['D77'],['D78'],['D91'],['D92'],['D93'],['D94'],['D95']]
		codeISO='FR-IDF'
	if code==28:
		listDpt=[['D14'],['D27'],['D50'],['D61'],['D76']]
		codeISO='FR-NOR'
	if code==75:
		listDpt=[['D16'],['D17'],['D19'],['D23'],['D24'],['D33'],['D40'],['D47'],['D64'],['D79'],['D86'],['D87']]
		codeISO='FR-NAQ'
	if code==76:
		listDpt=[['D09'],['D11'],['D12'],['D30'],['D31'],['D32'],['D34'],['D46'],['D48'],['D65'],['D66'],['D81'],['D82']]
		codeISO='FR-OCC'
	if code==52:
		listDpt=[['D44'],['D49'],['D53'],['D72'],['D85']]
		codeISO='FR-PDL'
	if code==93:
		listDpt=[['D04'],['D05'],['D06'],['D13'],['D83'],['D84']]
		codeISO='FR-PAC'

	if '+' in codeR:
		return [listDpt]
	else:
		return listDpt


def getPlace(element):
	
	if element[0]=='D': 
		dpt = [[element]]
		return dpt, dpt
	if element[0]=='R': 
		region = getRegionFromCodeInsee(element)
		if element[-1]=='+':
			return region, [[element]]
		else:
			return region, region
	if element=='MetropoleR+': 
		return [getRegionFromCodeInsee('R84+'), getRegionFromCodeInsee('R27+'), \
				getRegionFromCodeInsee('R53+'), getRegionFromCodeInsee('R24+'), \
				getRegionFromCodeInsee('R94+'), getRegionFromCodeInsee('R44+'), \
				getRegionFromCodeInsee('R32+'), getRegionFromCodeInsee('R11+'), \
				getRegionFromCodeInsee('R28+'), getRegionFromCodeInsee('R75+'), \
				getRegionFromCodeInsee('R76+'), getRegionFromCodeInsee('R52+'), \
				getRegionFromCodeInsee('R93+'), \
			], [['R84+'], ['R27+'], ['R53+'], ['R24+'], ['R94+'], ['R44+'], ['R32+'], ['R11+'], ['R28+'], ['R75+'], ['R76+'], ['R52+'], ['R93+']]
	if element=='MetropoleD' or element=='MetropoleD+': 
		listdpts = []
		for i in range(1,96):
			number_str = str(i)
			zero_filled_number = number_str.zfill(2)
			listdpts.append(['D'+zero_filled_number])
		listdpts.remove(['D20']) # Ce département n'est pas dans les données
		listdpts.append(['D2A'])
		listdpts.append(['D2B'])
		if element=='MetropoleD':
			return listdpts, listdpts
		else:
			return [listdpts], [['MetropoleD+']]

if __name__ == '__main__':
	# print(getRegionFromCodeInsee('R52'))
	# print(getRegionFromCodeInsee('R52+'))


	Chaine='FRANCE,R84+'
	ChaineSplit = Chaine.split(',')
	ListePlace = []
	if ChaineSplit[0]=='FRANCE':
		for el in ChaineSplit[1:]:
			ListePlace.extend(getPlace(el))

	# Chaine='FRANCE,MetropoleR+,R93+,D19'
	# ChaineSplit = Chaine.split(',')
	# ListePlace = []
	# if ChaineSplit[0]=='FRANCE':
	# 	for el in ChaineSplit[1:]:
	# 		ListePlace.extend(getPlace(el))

	# Chaine='FRANCE,MetropoleD'
	# ChaineSplit = Chaine.split(',')
	# ListePlace = []
	# if ChaineSplit[0]=='FRANCE':
	# 	for el in ChaineSplit[1:]:
	# 		ListePlace.extend(getPlace(el))

	print('Chaine=', Chaine)
	print('ListePlace', ListePlace)

	for el in ListePlace:
		print('el=', el)
