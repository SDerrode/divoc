#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

dpi           = 100    # plot resolution of saved figures
figsize       = (8, 4) # figure's size (width, height)


def midDateStr(startDate, endDate):
	d1 = datetime.strptime(datemin,"%Y-%m-%d")
	d2 = datetime.strptime(dateconf,"%Y-%m-%d")
	d  = d1.date() + (d2.date()-d1.date()) / 2

def readDataGouvFr(Plot=False):
	'''
		Lecture des données du gouvernement français (data.gouv.fr)
		Les données débutent à la date de confinement (pourquoi?)
	''' 

	url="https://static.data.gouv.fr/resources/donnees-hospitalieres-relatives-a-lepidemie-de-covid-19/20200504-190020/donnees-hospitalieres-covid19-2020-05-04-19h00.csv"
	url_stable="https://www.data.gouv.fr/fr/datasets/r/63352e38-d353-4b54-bfd1-f1b3ee1cabd7"
	covid_orig=pd.read_csv(url_stable, sep=';', parse_dates=[2])
	covid_orig.set_index('jour', inplace=True)

	covid = covid_orig.query(expr='sexe==0').drop(columns=['dep', 'sexe'])
	covid1 = covid.groupby(covid.index)[['hosp', 'rea', 'rad', 'dc']].sum()
	print(covid1.head())

	if Plot==True:
		ax=covid1.plot(figsize=(16, 4))
		#covid1.plot(subplots=True, figsize=(6, 6))
		style = dict(size=12, color='magenta')
		#ax.text('2020-03-17', 10000, "Confine date", **style)
		ax.annotate('Fr confine date: 2020-03-18', xy=('2020-03-18', 1), xycoords='data', xytext=('2020-03-30', 40000), bbox=dict(boxstyle="round4,pad=.5", fc="0.9"), arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=0,angleB=-90"));
		plt.savefig('toto.png')

	return  


def readDataEurope(country='France', dateMin=None, dateMax=None, Plot=False, fileLocalCopy=False, verbose=0):
	'''
		Lecture des données receuilli au niveau de  (data.gouv.fr)
		Remarque: ll semble qu'il y ai un décalage d'un jour avec les données françaises
	'''

	covid_orig = None
	if fileLocalCopy==True:
		name = './data/csv_2020-05-06'
		try:
			covid_orig=pd.read_csv(name, sep=',', parse_dates=[0], dayfirst=True)
		except:
			fileLocalCopy = False

	if fileLocalCopy==False:
		url="https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
		covid_orig=pd.read_csv(url, sep=',', parse_dates=[0], dayfirst=True)


	#covid_orig.dtypes
	covid_orig.set_index('dateRep', inplace=True)
	covid_orig.sort_index(inplace=True)

	covid_orig.drop(columns=['day', 'year', 'month', 'geoId', 'countryterritoryCode'], inplace=True)
	if verbose>1:
		print(covid_orig.head())

	covid_country = covid_orig.loc[covid_orig['countriesAndTerritories'] == country]
	#covid_country.dtypes
	#covid_country.tail()
	covid_country1   = covid_country[['cases', 'deaths']].cumsum()

	# extraction entre dateMin et dateMax
	if verbose>0:
		print('dateMin=', dateMin, ', dateMax=', dateMax)
	if dateMin==None:
		if dateMax==None:
			excerpt_country1 = covid_country1.loc[:]
		else:
			excerpt_country1 = covid_country1.loc[:dateMax]
	else:
		excerpt_country1 = covid_country1.loc[dateMin:dateMax]
	
	if verbose>0:
		print('TAIL=', excerpt_country1.tail())
	
	observ_label = ['cases']
	return excerpt_country1, observ_label

def drawAnnotation(ax, strin, date, color='black'):
	bbox=dict(boxstyle='round4,pad=.3', fc='0.9', ec=color, lw=0.5)
	arrowprops=dict(arrowstyle="->", color=color, lw=0.5)
	ax.annotate(strin+date, xy=(date, ax.get_ylim()[0]), xycoords='data', xytext=(date, ax.get_ylim()[0]-(ax.get_ylim()[1]-ax.get_ylim()[0])/6.), \
		 	fontsize=6, bbox=bbox, arrowprops=arrowprops, ha="center", va="center")

def Plot(pd, y, country, NameFig, Dates=None, z_observ=None):
		fig, ax = plt.subplots(figsize=figsize)
		pd.plot(ax=ax, y=y, title=country)
		if z_observ != None:
			pd.plot(ax=ax, y=z_observ, marker='x', ls='', color='green')
		ax.grid(True, which='major', axis='both')
		ax.grid(True, which='minor', axis='both')

		if Dates!=None:
			for d in Dates.listConfDates:
				drawAnnotation(ax, 'Conf. date\n', d, color='red')
			for d in Dates.listDeconfDates:
				drawAnnotation(ax, 'Deconf. date\n', d, color='green')
			for d in Dates.listOtherDates:
				drawAnnotation(ax, 'Other date\n', d)

		plt.tight_layout()
		plt.savefig(NameFig, dpi=dpi)
		plt.close()

