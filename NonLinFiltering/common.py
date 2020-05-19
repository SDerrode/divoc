#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from datetime           import datetime, timedelta
from Covid_SpecialDates import Covid_SpecialDates

dpi           = 120    # plot resolution of saved figures
figsize       = (8, 4) # figure's size (width, height)



def midDateStr(startDate, endDate):
	d1 = datetime.strptime(startDate,"%Y-%m-%d")
	d2 = datetime.strptime(endDate,"%Y-%m-%d")
	d  = d1.date() + (d2.date()-d1.date()) / 2
	return d.strftime("%Y-%m-%d")

def addDaystoStrDate(startDate, d):
	d1 = datetime.strptime(startDate,"%Y-%m-%d")
	d2  = d1.date() + timedelta(d)
	return d2.strftime("%Y-%m-%d")

def getLowerDateFromString(strdate1, strdate2):
	d1 = datetime.strptime(strdate1,"%Y-%m-%d")
	d2 = datetime.strptime(strdate1,"%Y-%m-%d")
	if d1<d2:
		return d1.strftime("%Y-%m-%d")
	else:
		return d2.strftime("%Y-%m-%d")

def getNbDaysBetweenDateFromString(strdate1, strdate2):
	d1 = datetime.strptime(strdate1,"%Y-%m-%d")
	d2 = datetime.strptime(strdate2,"%Y-%m-%d")
	return (d2-d1).days

def getDates(country, verbose):
	Dates = None
	if country == 'France':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates      ('2020-03-16')
		Dates.addDeconfDates    ('2020-05-11')
		Dates.setListOtherDates(['2020-03-06'])
	if country == 'Germany':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates      ('2020-03-22')
		Dates.addDeconfDates    ('2020-04-20')
	if country == 'Italy':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates      ('2020-03-09')
		Dates.addDeconfDates    ('2020-05-04')
	if country == 'Spain':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates      ('2020-03-15')
		Dates.addDeconfDates    ('2020-05-10')
	if country == 'United_Kingdom':
		Dates = Covid_SpecialDates(country=country)
		Dates.addConfDates      ('2020-03-23')
		Dates.addDeconfDates    ('2020-05-13')

	if verbose>1:
		print(Dates)
	
	return Dates

def readDataGouvFr(plot=False):
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


def readDataEurope(country='France', dateMin=None, dateMax=None, plot=False, fileLocalCopy=False, verbose=0):
	'''
		Lecture des données recueillies au niveau du site européen
		Remarque: ll semble qu'il y ai un décalage d'un jour avec les données françaises
	'''

	covid_orig = None
	if fileLocalCopy==True:
		name = './data/csv-2020-05-18'
		try:
			covid_orig=pd.read_csv(name, sep=',', parse_dates=[0], dayfirst=True)
		except:
			fileLocalCopy = False

	if fileLocalCopy==False:
		url="https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
		covid_orig=pd.read_csv(url, sep=',', parse_dates=[0], dayfirst=True)


	#covid_orig.dtypes
	observ_label = ['cases']
	covid_orig.set_index('dateRep', inplace=True)
	covid_orig.sort_index(inplace=True)

	if verbose>0:
		print('TAIL=', covid_orig.tail())

	covid_orig.drop(columns=['day', 'year', 'month', 'geoId', 'countryterritoryCode'], inplace=True)
	if verbose>1:
		print(covid_orig.head())

	covid_country  = covid_orig.loc[covid_orig['countriesAndTerritories'] == country]
	covid_country1 = covid_country[['cases', 'deaths']].cumsum()
	# récupération de la taille de la population
	pop_size = int(covid_country[['popData2018']].iloc[0])

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
	
	return excerpt_country1, observ_label, pop_size


def get_WE_indice(pd):
	WE_indices = []
	BoolWE = False
	for i in range(len(pd)):
		if pd.index[i].weekday() >= 5 and BoolWE == False:
			WE_indices.append(i) 
			BoolWE = True
		if pd.index[i].weekday() < 5 and BoolWE == True:
			WE_indices.append(i) 
			BoolWE = False
	# refermer si ouvert
	if BoolWE==True:
		WE_indices.append(i)
	# print('WE_indices=', WE_indices)
	# input('weekday')
	return WE_indices 

def drawAnnotation(ax, strin, date, color='black'):
	bbox=dict(boxstyle='round4,pad=.3', fc='0.9', ec=color, lw=0.5)
	arrowprops=dict(arrowstyle="->", color=color, lw=0.5)
	ax.annotate(strin+date, xy=(date, ax.get_ylim()[0]), xycoords='data', xytext=(date, ax.get_ylim()[0]-(ax.get_ylim()[1]-ax.get_ylim()[0])/6.), \
			fontsize=6, bbox=bbox, arrowprops=arrowprops, ha="center", va="center")


def Plot(pd, titre, NameFig, modele, y, Dates=None, data=None):

	if len(y)==0 or y is None: pass

	listeString, listeColor = [], []
	for p in y:
		listeString.append(modele.getString(p))
		listeColor.append(modele.getColor(p))

	fig = plt.figure(facecolor='w',figsize=figsize)
	ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

	# Dessin des courbes théoriques
	pd.plot(ax=ax, y=listeString, color=listeColor, title=titre)
	# Dessin des observations
	if data != None:
		pd.plot(ax=ax, y=data, marker='x', ls='', color=modele.getColorFromString('$R^1(t)$') )
	
	# ajout des dates spéciales
	if Dates!=None:
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

	# axes
	ax.grid(True, which='major', axis='both')
	ax.grid(True, which='minor', axis='both')
	ax.grid(b=True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
	ax.grid(b=True, which='minor', c='w', lw=0.5, ls='-')
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False, useLocale=False)

	# On enlève le label sur l'axe x
	x_label = ax.axes.get_xaxis().get_label().set_visible(False)

	# legende
	legend = ax.legend().get_frame().set_alpha(0.8)
	plt.legend(fontsize=7)
	
	plt.tight_layout()
	plt.savefig(NameFig, dpi=dpi)
	plt.close()
