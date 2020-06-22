#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from sklearn.metrics import mean_squared_error

from datetime           import datetime, timedelta
from Covid_SpecialDates import Covid_SpecialDates

dpi           = 120    # plot resolution of saved figures
figsize       = (8, 4) # figure's size (width, height)
strDate       = "%Y-%m-%d"

def getRepertoire(bool, rep1, rep2=''):
	if bool == True:
		repertoire = rep1
	else:
		repertoire = rep2
	if not os.path.exists(repertoire):
		os.makedirs(repertoire)
	return repertoire

def getMaxEQM(sol_edo_R1, data, T):
	dataLength = len(data)
	eqm=np.zeros(shape=(T-dataLength))
	for t in range(T-dataLength):
		eqm[t] = mean_squared_error(sol_edo_R1[t:t+dataLength], data)
	ts0 = np.argmin(eqm)
	return ts0

def midDateStr(startDate, endDate):
	d1 = datetime.strptime(startDate,"%Y-%m-%d")
	d2 = datetime.strptime(endDate,"%Y-%m-%d")
	d  = d1.date() + (d2.date()-d1.date()) / 2
	return d.strftime("%Y-%m-%d")

def addDaystoStrDate(startDate, d):
	d = int(d)
	d1 = datetime.strptime(startDate,"%Y-%m-%d")
	d2 = d1.date() + timedelta(d)
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

def readDates(place, verbose=0):
	'''
		Lecture des dates de confinement et deconfinement pour tous les pays enregistrés
	''' 
	dates_orig = None
	name = './data/dates.csv'
	try:
		dates_orig=pd.read_csv(name, sep=',')#, parse_dates=[1,2,3])
	except:
		print('PB readDates --> exit!')
		exit(1)

	if verbose>1:
		print('TAIL=', dates_orig.tail())

	dates_country = dates_orig.loc[dates_orig['country'] == place]

	Dates = Covid_SpecialDates(country=place)
	Dates.addFirstCaseDates(dates_country.iloc[0]['firstcasedate'])
	Dates.addConfDates     (dates_country.iloc[0]['confdate'])
	Dates.addDeconfDates   (dates_country.iloc[0]['deconfdate'])
	Dates.addOtherDates    (dates_country.iloc[0]['otherdate'])
	
	if verbose>1:
		print(Dates)
	
	return Dates


def GetPairListDates(readStartDate, readStopDate, DatesString, decalage, nbperiodes, recouvrement):

	ListDates = [readStartDate, readStopDate]
	if nbperiodes!=1:
		confin_decalage = None
		if DatesString.listConfDates != []:
			confin_decalage = datetime.strptime(addDaystoStrDate(DatesString.listConfDates[0], decalage-recouvrement), strDate)
			if confin_decalage>ListDates[-2] and confin_decalage<ListDates[-1]:
				ListDates.insert(-1, confin_decalage)
		
		deconfin_decalage = None
		if DatesString.listDeconfDates != []:
			deconfin_decalage = datetime.strptime(addDaystoStrDate(DatesString.listDeconfDates[0], decalage-recouvrement), strDate)
			if deconfin_decalage>ListDates[-2] and deconfin_decalage<ListDates[-1]:
				ListDates.insert(-1, deconfin_decalage)

	ListePairDates = []
	for i in range(len(ListDates)-1):
		if i>0:
			ListePairDates.append((ListDates[i]+timedelta(recouvrement), ListDates[i+1]))
		else:
			ListePairDates.append((ListDates[i], ListDates[i+1]))

	# Conversion date chaine
	ListePairDatesStr = [(date1.strftime(strDate), date2.strftime(strDate)) for date1, date2 in ListePairDates]
	
	return ListePairDates, ListePairDatesStr

def readDataFrance(place='69', dateMinStr=None, dateMaxStr=None, fileLocalCopy=False, verbose=0):
	'''
		Lecture des données du gouvernement français (data.gouv.fr)
		Les données débutent à la date de confinement (pourquoi?)
	''' 
	covid_orig = None
	if fileLocalCopy==True:
		name = './data/csvFrance_2020-06-21.csv'
		try:
			covid_orig=pd.read_csv(name, sep=';', parse_dates=[2], dayfirst=True)
		except:
			fileLocalCopy = False
	
	if fileLocalCopy==False:
		url="https://static.data.gouv.fr/resources/donnees-hospitalieres-relatives-a-lepidemie-de-covid-19/20200504-190020/donnees-hospitalieres-covid19-2020-05-04-19h00.csv"
		url_stable="https://www.data.gouv.fr/fr/datasets/r/63352e38-d353-4b54-bfd1-f1b3ee1cabd7"
		covid_orig=pd.read_csv(url_stable, sep=';', parse_dates=[2], dayfirst=True)
	
	covid_orig.set_index('jour', inplace=True)
	covid_orig.sort_index(inplace=True)

	if verbose>0:
		print('TAIL=', covid_orig.tail())

	covid_orig.drop(columns=['hosp', 'rea'], inplace=True)
	covid_country0 = covid_orig.query(expr='sexe==0').drop(columns=['sexe'])
	covid_country  = covid_country0.loc[covid_country0['dep'] == place]
	
	covid_country1 = covid_country[['rad', 'dc']] #.cumsum()
	if verbose>1:
		print('TAIL=', covid_country1.tail())

	# extraction entre dateMin et dateMaxStr
	if dateMinStr==None:
		dateMinStr = covid_country1.index[0].strftime("%Y-%m-%d")
	if dateMaxStr==None:
		dateMaxStr = addDaystoStrDate(covid_country1.index[-1].strftime("%Y-%m-%d"), 1)
	if verbose>1:
		print('dateMinStr=', dateMinStr, ', dateMaxStr=', dateMaxStr)
	excerpt_country1 = covid_country1.loc[dateMinStr:dateMaxStr].copy()
	# On rajoute la somme des cas et des morts
	excerpt_country1.loc[:, ('radplusdc')] = excerpt_country1.loc[:, ('rad','dc')].sum(axis=1)
	
	if verbose>0:
		print('HEAD=', excerpt_country1.head())
		print('TAIL=', excerpt_country1.tail())
	
	# On recherche la taille de la population en France estimée en 2020
	# site we dont est extrait le fichier local: https://www.insee.fr/fr/statistiques/1893198
	db_pop_size = pd.read_csv('./data/popsizedpt_2020.csv', sep=';')
	pop_size    = int(db_pop_size.loc[place]["popsize"])
	
	return excerpt_country1, list(excerpt_country1), pop_size, dateMinStr, dateMaxStr


def readDataEurope(country='France', dateMinStr=None, dateMaxStr=None, fileLocalCopy=False, verbose=0):
	'''
		Lecture des données recueillies au niveau du site européen
		Remarque: il semble qu'il y ai un décalage d'un jour avec les données françaises
	'''

	covid_orig = None
	if fileLocalCopy==True:
		name = './data/csvEurope_2020-06-11.csv'
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

	if verbose>0:
		print('TAIL=', covid_orig.tail())

	covid_orig.drop(columns=['day', 'year', 'month', 'geoId', 'countryterritoryCode'], inplace=True)
	if verbose>1:
		print('covid_orig.head()=', covid_orig.head())
		print('country=', country)

	covid_country = covid_orig.loc[covid_orig['countriesAndTerritories'] == country]
	# récupération de la taille de la population
	pop_size = int(covid_country[['popData2018']].iloc[1])
	# print('pop_size=', pop_size)
	# input('pause pop')
	covid_country1 = covid_country[['cases', 'deaths']].cumsum()
	# print(covid_country1.head())
	
	# extraction entre dateMin et dateMaxStr
	if dateMinStr==None:
		dateMinStr = covid_country1.index[0].strftime("%Y-%m-%d")
	if dateMaxStr==None:
		dateMaxStr = addDaystoStrDate(covid_country1.index[-1].strftime("%Y-%m-%d"), 1)
	if verbose>1:
		print('dateMinStr=', dateMinStr, ', dateMaxStr=', dateMaxStr)
	excerpt_country1 = covid_country1.loc[dateMinStr:dateMaxStr].copy()
	# On rajoute la somme des cas et des morts
	excerpt_country1.loc[:, ('casesplusdeaths')] = excerpt_country1.loc[:, ('cases','deaths')].sum(axis=1)

	# On rempli les données manquantes avec la données la plus proche
	idx = pd.date_range(start=excerpt_country1.index.min(), end=excerpt_country1.index.max())
	excerpt_country1 = excerpt_country1.reindex(idx, method='nearest')
	# On complète éventuellement au début avec des 0
	idx = pd.date_range(start=datetime.strptime(dateMinStr,"%Y-%m-%d"), end=datetime.strptime(dateMaxStr,"%Y-%m-%d")-timedelta(1))
	excerpt_country1 = excerpt_country1.reindex(idx, fill_value=0.)

	if verbose>0:
		print('TAIL=', excerpt_country1.tail())
	
	return excerpt_country1, list(excerpt_country1), pop_size, dateMinStr, dateMaxStr


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


def PlotData(pd, titre, filenameFig, y, color='black', Dates=None):

	if len(y)==0 or y is None: pass
	
	fig = plt.figure(facecolor='w',figsize=figsize)
	ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

	# Dessin des courbes théoriques
	pd.plot(ax=ax, y=y, color=color, title=titre, marker='x', ls='-', lw=0.5)
	
	# ajout des dates spéciales
	if Dates!=None:
		# for d in Dates.listFirstCaseDates:
		# 	drawAnnotation(ax, 'First case date\n', d, color='blue')
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
	ax.grid(True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
	ax.grid(True, which='minor', c='w', lw=0.5, ls='-')
	for spine in ('top', 'right', 'bottom', 'left'):
		ax.spines[spine].set_visible(False)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False, useLocale=False)

	# On enlève le label sur l'axe x
	x_label = ax.axes.get_xaxis().get_label().set_visible(False)

	# legende
	legend = ax.legend().get_frame().set_alpha(0.8)
	plt.legend(fontsize=7)
	plt.tight_layout()
	plt.savefig(filenameFig, dpi=dpi)
	plt.close()
