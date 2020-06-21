#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy         as np
from datetime        import datetime as dt
from dateutil.parser import parse


class Covid_SpecialDates:

	def __init__(self, country):

		self.country = country

		# Empty lists of dates
		self.listFirstCaseDates = []
		self.listConfDates      = []
		self.listDeconfDates    = []
		self.listOtherDates     = []

	def __str__(self):
		S =  'Interesting dates for ' + self.country + ':' 
		S += '\n  First case date(s)   :' + str(self.listFirstCaseDates)
		S += '\n  Confinement date(s)  :' + str(self.listConfDates)
		S += '\n  Deconfinement date(s):' + str(self.listDeconfDates)
		S += '\n  Other date(s)        :' + str(self.listOtherDates)
		return S

	def setListFirstCaseDates(self, listFirstCaseDates):
		for date in listFirstCaseDates:
			self.addFirstCaseDates(date)
	def addFirstCaseDates(self, aStrDate):
		if type(aStrDate) == str:
			try: 
				parse(aStrDate, fuzzy=False)
				d1 = dt.strptime(aStrDate, '%Y-%m-%d')
				self.listFirstCaseDates.append(aStrDate)
			except ValueError:
				print('-->This string ', aStrDate, ' is not a date')

	def setListConfDates(self, listConfDates):
		for date in listConfDates:
			self.addConfDates(date)
	def addConfDates(self, aStrDate):
		if type(aStrDate) == str:
			try: 
				parse(aStrDate, fuzzy=False)
				d1 = dt.strptime(aStrDate, '%Y-%m-%d')
				self.listConfDates.append(aStrDate)
			except ValueError:
				print('-->This string ', aStrDate, ' is not a date')

	def setListDeconfDates(self, listDeconfDates):
		for date in listDeconfDates:
			self.addDeconfDates(date)
	def addDeconfDates(self, aStrDate):
		if type(aStrDate) == str:
			try: 
				parse(aStrDate, fuzzy=False)
				d1 = dt.strptime(aStrDate, '%Y-%m-%d')
				self.listDeconfDates.append(aStrDate)
			except ValueError:
				print('-->This string ', aStrDate, ' is not a date')

	def setListOtherDates(self, listOtherDates):
		for date in listOtherDates:
			self.addOtherDates(date)
	def addOtherDates(self, aStrDate):
		if type(aStrDate) == str:	
			try: 
				parse(aStrDate, fuzzy=False)
				d1 = dt.strptime(aStrDate, '%Y-%m-%d')
				self.listOtherDates.append(aStrDate)
			except ValueError:
				print('-->This string ', aStrDate, ' is not a date')

