#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
			], \
			[['R84+'], ['R27+'], ['R53+'], ['R24+'], ['R94+'], ['R44+'], ['R32+'], ['R11+'], ['R28+'], ['R75+'], ['R76+'], ['R52+'], ['R93+']]
	if element=='MetropoleD': 
		listdpts = []
		for i in range(1,96):
			number_str = str(i)
			zero_filled_number = number_str.zfill(2)
			listdpts.append(['D'+zero_filled_number])
		listdpts.remove(['D20']) # Ce département n'est pas dans les données
		listdpts.append(['D2A'])
		listdpts.append(['D2B'])
		return listdpts, listdpts

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
