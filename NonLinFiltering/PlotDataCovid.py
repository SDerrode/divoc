#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from datetime  import datetime, timedelta
from common    import readDataEurope, readDataFrance, getDates, PlotData, addDaystoStrDate

# constante
fileLocalCopy = False         # if we upload the file from the url (to get latest results) or from a local copy file

if __name__ == '__main__':
    """
        Program to perform UKF filtering on Covid Data.
 
        :Example:

        For country (European database)
        >> python3 PlotDataCovid.py 
        >> python3 PlotDataCovid.py United_Kingdom
        >> python3 PlotDataCovid.py Italy 1 
        >> python3 PlotDataCovid.py France,Germany 1 # Shortcut for processing the two countries successively
        >> python3 PlotDataCovid.py France,Spain,Italy,United_Kingdom,Germany,Belgium 0

        For French Regions (French database)
        >> python3 ProcessSEIR1R2.py France,69       # Dpt 69 (Rhône)
        >> python3 ProcessSEIR1R2.py France,69,75,01 # Dpt 69 (Rhône) + Dpt 75 + Dpt 01
       
        argv[1] : Country (or list separeted by ','), or 'France' followed by a list of departments. Default: France 
        argv[2] : Verbose level (debug: 3, ..., almost mute: 0).                                     Default: 1
    """

    print('Command line : ', sys.argv, flush=True)
    if len(sys.argv) > 3:
        print('  CAUTION : bad number of arguments - see help')
        exit(1)

    # Default value for parameters
    listplaces = ['France']
    verbose    = 2
    
    # Parameters from argv
    if len(sys.argv)>1: listplaces = list(sys.argv[1].split(','))
    FrDatabase = False
    if listplaces[0]=='France' and len(listplaces)>1:
        try:
            int(listplaces[1])
        except Exception as e:
            FrDatabase=False
        else:
            FrDatabase = True
            listplaces = listplaces[1:]
            France     = 'France'

    if len(sys.argv)>2: verbose = int(sys.argv[2])

    # Constantes
    dt = 1
    readStartDateStr = "2020-03-01" 
    readStopDateStr  = None
    
    # Loop for all places
    for place in listplaces:

        placefull = place
        # These are the date of confinement and deconfinement + other. 
        # See function getDates on how to add or delete dates to put the focus on
        DatesString = getDates(place, verbose)
        if FrDatabase==True: 
            placefull   = France + place
            DatesString = getDates(France, verbose)

        print('PROCESSING of', placefull, 'in', listplaces)

        # Constantes
        import os
        repertoire = './figures/'+ placefull
        if not os.path.exists(repertoire):
            os.makedirs(repertoire)
        prefFig = repertoire + '/Data_' + placefull
        
        # Lecture des données et copy of the observation
        #############################################################################
        if FrDatabase==True:
            pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
        else:
            pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)

        readStartDate = datetime.strptime(readStartDateStr, "%Y-%m-%d")
        readStopDate  = datetime.strptime(readStopDateStr,  "%Y-%m-%d")
        dataLength = pd_exerpt.shape[0]
        if verbose>1:
            print('readStartDateStr=', readStartDateStr, ', readStopDateStr=', readStopDateStr)
            print('readStartDate   =', readStartDate,    ', readStopDate   =', readStopDate)
            #input('pause')

        # On ajoute le gradient
        pd_exerpt['Instant cases']  = pd_exerpt[HeadData[0]].diff()
        pd_exerpt['Instant deaths'] = pd_exerpt[HeadData[1]].diff()
        # print('Head=', pd_exerpt.head())
        # print('tail=', pd_exerpt.tail())
        # print('HeadData=', HeadData)
        # print('liste=', list(pd_exerpt))

        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'_'    +HeadData[0]+'.png', y=HeadData[0],        color='green', Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'_'    +HeadData[1]+'.png', y=HeadData[1],        color='red',   Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'_Diff'+HeadData[0]+'.png', y=['Instant cases'],  color='green', Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'_Diff'+HeadData[1]+'.png', y=['Instant deaths'], color='red',   Dates=DatesString)
