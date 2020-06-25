#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from datetime         import datetime, timedelta

from filterpy.kalman  import UnscentedKalmanFilter as UKF
from filterpy.kalman  import JulierSigmaPoints, MerweScaledSigmaPoints, rts_smoother
from filterpy.common  import Q_discrete_white_noise

from common           import readDataEurope, readDataFrance, readDates, PlotData
from common           import getRepertoire, addDaystoStrDate
from France           import getPlace

# constante
fileLocalCopy = False # if we upload the file from the url (to get latest results) or from a local copy file


def fR1(r1, dt):
    return r1 # on renvoie R1
def hR1(r1):
    return r1 # on renvoie R1
    
def fF(f, dt):
    return f # on renvoie F
def hF(f):
    return f # on renvoie F

def fR1F(r1f, dt):
    return r1f # on renvoie R1 et F
def hR1F(r1f):
    return r1f # on renvoie R1 et F

def main(sysargv):
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
        >> python3 PlotDataCovid.py FRANCE,D69         # Code Insee Dpt 69 (Rhône)
        >> python3 PlotDataCovid.py FRANCE,R84         # Tous les dpts de la Région dont le code Insee est 84 (Auvergne-Rhone-Alpes)
        >> python3 PlotDataCovid.py FRANCE,R32+        # Somme de tous les dpts de la Région 32 (Hauts de France)
        >> python3 PlotDataCovid.py FRANCE,MetropoleD  # Tous les départements de la France métropolitaies
        >> python3 PlotDataCovid.py FRANCE,MetropoleR+ # Somme des dpts de toutes les régions françaises
        Toute combinaison est possible : exemple FRANCE,R32+,D05,R84

        argv[1] : List of countries (ex. France,Germany,Italy), or see above. Default: France 
        argv[2] : Verbose level (debug: 3, ..., almost mute: 0).              Default: 1
    """

    #Austria,Belgium,Croatia,Czechia,Finland,France,Germany,Greece,Hungary,Ireland,Italy,Lithuania,Poland,Portugal,Romania,Serbia,Spain,Switzerland,Ukraine
    

    print('Command line : ', sysargv, flush=True)
    if len(sysargv) > 3:
        print('  CAUTION : bad number of arguments - see help')
        exit(1)

    # Default value for parameters
    listplaces = ['France']
    verbose    = 1
    
    # Parameters from argv
    if len(sysargv)>1: liste   = list(sysargv[1].split(','))
    if len(sysargv)>2: verbose = int(sysargv[2])

    listplaces = []
    listnames = []
    if liste[0]=='FRANCE':
        FrDatabase = True
        liste = liste[1:]
        for el in liste:
            l,n=getPlace(el)
            if el=='MetropoleR+':
                for l1,n1 in zip(l,n):
                    listplaces.extend(l1)
                    listnames.extend([n1])
            else:
                listplaces.extend(l)
                listnames.extend(n)
    else:
        listplaces = liste[:]
        FrDatabase = False

    # print('liste=', liste)
    # print('listplaces=', listplaces)
    # print('listnames=', listnames)
    # #input(';sjbvcbj')

    # Constantes
    dt = 1
    readStartDateStr = "2020-03-01" 
    readStopDateStr  = None
    
    # Loop for all places
    for index,place in enumerate(listplaces):
        # print('index,place=', index,place)
        # print('listnames=', listnames)

        # Get the full name of the place to process, and the special dates corresponding to the place
        France='France'
        if FrDatabase == True: 
            placefull   = 'France-' + listnames[index][0]
            DatesString = readDates(France, verbose)
        else:
            placefull   = place
            DatesString = readDates(place, verbose)

        print('PROCESSING of', placefull, 'in', listnames)
        
        # Repertoire figures
        repertoire = getRepertoire(True, './figures/data/'+placefull)
        prefFig = repertoire + '/'
        
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
        pd_exerpt['Diff ' + HeadData[0]] = pd_exerpt[HeadData[0]].diff()
        pd_exerpt['Diff ' + HeadData[1]] = pd_exerpt[HeadData[1]].diff()
        pd_exerpt['Diff ' + HeadData[2]] = pd_exerpt[HeadData[2]].diff()
        # print('Head=', pd_exerpt.head())
        # print('tail=', pd_exerpt.tail())
        # print('HeadData=', HeadData)
        # print('liste=', list(pd_exerpt))

        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+       HeadData[0]+'.png',                         y=HeadData[0],                 color='red',   Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+       HeadData[1]+'.png',                         y=HeadData[1],                 color='black', Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+       HeadData[2]+'.png',                         y=HeadData[2],                 color='black', Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+       HeadData[0]+HeadData[1]+'.png',             y=[HeadData[0], HeadData[1]],  color=['red', 'black'], Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'Diff'+HeadData[0]+'.png',                         y=['Diff ' + HeadData[0]], color='red',   Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'Diff'+HeadData[1]+'.png',                         y=['Diff ' + HeadData[1]], color='black', Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'Diff'+HeadData[0]+HeadData[1]+'.png',             y=['Diff ' + HeadData[0], 'Diff ' + HeadData[1]], color=['red', 'black'], Dates=DatesString)
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'Diff'+HeadData[0]+HeadData[1]+HeadData[2]+'.png', y=['Diff ' + HeadData[0], 'Diff ' + HeadData[1], 'Diff ' + HeadData[2]], color=['red', 'black', 'blue'], Dates=DatesString)
        
        # on filtre R1plusD par UKF
        #############################################################################
        data   = pd_exerpt[HeadData[2]].tolist()
        dt     = 1
        sigmas = MerweScaledSigmaPoints(n=1, alpha=.5, beta=2., kappa=0.) #1-3.)
        ukf    = UKF(dim_x=1, dim_z=1, fx=fR1, hx=hR1, dt=dt, points=sigmas)
        # Filter init
        ukf.x[0] = data[0]
        ukf.Q    = np.diag([30.])
        ukf.R    = np.diag([170.])
        if verbose>1:
            print('ukf.x[0]=', ukf.x[0])
            print('ukf.R   =', ukf.R)
            print('ukf.Q   =', ukf.Q)
            
        # UKF filtering and smoothing, batch mode
        R1filt, _ = ukf.batch_filter(data)

        # plotting
        pd_exerpt[HeadData[2] + ' filt'] = R1filt
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'filt'+HeadData[2]+'.png', y=[HeadData[2], HeadData[2] + ' filt'], color=['red', 'darkred'], Dates=DatesString)

        pd_exerpt['Diff ' + HeadData[2] + ' filt']  = pd_exerpt[HeadData[2] + ' filt'].diff()
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'diff_filt'+HeadData[2]+'.png', y=['Diff ' + HeadData[2], 'Diff ' + HeadData[2] + ' filt'], color=['red', 'darkred'], Dates=DatesString)


        # on filtre Diff R1 par UKF
        # Ca marche mais identique au précédent
        #############################################################################
        # data    = pd_exerpt['Diff cases'].tolist()
        # data[0] = data[1]
        # print('data=', data)
        # dt     = 1
        # sigmas = MerweScaledSigmaPoints(n=1, alpha=.5, beta=2., kappa=1.) #1-3.)
        # ukf    = UKF(dim_x=1, dim_z=1, fx=fR1, hx=hR1, dt=dt, points=sigmas)
        # # Filter init
        # ukf.x[0] = data[0]
        # ukf.Q    = np.diag([30.])
        # ukf.R    = np.diag([170.])
        # if verbose>1:
        #     print('ukf.x[0]=', ukf.x[0])
        #     print('ukf.R   =', ukf.R)
        #     print('ukf.Q   =', ukf.Q)
            
        # # UKF filtering and smoothing, batch mode
        # diffR1filt, _ = ukf.batch_filter(data)
        # pd_exerpt['diffR1 filt'] = diffR1filt
        # PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'diffcases_filt'+HeadData[0]+'.png', y=['Diff cases', 'diffR1 filt'], color=['red', 'darkred'], Dates=DatesString)

        # on filtre F par UKF
        #############################################################################
        data   = pd_exerpt[HeadData[1]].tolist()
        dt     = 1
        sigmas = MerweScaledSigmaPoints(n=1, alpha=.5, beta=2., kappa=0.) #1-3.)
        ukf    = UKF(dim_x=1, dim_z=1, fx=fF, hx=hF, dt=dt, points=sigmas)
        # Filter init
        ukf.x[0] = data[0]
        ukf.Q    = np.diag([15.])
        ukf.R    = np.diag([100.])
        if verbose>1:
            print('ukf.x[0]=', ukf.x[0])
            print('ukf.R   =', ukf.R)
            print('ukf.Q   =', ukf.Q)
            
        # UKF filtering and smoothing, batch mode
        Ffilt, _ = ukf.batch_filter(data)

        # plotting
        pd_exerpt[HeadData[1] + ' filt'] = Ffilt
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'filt'+HeadData[1]+'.png', y=[HeadData[1], HeadData[1] + ' filt'], color=['gray', 'black'], Dates=DatesString)

        pd_exerpt['Diff ' + HeadData[1] + ' filt']  = pd_exerpt[HeadData[1] + ' filt'].diff()
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'diff_filt'+HeadData[1]+'.png', y=['Diff ' + HeadData[1], 'Diff ' + HeadData[1] + ' filt'], color=['gray', 'black'], Dates=DatesString)


        # on filtre R1 et F simultanément par UKF
        #############################################################################
        data   = pd_exerpt[[HeadData[0], HeadData[1]]].to_numpy(copy=True)
        dt     = 1
        sigmas = MerweScaledSigmaPoints(n=2, alpha=.5, beta=2., kappa=1.) #1-3.)
        ukf    = UKF(dim_x=2, dim_z=2, fx=fR1F, hx=hR1F, dt=dt, points=sigmas)
        # Filter init
        ukf.x[:] = data[0, :]
        ukf.Q    = np.diag([30., 15.])
        ukf.R    = np.diag([170., 100.])
        if verbose>1:
            print('ukf.x[:]=', ukf.x[:])
            print('ukf.R   =', ukf.R)
            print('ukf.Q   =', ukf.Q)
            
        # UKF filtering and smoothing, batch mode
        R1Ffilt, _ = ukf.batch_filter(data)

        # plotting
        pd_exerpt[HeadData[0]+' filtboth'] = R1Ffilt[:, 0]
        pd_exerpt[HeadData[1]+' filtboth'] = R1Ffilt[:, 1]
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'filtboth'+HeadData[0]+HeadData[1]+'.png', \
                y=[HeadData[0], HeadData[0]+' filtboth', HeadData[1], HeadData[1]+' filtboth'], color=['red', 'darkred', 'gray', 'black'], Dates=DatesString)

        pd_exerpt['Diff '+HeadData[0]+' filtboth']  = pd_exerpt[HeadData[0]+' filtboth'].diff()
        pd_exerpt['Diff '+HeadData[1]+' filtboth']  = pd_exerpt[HeadData[1]+' filtboth'].diff()
        PlotData(pd_exerpt, titre=placefull, filenameFig=prefFig+'diff_filt'+HeadData[0]+HeadData[1]+'.png', \
                y=['Diff '+HeadData[0], 'Diff '+HeadData[0]+' filtboth', 'Diff '+HeadData[1], 'Diff '+HeadData[1]+' filtboth', ], color=['red', 'darkred', 'gray', 'black'], Dates=DatesString)


if __name__ == '__main__':
    main(sys.argv)