#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np

from datetime            import datetime, timedelta
from sklearn.metrics     import mean_squared_error

from common              import readDataEurope, getDates, Plot, addDaystoStrDate, getLowerDateFromString, getNbDaysBetweenDateFromString
from SolveDiff_SEIR1R2   import SolveDiff_SEIR1R2

# constante
fileLocalCopy    = True  # if we upload the file from the url (to get latest results) or from a local copy file
readStartDateStr = "2020-03-01" 
readStopDateStr  = None

if __name__ == '__main__':

    verbose     = 2
    plot        = True
    country     = 'France' #United_Kingdom
    DatesString = getDates(country, verbose)

    # These are the date of confinement and deconfinement + other. See function getDates on how to add or delete dates to put the focus on
    fitStartDate = addDaystoStrDate(DatesString.listConfDates[0],   10) # ajout de 10 jours à la date de confinement
    fitStopDate  = addDaystoStrDate(DatesString.listDeconfDates[0], 10) # ajout de 10 jours à la date de déconfinement

    # Lecture des données et copy of the observations
    #############################################################################
    pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = \
                readDataEurope(country, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
    readStartDate = datetime.strptime(readStartDateStr, "%Y-%m-%d")
    readStopDate  = datetime.strptime(readStopDateStr,  "%Y-%m-%d")
    fitStopDate   = getLowerDateFromString(readStopDateStr, fitStopDate)
    dataLength    = pd_exerpt.shape[0]
    
    if verbose>1:
        print(readStartDateStr, readStopDateStr)
        print(readStartDate, readStopDate)
        print('dataLength=', dataLength)
        #input('pause')
    
    # Modele d'eq. diff non lineaires
    #############################################################################
    # Solveur eq. diff.
    dt      = 1
    solveur = SolveDiff_SEIR1R2(N, dt, verbose)
    #solveur.setParamInit(N=N, E0=0, I0=N*0.005, R10=0, R20=0) # 1% de la population infectée
    solveur.setParamInit(N=N, E0=183062, I0=86795, R10=5692, R20=43778) # paramètre de la précédente pahse
    solveur.modele.setParam(N=N, a=0.100, b=0.20, c=0.15, f=0.19)       # normalement f=0.115
    if verbose>0:
        print('solveur apres=', solveur)

    # Récupérations des données observées
    nbdata = getNbDaysBetweenDateFromString(fitStartDate, fitStopDate)
    timefocusdata = slice(0, nbdata)
    data   = np.zeros(shape=(nbdata))
    for i, z in enumerate(pd_exerpt.loc[fitStartDate:addDaystoStrDate(fitStopDate, -1), (HeadData[0])]):
        data[i]=z

    prefFig = './figures/' + solveur.modele.modelShortName + '_' + country

    # Integration time grid
    T    = 200
    time = np.linspace(0, T-1, T)

    # Solve ode avant optimization
    sol_ode = solveur.solve_SEIR1R2(time)
    # calcul time shift initial (ts0) with respect to data
    ts0 = solveur.compute_tsfromEQM(data, T)
    if verbose>0:
        print(solveur)
        print('  ts0='+str(ts0))

    # Panda dataframe
    columns = ['$S(t)$', '$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$', '$R(t)=R^1(t)+R^2(t)$']

    # plot
    if plot==True:
        listePlot    = [3]
        filename     = prefFig + '_Fitinit_' + ''.join(map(str, listePlot)) + '.png'
        timefocusedo = slice(ts0, ts0+nbdata+5)
        solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data, text=solveur.getTextParam())
        listePlot    = [1,2,3]
        filename     = prefFig + '_Fitinit_'+ ''.join(map(str, listePlot)) + '.png'
        timefocusedo = slice(ts0, ts0+nbdata+5)
        solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data, text=solveur.getTextParam())

        # begDate = addDaystoStrDate(fitStartDate, -int(ts0))
        # #print('begDate=', begDate)
        # for j in range(5):
        #     pd_exerpt.loc[begDate:, (columns[j])]=solut_edo[0:nbdata+ts0, j]
        # pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(solut_edo[0:nbdata+ts0, 0:4], axis=1)
        # pd_exerpt.loc[begDate:, (columns[6])] = np.sum(solut_edo[0:nbdata+ts0, 3:5], axis=1)
        # if verbose>1:
        #     print(pd_exerpt.tail())

        # titre = country + ' - ' + solveur.modele.modelName
        # # Plot de E, I, R^1, R^2
        # Plot(pd_exerpt, titre, prefFig+'_FitInit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)


    # Parameters optimization
    ############################################################################

    solveur.paramOptimization(data, time)
    if verbose>0:
        print('Solver''s state after optimization=', solveur)

     # After optimization
    ###############################
    
    # Solve equa diff apres optimization
    sol_ode = solveur.solve_SEIR1R2(time)
    # calcul time shift  (ts) with respect to data
    ts = solveur.compute_tsfromEQM(data, T)
    if verbose>0:
        print(solveur)
        print('  ts='+str(ts))

    solut_edo = solveur.solve_SEIR1R2(time)
    if plot==True:
        listePlot    = [3]
        filename     = prefFig + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
        timefocusedo = slice(ts, ts+nbdata+5)
        solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data, text=solveur.getTextParam())

        listePlot    = [1,2,3]
        filename     = prefFig + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
        timefocusedo = slice(ts, ts+nbdata+5)
        solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data, text=solveur.getTextParam())


        # begDate = addDaystoStrDate(fitStartDate, -int(ts))
        # #print('begDate=', begDate)
        # for j in range(5):
        #     pd_exerpt.loc[begDate:, (columns[j])]=solut_edo[0:nbdata+ts, j]
        # pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(solut_edo[0:nbdata+ts, 0:4], axis=1)
        # pd_exerpt.loc[begDate:, (columns[6])] = np.sum(solut_edo[0:nbdata+ts, 3:5], axis=1)
        # if verbose>1:
        #     print(pd_exerpt.tail())

        # titre = country + ' - ' + solveur.modele.modelName
        # # Plot de E, I, R^1, R^2
        # Plot(pd_exerpt, titre, prefFig+'_Fit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)


    
   