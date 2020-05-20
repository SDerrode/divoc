#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from datetime        import datetime, timedelta
from sklearn.metrics import mean_squared_error

from common            import readDataEurope, getDates, Plot, addDaystoStrDate, getLowerDateFromString, getNbDaysBetweenDateFromString
from SolveDiff_SEIR1R2 import SolveDiff_SEIR1R2

# constante
fileLocalCopy    = True  # if we upload the file from the url (to get latest results) or from a local copy file
readStartDateStr = "2020-03-01" 
readStopDateStr  = None 

if __name__ == '__main__':

    verbose       = 1
    plot          = True
    country, tabf = 'United_Kingdom', [0.032, 0.6, 0.8]
    DatesString   = getDates(country, verbose)

    # paramètres initiaux
    E0, I0, R10, R20 = 0, 1, 0, 0
    #a0, b0, c0, f0   = 0.155, 1./5.2, 1./12.39, 0.02

    l, b0, c0, f0 = 0.255, 1./5.2, 1./12.68, 0.08
    a0 = (l+c0)*(1.+l/b0)
    print('Reproductivité avant: ', a0/c0)

    # Lecture des données et copy of the observations
    #############################################################################
    pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = \
                readDataEurope(country, readStartDateStr, readStopDateStr, plot, fileLocalCopy, verbose=0)
    readStartDate = datetime.strptime(readStartDateStr, "%Y-%m-%d")
    readStopDate  = datetime.strptime(readStopDateStr,  "%Y-%m-%d")
    if verbose>1:
        print(readStartDateStr, readStopDateStr)
        print(readStartDate, readStopDate)
        #input('pause')
    
    # Liste des dates pour les approximations piecewise
    ListDatesStr = [readStartDateStr, DatesString.listConfDates[0], DatesString.listDeconfDates[0], readStopDateStr]
    ListDates    = [datetime.strptime(dateString,"%Y-%m-%d") for dateString in ListDatesStr]
    if verbose>0:
        print('ListDates=', ListDates)
        print('ListDatesStr=', ListDatesStr)
    
    # Modele d'eq. diff non lineaires
    #############################################################################
    # Solveur eq. diff.
    dt      = 1
    solveur = SolveDiff_SEIR1R2(N, dt, verbose)

    # Constantes
    prefFig = './figures/' + solveur.modele.modelShortName + '_' + country
    columns = ['$S(t)$', '$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$', '$R(t)=R^1(t)+R^2(t)$']


    # Boucle pour traiter successivement les différentes fenêtres
    #############################################################################
    for i in range(len(ListDatesStr)-1):
        fitStartDate    = ListDates[i]
        fitStopDate     = ListDates[i+1]
        fitStartDateStr = ListDatesStr[i]
        fitStopDateStr  = ListDatesStr[i+1]

        # Récupérations des données observées
        data = []
        for z in pd_exerpt.loc[fitStartDateStr:fitStopDateStr, (HeadData[0])]:
            #data.append(np.array([z]))
            data.append(z)
        dataLength = len(data)
        if verbose>0:
            print('dataLength=', dataLength)

        # Set initialisation data for the solveur
        solveur.modele.setParam(N=N, a=a0, b=b0, c=c0, f=f0)
        solveur.setParamInit(N=N, E0=E0, I0=I0, R10=R10, R20=R20) #data[0].item()*(1.-f0)/f0)
    
        # Integration time grid for eq. diff.
        T    = 350
        time = np.linspace(0, T-1, T)


        # Before optimization
        ###############################@

        # Solve equa diff avant optimization
        solut_eqdiff = solveur.solve_SEIR1R2_sol1(time)
        if verbose>0:
            print(solveur)

        # calcul time shift initial (ts0)
        eqm=np.zeros(shape=(T-dataLength))
        for t in range(T-dataLength):
            eqm[t] = mean_squared_error(solut_eqdiff[t:t+dataLength, 3], data)
        ts0 = np.argmin(eqm)
        if verbose>0:
            print('ts0=', ts0)

        # plot
        if plot==True:
            listePlot=[3]
            filename = prefFig + '_Period' + str(i) + '_Fitinit.png'
            timefocus = slice(ts0, ts0+dataLength)
            solveur.plot_SEIR1R2(filename, timefocus, plot=listePlot, data=data)

            # begDate = (fitStartDate+ timedelta(-int(ts0))).strftime("%Y-%m-%d")
            # print('begDate=', begDate)
            # for j in range(5):
            #     pd_exerpt.loc[begDate:, (columns[j])]=solut_eqdiff[0:dataLength+ts0, j]
            # pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(solut_eqdiff[0:dataLength+ts0, 0:4], axis=1)
            # pd_exerpt.loc[begDate:, (columns[6])] = np.sum(solut_eqdiff[0:dataLength+ts0, 3:5], axis=1)
            # if verbose>1:
            #     print(pd_exerpt.tail())

            # titre = country + ' - ' + solveur.modele.modelName + ' - Period ' + str(i)
            # # Plot de E, I, R^1, R^2
            # Plot(pd_exerpt, titre, prefFig + '_Period' + str(i) + '_FitInit_EIR1R2.png', solveur.modele, y=[3], Dates=DatesString, data=data)


        # Parameters optimization
        ############################################################################

        ts = solveur.paramOptimization(ts0, data, time)
        if verbose>0:
            print('Solver''s state after optimization=', solveur)
        _, a1, b1, c1, f1 = solveur.modele.getParam()
        print('Reproductivité après: ', a1/c1)


        # After optimization
        ###############################
        
        # Solve equa diff apres optimization
        solut_eqdiff = solveur.solve_SEIR1R2_sol1(time)

        # plot
        if plot==True:
            listePlot=[3]
            filename = prefFig + '_Period' + str(i) + '_Fit.png'
            timefocus = slice(ts, ts+dataLength)
            solveur.plot_SEIR1R2(filename, timefocus, plot=listePlot, data=data)

            # begDate = addDaystoStrDate(fitStartDate, -int(ts))
            # #print('begDate=', begDate)
            # for j in range(5):
            #     pd_exerpt.loc[begDate:, (columns[j])]=solut_eqdiff[0:dataLength+ts, j]
            # pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(solut_eqdiff[0:dataLength+ts, 0:4], axis=1)
            # pd_exerpt.loc[begDate:, (columns[6])] = np.sum(solut_eqdiff[0:dataLength+ts, 3:5], axis=1)
            # if verbose>1:
            #     print(pd_exerpt.tail())

            # titre = country + ' - ' + solveur.modele.modelName
            # # Plot de E, I, R^1, R^2
            # Plot(pd_exerpt, titre, prefFig+'_Fit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)

        exit()