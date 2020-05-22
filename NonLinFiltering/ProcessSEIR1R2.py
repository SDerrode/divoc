#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from datetime          import datetime, timedelta
from common            import readDataEurope, getDates, Plot, addDaystoStrDate, getLowerDateFromString, getNbDaysBetweenDateFromString
from SolveDiff_SEIR1R2 import SolveDiff_SEIR1R2

import matplotlib.pyplot as plt

# constante
fileLocalCopy    = True  # if we upload the file from the url (to get latest results) or from a local copy file
readStartDateStr = "2020-03-01" 
readStopDateStr  = None

if __name__ == '__main__':

    verbose     = 2
    plot        = True
    surplus     = False
    decalage    = 10       # 10 jours de decalage
    country     = 'France' # United_Kingdom
    DatesString = getDates(country, verbose)

    # Lecture des données et copy of the observations
    #############################################################################
    pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = \
                readDataEurope(country, readStartDateStr, readStopDateStr, plot, fileLocalCopy, verbose=0)
    readStartDate = datetime.strptime(readStartDateStr, "%Y-%m-%d")
    readStopDate  = datetime.strptime(readStopDateStr,  "%Y-%m-%d")
    dataLength = pd_exerpt.shape[0]
    if verbose>1:
        print(readStartDateStr, readStopDateStr)
        print(readStartDate, readStopDate)
        #input('pause')
    
    # Liste des dates pour les approximations piecewise
    fitStartDate = addDaystoStrDate(DatesString.listConfDates[0],   10) # ajout de 10 jours à la date de confinement
    fitStopDate  = addDaystoStrDate(DatesString.listDeconfDates[0], 10) # ajout de 10 jours à la date de déconfinement

    confin_decalage   = addDaystoStrDate(DatesString.listConfDates[0],   decalage)
    deconfin_decalage = addDaystoStrDate(DatesString.listDeconfDates[0], decalage)
    ListDatesStr = [readStartDateStr, confin_decalage, deconfin_decalage, readStopDateStr]
    ListDatesStr = [readStartDateStr, confin_decalage, deconfin_decalage, readStopDateStr]
    ListDates    = [datetime.strptime(dateString,"%Y-%m-%d") for dateString in ListDatesStr]
    if ListDates[-1]<ListDates[-2]:
        ListDatesStr.pop()
        ListDates.pop()
    if verbose>2:
        print('ListDates=', ListDates)
        print('ListDatesStr=', ListDatesStr)
        #input('pause')
    
    # Modele d'eq. diff non lineaires
    #############################################################################
    # Solveur eq. diff.
    dt      = 1
    solveur = SolveDiff_SEIR1R2(N, dt, verbose)

    # Constantes
    prefFig = './figures/' + solveur.modele.modelShortName + '_' + country
    columns = ['$S(t)$', '$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$', '$R(t)=R^1(t)+R^2(t)$']

    # Surplus de la précédente période (nulle pour le début)
    data_surplus   = np.zeros(shape=(dataLength))
    data_corrected = np.zeros(shape=(dataLength))
    indMaxPeriod   = 0

    # Boucle pour traiter successivement les différentes fenêtres
    #############################################################################
    for i in range(len(ListDatesStr)-1):

        # paramètres initiaux
        if i==0:
            E0, I0, R10, R20 = 0, 1, 0, 0
            l, b0, c0, f0    = 0.255, 1./5.2, 1./10, 0.10
            a0 = (l+c0)*(1.+l/b0)
            T  = 150
        if i==1:
            if surplus==True:
                a0, b0, c0, f0   = 0.100, 0.18, 0.09, 0.14 # en enlevant le surplus
            else:
                a0, b0, c0, f0   = 0.03, 0.14, 0.10, 0.10 # sans enlever le surplus
            T = 100
        if i==2:
            if surplus==True:
                a0, b0, c0, f0   = 0.100, 0.18, 0.09, 0.14 # en enlevant le surplus
            else:
                a0, b0, c0, f0   = 0.05, 0.10, 0.05, 0.04 # en enlevant le surplus

            T = 100
        #print('Reproductivité avant: ', a0/c0)
        time = np.linspace(0, T-1, T)

        # date 
        fitStartDate    = ListDates[i]
        fitStopDate     = ListDates[i+1]
        fitStartDateStr = ListDatesStr[i]
        fitStopDateStr  = ListDatesStr[i+1]

        # Récupérations des données observées
        indMinPeriod = indMaxPeriod
        dataLengthPeriod = 0
        for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), (HeadData[0])]):
            data_corrected[indMinPeriod+j] = z - data_surplus[indMinPeriod+j]
            dataLengthPeriod +=1
        indMaxPeriod += dataLengthPeriod
        timefocusdata = slice(indMinPeriod, indMaxPeriod)
        if verbose>1:
            print('dataLength      =', dataLength)
            print('indMinPeriod    =', indMinPeriod)
            print('indMaxPeriod    =', indMaxPeriod)
            print('dataLengthPeriod=', dataLengthPeriod)
            print('fitStartDateStr =', fitStartDateStr)
            print('fitStopDateStr  =', fitStopDateStr)


        # Set initialisation data for the solveur
        solveur.modele.setParam(N=N, a=a0, b=b0, c=c0, f=f0)
        solveur.setParamInit(N=N, E0=E0, I0=I0, R10=R10, R20=R20)

        # Before optimization
        ###############################@

        # Solve ode avant optimization
        sol_ode = solveur.solve_SEIR1R2(time)
        # calcul time shift initial (ts0) with respect to data
        ts0 = solveur.compute_tsfromEQM(data_corrected[indMinPeriod:indMaxPeriod], T)
        if verbose>0:
            print(solveur)
            print('  ts0=', ts0)
        
        # plot
        if plot==True:
            listePlot    = [3]
            filename     = prefFig + '_Period' + str(i) + '_Fitinit_' + ''.join(map(str, listePlot)) + '.png'
            timefocusedo = slice(ts0, ts0+dataLengthPeriod+5)
            solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam())

        # Parameters optimization
        ############################################################################

        solveur.paramOptimization(data_corrected[indMinPeriod:indMaxPeriod], time)
        if verbose>0:
            print('Solver''s state after optimization=', solveur)
            _, a1, b1, c1, f1 = solveur.modele.getParam()
            print('  Reproductivité après: ', a1/c1)


        # After optimization
        ###############################
        
        # Solve equa diff apres optimization
        sol_ode = solveur.solve_SEIR1R2(time)
        # calcul time shift initial (ts0) with respect to data
        ts = solveur.compute_tsfromEQM(data_corrected[indMinPeriod:indMaxPeriod], T)
        if verbose>0:
            print(solveur)
            print('  ts=', ts)

        timefocusedo = slice(ts, min(ts+dataLengthPeriod+5, T))
        
        if plot==True:
            listePlot=[1,2,3]
            filename = prefFig + '_Period' + str(i) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
            solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam())
            listePlot=[3]
            filename = prefFig + '_Period' + str(i) + '_Fit_' + ''.join(map(str, listePlot)) + '.png'
            solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam())

        sol_ode_withSwitch = solveur.solve_SEIR1R2_withSwitch(T, timeswitch=ts+dataLengthPeriod)
        if plot==True:
            listePlot=[1,2,3]
            filename = prefFig + '_Period' + str(i) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
            solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam())
            listePlot=[3]
            filename = prefFig + '_Period' + str(i) + '_Fi_withSwitch_' + ''.join(map(str, listePlot)) + '.png'
            solveur.plot_SEIR1R2(filename, timefocusedo, timefocusdata, plot=listePlot, data=data_corrected, text=solveur.getTextParam())

            # begDate = addDaystoStrDate(fitStartDate, -int(ts))
            # #print('begDate=', begDate)
            # for j in range(5):
            #     pd_exerpt.loc[begDate:, (columns[j])]=sol_ode_withSwitch[0:dataLengthPeriod+ts, j]
            # pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(sol_ode_withSwitch[0:dataLengthPeriod+ts, 0:4], axis=1)
            # pd_exerpt.loc[begDate:, (columns[6])] = np.sum(sol_ode_withSwitch[0:dataLengthPeriod+ts, 3:5], axis=1)
            # if verbose>1:
            #     print(pd_exerpt.tail())

            # titre = country + ' - ' + solveur.modele.modelName
            # # Plot de E, I, R^1, R^2
            # Plot(pd_exerpt, titre, prefFig+'_Fit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)
        
        # nouvelles valeurs d'init de l'edo
        S0, E0, I0, R10, R20 = map(int, sol_ode_withSwitch[ts+dataLengthPeriod, :])
        
        if surplus == True:
            # get the residuals
            bout  = min(T, ts+dataLengthPeriod +(dataLength-indMaxPeriod))
            debut = ts+dataLengthPeriod
            # print('bout=', bout)
            # print('debut=', debut)
            # print('data_corrected[indMaxPeriod-1]=', data_corrected[indMaxPeriod-1])
            data_surplus[indMaxPeriod:indMaxPeriod+(bout-debut)] += sol_ode_withSwitch[debut:bout, 3]-data_corrected[indMaxPeriod-1]
        
        

        if surplus == True:
            fig = plt.figure(facecolor='w', figsize=(8, 4))
            ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
            ax.plot(data_surplus, label='data_surplus')
            plt.legend()
            #plt.show()
            plt.savefig(prefFig + '_Period' + str(i) + '_surplus.png')
            #plt.savefig(prefFig + 'surplus.png')
            plt.close()

        #input('attente')

