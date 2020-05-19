#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
from sklearn.metrics import mean_squared_error

from common              import readDataEurope, getDates, Plot, addDaystoStrDate, getLowerDateFromString, getNbDaysBetweenDateFromString
from SolveDiff_SEIR1R2   import SolveDiff_SEIR1R2

# constante
fileLocalCopy = True  # if we upload the file from the url (to get latest results) or from a local copy file
startDate     = "2020-03-01" 
stopDate      = None 

if __name__ == '__main__':

    verbose = 1
    plot    = True
    #country, f = 'United_Kingdom', 0.45 #0.007 
    #country, f = 'Italy', 0.11 #0.007 
    #country, f = 'Spain', 0.14 #0.007 
    country, f = 'France',         0.032 # 0.0032 qd R10 = 0. ; 0.093 quand on met R10 = zs[0].item()

    # These are the date of confinement and deconfinement + other. See function getDates on how to add or delete dates to put the focus on
    Dates        = getDates(country, verbose)
    fitStartDate = addDaystoStrDate(Dates.listConfDates[0],   10) # ajout de 10 jours à la date de confinement
    fitStopDate  = addDaystoStrDate(Dates.listDeconfDates[0], 10) # ajout de 10 jours à la date de déconfinement

    # Lecture des données et copy of the observations
    #############################################################################
    pd_exerpt, z_observ, N = readDataEurope(country=country, dateMin=startDate, dateMax=stopDate, \
                            plot=plot, fileLocalCopy=fileLocalCopy, verbose=verbose)

    # stopDates is corrected to fit the data dates
    stopDate    = pd_exerpt.index[-1].strftime("%Y-%m-%d")
    fitStopDate = getLowerDateFromString(stopDate, fitStopDate)

    # Récupérations des données observées
    zs = []
    for z in pd_exerpt.loc[fitStartDate:fitStopDate, (z_observ[0])]:
        zs.append(np.array([z]))
    # On remet le premier à 0
    R1_shift=zs[0]
    zs -= R1_shift
    # print('zs=', zs)
    # input('temp')

    # Modele d'eq. diff non lineaires
    #############################################################################
    # Solveur eq. diff.
    dt      = 1
    solveur = SolveDiff_SEIR1R2(N, dt, verbose)
    solveur.set_paramf(f)
    solveur.setParamInit(N=N, E0=0, I0=N*0.01, R10=0, R20=0) # 1% de la population infectée
    
    prefFig = './figures/' + solveur.modele.modelShortName + '_' + country

    # Solveur
    ############################################################################
    
    # Integration time grid
    simulLenght = 130
    vectTime    = np.linspace(0, simulLenght-1, simulLenght)

    # equa diff with the initial values
    solut_eqdiff = solveur.solve_SEIR1R2_sol1(vectTime)
    #print(solveur)

    # calcul timeshift initial (ts0)
    eqm=[]
    for t in range(0, simulLenght-len(zs)):
        eqm.append(mean_squared_error(solut_eqdiff[t:t+len(zs), 3], zs))
    ts0 = np.argmin(eqm)

    # Panda dataframe
    nbdata  = getNbDaysBetweenDateFromString(fitStartDate, fitStopDate)+1
    columns = ['$S(t)$', '$E(t)$', '$I(t)$', '$R^1(t)$', '$R^2(t)$', '$R^2(t)=N-\sum(SEIR^1)$', '$R(t)=R^1(t)+R^2(t)$']

    # plot
    if plot==True:
        listePlot=[3]
        print('ts0=', ts0)
        solveur.plot_SEIR1R2(prefFig+'_Fitinit.png', vectTime, plot=listePlot, zs=(zs, ts0))

        begDate = addDaystoStrDate(fitStartDate, -int(ts0))
        #print('begDate=', begDate)
        for j in range(5):
            pd_exerpt.loc[begDate:, (columns[j])]=solut_eqdiff[0:nbdata+ts0, j]
        pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(solut_eqdiff[0:nbdata+ts0, 0:4], axis=1)
        pd_exerpt.loc[begDate:, (columns[6])] = np.sum(solut_eqdiff[0:nbdata+ts0, 3:5], axis=1)
        if verbose>1:
            print(pd_exerpt.tail())

        titre = country + ' - ' + solveur.modele.modelName
        # Plot de E, I, R^1, R^2
        Plot(pd_exerpt, titre, prefFig+'_FitInit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)


    # Parameters optimization
    ############################################################################

    ts = solveur.paramOptimization(ts0, zs, vectTime)
    if verbose>0:
        print('Solver''s state after optimization=', solveur)

    
    solut_eqdiff = solveur.solve_SEIR1R2_sol1(vectTime)
    if plot==True:
        listePlot=[3]
        solveur.plot_SEIR1R2(prefFig+'_Fit.png', vectTime, plot=listePlot, zs=(zs, ts))

        begDate = addDaystoStrDate(fitStartDate, -int(ts))
        #print('begDate=', begDate)
        for j in range(5):
            pd_exerpt.loc[begDate:, (columns[j])]=solut_eqdiff[0:nbdata+ts, j]
        pd_exerpt.loc[begDate:, (columns[5])] = solveur.modele.N-np.sum(solut_eqdiff[0:nbdata+ts, 0:4], axis=1)
        pd_exerpt.loc[begDate:, (columns[6])] = np.sum(solut_eqdiff[0:nbdata+ts, 3:5], axis=1)
        if verbose>1:
            print(pd_exerpt.tail())

        titre = country + ' - ' + solveur.modele.modelName
        # Plot de E, I, R^1, R^2
        Plot(pd_exerpt, titre, prefFig+'_Fit_EIR1R2.png', solveur.modele, y=[3], Dates=Dates, z_observ=z_observ)


    
   