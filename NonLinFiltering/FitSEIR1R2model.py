#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
from sklearn.metrics import mean_squared_error

from common              import readDataEurope, getDates, addDaystoStrDate
from SolveDiff_SEIR1R2   import SolveDiff_SEIR1R2

# constante
fileLocalCopy = False  # if we upload the file from the url (to get latest results) or from a local copy file
startDate     = None 
stopDate      = None 

if __name__ == '__main__':

    verbose = 1
    plot    = True
    country, f = 'United_Kingdom', 0.45 #0.007 
    country, f = 'Italy', 0.11 #0.007 
    #country, f = 'Spain', 0.14 #0.007 
    country, f = 'France',         0.064 # 0.0032 qd R10 = 0. ; 0.093 quand on met R10 = zs[0].item()

    # These are the date of confinement and deconfinement + other. See function getDates on how to add or delete dates to put the focus on
    Dates     = getDates(country, verbose)
    startDate = addDaystoStrDate(Dates.listConfDates[0],   10) # ajout de 10 jours à la date de confinement
    stopDate  = addDaystoStrDate(Dates.listDeconfDates[0], 10) # ajout de 10 jours à la date de déconfinement
    

    # Lecture des données et copy of the observation
    #############################################################################
    pd, z_observ, N = readDataEurope(country=country, dateMin=startDate, dateMax=stopDate, \
                            plot=plot, fileLocalCopy=fileLocalCopy, verbose=verbose)
    if pd.empty==True:
        print('pd is empty - No data! --> exit')
        exit(1)

    zs = []
    for z in pd[z_observ[0]]:
        zs.append(np.array([z]))

    # Modele d'eq. diff non lineaires
    #############################################################################
    # Solveur eq. diff.
    dt      = 1
    solveur = SolveDiff_SEIR1R2(N, dt, verbose)
    solveur.set_paramf(f)
    solveur.setParamInit(N, 0, N*0.01, 0, 0) # 1% de la population infectée
    
    prefFig = './figures/' + solveur.modele.modelShortName + '_' + country

    # Solveur
    ############################################################################
    
    # Integration time grid
    simulLenght = 80
    vectTime    = np.linspace(0, simulLenght-1, simulLenght)

    # equa diff with the initial values
    solut = solveur.solve_SEIR1R2_sol1(vectTime)
    print(solveur)

    # calcul timeshift initial (ts0)
    eqm=[]
    for t in range(0, simulLenght-len(zs)):
        eqm.append(mean_squared_error(solut[t:t+len(zs), 3], zs))
    ts0 = np.argmin(eqm)

    # plot
    if plot==True:
        listePlot=[3]
        print('ts0=', ts0)
        solveur.plot_SEIR1R2(prefFig+'_Fitinit.png', vectTime, plot=listePlot, zs=(zs, ts0))

    # Parameters optimization
    ############################################################################

    ts = solveur.paramOptimization(ts0, zs, vectTime)
    
    solveur.solve_SEIR1R2_sol1(vectTime)
    if plot==True:
        listePlot=[3]
        solveur.plot_SEIR1R2(prefFig+'_Fit.png', vectTime, plot=listePlot, zs=(zs, ts))

    
   