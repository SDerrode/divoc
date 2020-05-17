#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
from lmfit               import minimize, Parameters, Parameter, report_fit
from scipy.integrate     import odeint

from common              import readDataEurope
from SolveDiff_SEIR1R2   import SolveDiff_SEIR1R2

# constante
fileLocalCopy = False         # if we upload the file from the url (to get latest results) or from a local copy file
startDate     = '2020-02-25'  # whatever the upload data starts, this sets the start date to be processed

def residual(paras, t, data, solveur):
    """
    compute the residual between actual data and fitted data
    """

    solveur.setParamInit(paras['E'].value, paras['I'].value, paras['R1'].value, paras['R2'].value)
    solveur.modele.setParam(paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
    model = solveur.solve_SEIR1R2_sol1(t)
    # print('model=', model)
    # input('apuse')

    # Only R1 to calculate the residual
    return (model[:, 3] - data).ravel()


if __name__ == '__main__':

    verbose = 1
    plot    = True
    country = 'France' #'United_Kingdom'

    # Modele d'eq. diff non lineaires
    #############################################################################
    # Solveur eq. diff.
    N       = 65.E6
    dt      = 1
    solveur = SolveDiff_SEIR1R2(N, dt, verbose)
    #solveur.setParamInit(1000, 10000, 5000*solveur.modele.f, 5000*(1.-solveur.modele.f))
    if verbose>0:
        print(solveur)

    prefixFig = './figures/' + solveur.modele.modelShortName + '_' + country

    # Lecture des données et copy of the observation
    #############################################################################

    pd, z_observ = readDataEurope(country=country, dateMin=startDate, dateMax=None, \
                            plot=plot, fileLocalCopy=fileLocalCopy, verbose=verbose)
    if pd.empty==True:
        print('pd is empty - No data! --> exit')
        exit(1)

    zs = []
    for z in pd[z_observ[0]]:
        zs.append(np.array([z]))


    # Solveur
    ############################################################################
    
    # Integration time grid
    simulLenght = 600 #len(zs) #600
    vectTime    = np.linspace(0, simulLenght-1, simulLenght)

    # with the inital values
    solveur.solve_SEIR1R2_sol1(vectTime)
    if plot==True:
        listePlot=[3]
        solveur.plot_SEIR1R2(prefixFig+'_Fitinit.png', vectTime, plot=listePlot, zs=(zs,270))
    input('attente')

    # set parameters including bounds; you can also fix parameters (use vary=False)
    S0, E0, I0, R10, R20 = solveur.getParamInit()
    a0, b0, c0, f0       = solveur.modele.getParam()
 
    params = Parameters()
    params.add('N',  value=N,   vary=False)
    params.add('S',  value=S0,  vary=False) #, min=1E7,    max=6.5E7)
    params.add('E',  value=E0,  vary=False) #, min=0.,     max=200)
    params.add('I',  value=I0,  vary=False) #, min=0.,     max=200)
    params.add('R1', value=R10, vary=False) #, min=0.,     max=200)
    params.add('R2', value=R20, vary=False) #, min=0.,     max=200)
    params.add('a',  value=a0,  vary=False) #, min=0.0001, max=0.9999)
    params.add('b',  value=b0,  vary=False) #, min=0.0001, max=0.9999)
    params.add('c',  value=c0,  vary=True,  min=0.0001, max=0.9999)
    params.add('f',  value=f0,  vary=True,  min=0.0001, max=0.9999)
    
    # fit model
    result = minimize(residual, params, args=(vectTime, zs, solveur), method='leastsq')  # leastsq nelder
    
    # with the optimized values
    solveur.setParamInit(result.params['E'], result.params['I'], result.params['R1'], result.params['R2'])
    solveur.modele.setParam(result.params['a'], result.params['c'], result.params['b'], result.params['f'])
    solveur.solve_SEIR1R2_sol1(vectTime)
    if plot==True:
        solveur.plot_SEIR1R2(prefixFig+'_Fit.png', vectTime)
    
    # # call main function
    # S_lmfit, E_lmfit, I_lmfit, R1_lmfit, R2_lmfit = solve_SEIR1R2_sol1(y0_lmfit, t, N, a_lmfit, b_lmfit, c_lmfit, f_lmfit)
    # # Plot the data on three separate curves for S(t), E(t), I(t), R1(t) and R2(t)
    # plot_SEIR1R2('prefixFig+fitsimSEIR1R2model_lmfit.png', timestep, N, S_lmfit, E_lmfit, I_lmfit, R1_lmfit, R2_lmfit, zs)
