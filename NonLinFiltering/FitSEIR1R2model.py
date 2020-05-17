#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
from lmfit               import minimize, Parameters, Parameter, report_fit
from scipy.integrate     import odeint
from scipy import signal
from sklearn.metrics import mean_squared_error


from common              import readDataEurope
from SolveDiff_SEIR1R2   import SolveDiff_SEIR1R2

# constante
fileLocalCopy = True         # if we upload the file from the url (to get latest results) or from a local copy file
startDate     = '2020-02-25'  # whatever the upload data starts, this sets the start date to be processed



def residual(paras, t, data, solveur):
    """
    compute the residual between actual data and fitted data
    """

    # print('current value', paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value, ' -- ', paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
    # print('current ts', int(paras['ts'].value))
    
    solveur.setParamInit(paras['E0'].value, paras['I0'].value, paras['R10'].value, paras['R20'].value)
    solveur.modele.setParam(paras['N'].value, paras['a'].value, paras['b'].value, paras['c'].value, paras['f'].value)
    model = solveur.solve_SEIR1R2_sol1(t)
    # print('model=', model)
    # input('apuse')

    # Only R1 to calculate the residual
    eqm=[]
    # print(paras['ts'].min, paras['ts'].max)
    for t in range(paras['ts'].min, paras['ts'].max):
        eqm.append(mean_squared_error(model[t:t+len(data), 3], data))
    delay = np.argmin(eqm)
    ts    = paras['ts'].min+delay
    paras['ts'].set(ts)
    # print('delay=', delay)
    # print('next ts', int(paras['ts'].value))

    # Only R1 to calculate the residual, only on the window's size of the data
    result = (model[ts:ts+len(data), 3] - data).ravel()
    # paras.pretty_print()
    print('sum=', np.sum(np.abs(result)))

    return result


if __name__ == '__main__':

    verbose = 1
    plot    = True
    country = 'France' # 'United_Kingdom'

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
    simulLenght = 550 #len(zs) #600
    vectTime    = np.linspace(0, simulLenght-1, simulLenght)


    # with the inital values
    solut = solveur.solve_SEIR1R2_sol1(vectTime)
    print('solveur=', solveur)

    # calcul timeshift initial
    eqm=[]
    for t in range(0, 550-len(zs)):
        eqm.append(mean_squared_error(solut[t:t+len(zs), 3], zs))
    ts0 = np.argmin(eqm)

    # plot
    if plot==True:
        listePlot=[3]
        print('new ts=', ts0)
        solveur.plot_SEIR1R2(prefixFig+'_Fitinit.png', vectTime, plot=listePlot, zs=(zs, ts0))

    # Parameters optimization
    ############################################################################

    # set parameters including bounds; you can also fix parameters (use vary=False)
    S0, E0, I0, R10, R20 = solveur.getParamInit()
    _, a0, b0, c0, f0    = solveur.modele.getParam()
 
    params = Parameters()
    params.add('N',   value=N,   vary=False)
    params.add('S0',  value=S0,  vary=False) #, min=1E7,    max=6.5E7)
    params.add('E0',  value=E0,  vary=False) #, min=0.,     max=200)
    params.add('I0',  value=I0,  vary=False) #, min=0.,     max=200)
    params.add('R10', value=R10, vary=False) #, min=0.,     max=200)
    params.add('R20', value=R20, vary=False) #, min=0.,     max=200)
    params.add('a',   value=a0,  vary=True,  min=0.001, max=0.999)  # min=0.7*a0,        max=1.3*a0)
    params.add('b',   value=b0,  vary=True,  min=0.001, max=0.999)  # min=0.7*b0,        max=2.1*b0)
    params.add('c',   value=c0,  vary=True,  min=0.001, max=0.999)  #=0.9*c0,        max=1.8*c0)
    params.add('f',   value=f0,  vary=False, min=0.001, max=0.005)
    params.add('ts',  value=ts0, vary=True,  min=ts0-230, max=ts0+30)
    
    # fit model
    result = minimize(residual, params, args=(vectTime, zs, solveur), method='dual_annealing')  # leastsq nelder
    if verbose>0:
        result.params.pretty_print()
    
    # with the optimized values
    solveur1 = SolveDiff_SEIR1R2(N, dt, verbose)
    solveur1.setParamInit(result.params['E0'].value, result.params['I0'].value, result.params['R10'].value, result.params['R20'].value)
    solveur1.modele.setParam(result.params['N'].value, result.params['a'].value, result.params['b'].value, result.params['c'].value, result.params['f'].value)
    print('solveur1=', solveur1)
    solveur1.solve_SEIR1R2_sol1(vectTime)
    if plot==True:
        listePlot=[3]
        ts = result.params['ts'].value
        print('new ts=', ts)
        solveur1.plot_SEIR1R2(prefixFig+'_Fit.png', vectTime, plot=listePlot, zs=(zs, ts))

    
    # # call main function
    # S_lmfit, E_lmfit, I_lmfit, R1_lmfit, R2_lmfit = solve_SEIR1R2_sol1(y0_lmfit, t, N, a_lmfit, b_lmfit, c_lmfit, f_lmfit)
    # # Plot the data on three separate curves for S(t), E(t), I(t), R1(t) and R2(t)
    # plot_SEIR1R2('prefixFig+fitsimSEIR1R2model_lmfit.png', timestep, N, S_lmfit, E_lmfit, I_lmfit, R1_lmfit, R2_lmfit, zs)
