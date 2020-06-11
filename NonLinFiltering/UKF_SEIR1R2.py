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

from common            import readDataEurope, readDataFrance, getDates, addDaystoStrDate
from common            import getLowerDateFromString, getNbDaysBetweenDateFromString, GetPairListDates
from common            import drawAnnotation, get_WE_indice
from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)
strDate = "%Y-%m-%d"

def main(argv):
    """
        Program to perform UKF filtering on Covid Data.
 
        :Example:

        >> python3 UKF_SEIR.py 
        >> python3 UKF_SEIR.py Italy
        >> python3 UKF_SEIR.py Italy 1 1
        >> python3 UKF_SEIR.py France,Germany 1 1 # Shortcut for processing the two countries successively

        argv[1] : Country (or list separeted by ','), or Continent, or 'World'. Default: France 
        argv[2] : Verbose level (debug: 3, ..., almost mute: 0). Default: 1
        argv[3] : Plot graphique (0/1).                          Default: 1
    """

    # Interpetation of arguments - reparation
    ######################################################@

    print('Command line : ', argv, flush=True)
    if len(argv) > 4:
        print('  CAUTION : bad number of arguments - see help')
        exit(1)

    # Default value for parameters
    listplaces = ['France']
    nbperiodes = 1
    decalage   = 0
    verbose    = 1
    plot       = True

    # Parameters from argv
    if len(argv)>1: listplaces = list(argv[1].split(','))
    FrDatabase = False
    if listplaces[0]=='France' and len(listplaces)>1:
        try:
            int(listplaces[1]) # If this is a number, then it is a french dpt
        except Exception as e:
            FrDatabase=False
        else:
            FrDatabase = True
            listplaces = listplaces[1:]

    if len(argv)>2: verbose = int(argv[2])
    if len(argv)>3 and int(argv[3])==0: plot = False

    # Constantes
    ######################################################@
    fileLocalCopy    = True  # if we upload the file from the url (to get latest data) or from a local copy file
    readStartDateStr = "2020-03-01" #"2020-03-01" 
    readStopDateStr  = None
    recouvrement     = 0
    dt               = 1
    indexdata        = [3]

    # Data reading to get first and last date available in the data set
    ######################################################
    if FrDatabase == True: 
        pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance('69', readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
    else:
        pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope('France', readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
    dataLength = pd_exerpt.shape[0]

    readStartDate = datetime.strptime(readStartDateStr, strDate)
    if readStartDate<pd_exerpt.index[0]:
        readStartDate    = pd_exerpt.index[0]
        readStartDateStr = pd_exerpt.index[0].strftime(strDate)
    readStopDate = datetime.strptime(readStopDateStr, strDate)
    if readStopDate<pd_exerpt.index[-1]:
        readStopDate     = pd_exerpt.index[-1]
        readStopDateStr  = pd_exerpt.index[-1].strftime(strDate)

    dataLength = pd_exerpt.shape[0]
    if verbose>1:
        print('readStartDateStr=', readStartDateStr, ', readStopDateStr=', readStopDateStr)
        print('readStartDate   =', readStartDate,    ', readStopDate   =', readStopDate)
        #input('pause')

    # Surplus and correction from previou period
    data = np.zeros(shape=(dataLength, len(indexdata)))

    # Loop on the places to process
    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
            DatesString = getDates('France', verbose)
        else:
            placefull   = place
            DatesString = getDates(place, verbose)

        print('PROCESSING of', placefull, 'in', listplaces)

        # data reading of the observations
        #############################################################################
        if FrDatabase == True:
            pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataFrance(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
        else:
            pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(place, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)

        # Get the list of dates to process
        ListDates, ListDatesStr = GetPairListDates(readStartDate, readStopDate, DatesString, decalage, nbperiodes, recouvrement)
        if verbose>1:
            print('ListDates   =', ListDates)
            print('ListDatesStr=', ListDatesStr)
            print('HeadData=', HeadData)
    
        # Solveur edo
        solveur = SolveEDO_SEIR1R2(N, dt, verbose)
        E0, I0, R10, R20 = 0, 1, 0, 0

        # Constantes
        import os
        repertoire = './figures/SEIR1R2/' + placefull
        if not os.path.exists(repertoire):
            os.makedirs(repertoire)
        prefFig = repertoire + '/' + placefull + '_UKF'
        
        # Remise à 0 de la correction
        data.fill(0.)

        # Boucle pour traiter successivement les différentes fenêtres
        ###############################################################

        for i in range(len(ListDatesStr)):

            # dates of the current period
            fitStartDate,    fitStopDate    = ListDates[i]
            fitStartDateStr, fitStopDateStr = ListDatesStr[i]

            if i>0:
                DatesString.addOtherDates(fitStartDateStr)

            # Récupérations des données observées
            dataLengthPeriod = 0
            indMinPeriod     = (fitStartDate-readStartDate).days

            for j, z in enumerate(pd_exerpt.loc[fitStartDateStr:addDaystoStrDate(fitStopDateStr, -1), (HeadData[0])]):
                data[indMinPeriod+j, 0] = z
                dataLengthPeriod +=1
            slicedata      = slice(indMinPeriod, indMinPeriod+dataLengthPeriod)
            slicedataderiv = slice(slicedata.start+1, slicedata.stop)
            sliceedo       = slice(0, dataLengthPeriod)
            sliceedoderiv  = slice(sliceedo.start+1, sliceedo.stop)
            if verbose>0:
                print('  dataLength      =', dataLength)
                print('  indMinPeriod    =', indMinPeriod)
                print('  dataLengthPeriod=', dataLengthPeriod)
                print('  fitStartDateStr =', fitStartDateStr)
                print('  fitStopDateStr  =', fitStopDateStr)
                #input('attente')

            # Set initialisation values for the solver
            ############################################################################
            # pour une période
            a0, b0, c0, f0  = 0.10, 0.29, 0.10, 0.0022
            
            # MAJ des parametres
            solveur.modele.setParam(N=N,  a=a0,  b=b0,   c=c0,    f=f0)
            solveur.setParamInit   (N=N, E0=E0, I0=I0, R10=R10, R20=R20)

            # Unscented KF
            #############################################################################
            sigmas = MerweScaledSigmaPoints(n=solveur.modele.dimState, alpha=.5, beta=2., kappa=solveur.modele.dimState-3.)
            ukf    = UKF(dim_x=solveur.modele.dimState, dim_z=solveur.modele.dimObs, fx=solveur.modele.fx, hx=solveur.modele.hx, dt=solveur.modele.dt, points=sigmas)
            # Filter init
            ukf.x[1] = E0
            ukf.x[2] = I0
            ukf.x[3] = R10
            ukf.x[4] = R20
            ukf.x[0] = solveur.modele.N-np.sum(ukf.x[1:])
            # Filter noise, starting with the measurment noise (R), and then the process noise (Q)
            ukf.Q[0:1, 0:1] = 50.*50.   #Q_discrete_white_noise(1, dt=dt, var=1000)
            ukf.Q[1:2, 1:2] = 15.*15.   #Q_discrete_white_noise(1, dt=dt, var=100)
            ukf.Q[2:3, 2:3] = 20.*20.   #Q_discrete_white_noise(1, dt=dt, var=50)
            ukf.Q[3:4, 3:4] = 10.*10.   #Q_discrete_white_noise(1, dt=dt, var=10)
            ukf.Q[4:5, 4:5] = 10.*10.   #Q_discrete_white_noise(1, dt=dt, var=10)
            # ukf.Q           = np.diag([1., 1., 1., 1., 1.])
            ukf.R           = np.diag([15.*15.])
            if verbose>1:
                print('ukf.x[1:]=', ukf.x[1:])
                print('ukf.R=', ukf.R)
                print('ukf.Q=', ukf.Q)
                # input('pause')

            # print('sigmas=', sigmas)
            # input('apuse')
            # print('ukf=', ukf)
            # input('apuse')


            # UKF filtering, day by day
            # xMean_est = []
            # for z in data:
            #   ukf.predict()
            #   ukf.update(z)
            #   xMean_est.append(ukf.x.copy())
            # xMean_est = np.array(xMean_est)

            # UKF filtering and smoothing, batch mode
            xMean_est, _ = ukf.batch_filter(data)
            # UKF filtering and smoothing, batch mode -- DOESN'T WORK!!!
            # uxs, ucovs = ukf.batch_filter(data)
            # xMean_est, _, _ = ukf.rts_smoother(uxs, ucovs)
            
            # Panda dataframe
            columns = [r'$S(t)$', r'$E(t)$', r'$I(t)$', r'$R^1(t)$', r'$R^2(t)$', r'$R^2(t)=N-\sum(SEIR^1)$', r'$R(t)=R^1(t)+R^2(t)$']
            for j in range(5):
                pd_exerpt[columns[j]]=xMean_est[:, j]
            pd_exerpt[columns[5]] = solveur.modele.N-np.sum(xMean_est[:, 0:4], axis=1)
            pd_exerpt[columns[6]] = np.sum(xMean_est[:, 3:5], axis=1)
            if verbose>1:
                print(pd_exerpt.tail())
            
            # Les dérivées numériques
            data_deriv_period    = (data[slicedataderiv, :] - data[slicedataderiv.start-1:slicedataderiv.stop-1, :])        / dt
            modelR1_deriv_period = (xMean_est[sliceedoderiv, indexdata] - xMean_est[sliceedoderiv.start-1:sliceedoderiv.stop-1, indexdata]) / dt

            
            # save the generated data in csv form
            pd_exerpt.to_csv(prefFig + '_all.csv', header=True, index=False)

            if plot == True:

                titre = placefull

                # listePlot =[0,1,2,3,4]
                # filename  = prefFig + '_' + ''.join(map(str, listePlot)) + '.png'
                # Plot(pd_exerpt, titre, filename, solveur.modele, y=listePlot, Dates=DatesString)
                # listePlot =[1,2,3]
                # filename  = prefFig + '_' + ''.join(map(str, listePlot)) + '.png'
                # Plot(pd_exerpt, titre, filename, solveur.modele, y=listePlot, Dates=DatesString, data='cases')
                listePlot =[3]
                filename  = prefFig + '_' + ''.join(map(str, listePlot)) + '.png'
                Plot(pd_exerpt, titre, filename, solveur.modele, y=listePlot, Dates=DatesString, data='cases')

                listePlot =[3]
                filename  = prefFig + '_' + ''.join(map(str, listePlot)) + '_bis.png'
                solveur.plotEDO(filename, titre, sliceedo, slicedata, plot=listePlot, data=data, text=solveur.getTextParam(readStartDateStr))

                #plotEDO(self, name, title, sliceedo, slicedata, plot, data='', text=''):
                

def Plot(pd, titre, filenameFig, modele, y, Dates=None, data=''):

    if len(y)==0 or y is None: pass

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    listeString, listeColor, lineStyles, listeMarkers, listeLW, listeAlphas = [], [], [], [], [], []
    for p in y:
        listeString.append(modele.getString(p))
        listeColor.append(modele.getColor(p))
        lineStyles.append('-')
        listeMarkers.append('')
        listeLW.append(1.5)
        listeAlphas.append(0.7)
    # pour les données numériques
    if data != '':
        listeString.append(data)
        listeColor.append(modele.getColor(3))
        lineStyles.append('-')
        listeMarkers.append('x')
        listeLW.append(0.5)
        listeAlphas.append(1.)


    for col, ls, lw, l, a, c, m in zip(listeString, lineStyles, listeLW, listeString, listeAlphas, listeColor, listeMarkers):
        pd[col].plot(title=titre, ax=ax, ls=ls, lw=lw, label=l, alpha=a, color=c, marker=m)

    # ajout des dates spéciales
    if Dates!=None:
        for d in Dates.listConfDates:
            drawAnnotation(ax, 'Conf. date\n', d, color='red')
        for d in Dates.listDeconfDates:
            drawAnnotation(ax, 'Deconf. date\n', d, color='green')
        for d in Dates.listOtherDates:
            drawAnnotation(ax, 'Other date\n', d)

    # surlignage des jours de WE
    WE_indices = get_WE_indice(pd)
    i = 0
    while i < len(WE_indices)-1:
        ax.axvspan(pd.index[WE_indices[i]], pd.index[WE_indices[i+1]], facecolor='gray', edgecolor='none', alpha=.15, zorder=-100)
        i += 2

    # axes
    ax.grid(True, which='major', axis='both')
    ax.grid(True, which='minor', axis='both')
    ax.grid(True, which='major', c='k', lw=0.5, ls='-', alpha=0.3)
    ax.grid(True, which='minor', c='w', lw=0.5, ls='-')
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False, useLocale=False)

    # On enlève le label sur l'axe x
    x_label = ax.axes.get_xaxis().get_label().set_visible(False)

    # legende
    legend = ax.legend().get_frame().set_alpha(0.8)
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(filenameFig, dpi=dpi)
    plt.close()

if __name__ == '__main__':
    main(sys.argv)
