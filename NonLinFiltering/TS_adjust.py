#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
from datetime          import datetime, timedelta
  
from common           import getDates, Plot, addDaystoStrDate, get_WE_indice, drawAnnotation
from common           import getLowerDateFromString, getNbDaysBetweenDateFromString
from ProcessSEIR1R2   import fit
from SolveEDO_SEIR1R2 import SolveEDO_SEIR1R2

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)


def main():
    verbose       = 0
    plot          = 0

    #places     = 'France,Spain,Italy,United_Kingdom,Germany,Belgium'
    #places     = 'France,Spain,Italy,Germany'
    #places     = 'France,69,75,01'
    places     = 'Italy'#,South_Korea'
    listplaces = list(places.split(','))

    # Detect si pays européen ou département français
    FrDatabase = False
    if listplaces[0]=='France' and len(listplaces)>1:
        try:
            int(listplaces[1])
        except Exception as e:
            FrDatabase=False
        else:
            FrDatabase = True
            listplaces = listplaces[1:]
    

    # fit avec 3 périodes + décalage
    ##################################
    surplus    = 0
    nbperiodes = -1

    TAB_decalage_corrige = []
    TAB_param_model      = []
    TAB_ListeEQM         = []
    
    decalagemini, decalagemaxi = 4,14 #1, 11
    for decalage in range(decalagemini, decalagemaxi):

        print('TIME-SHIFT', str(decalage), 'OVER', str(decalagemaxi))
    
        _, _, _, _, _, decalage_corrige, tabParamModel, ListeEQM = fit([places, nbperiodes, decalage, surplus, verbose, plot])

        TAB_decalage_corrige.append(float(decalage_corrige))
        TAB_param_model.append(tabParamModel)
        TAB_ListeEQM.append(ListeEQM)

    #TAB_param_model[decalage][place][nbperiodes]

    
    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
        else:
            placefull   = place

        # Constantes
        import os
        repertoire = './figures/'+ placefull
        if not os.path.exists(repertoire):
            os.makedirs(repertoire)
        prefFig = repertoire + '/' + placefull

        nbperiodes = len(TAB_param_model[0][indexplace][:])
        X = np.linspace(decalagemini, decalagemaxi-1, decalagemaxi-decalagemini)
        labelsparam  = [r'$a$', r'$b$', r'$c$', r'$f$', r'$R_O$']
        labelsperiod = ['Period 1', 'Period 2', 'Period 3']

        # plot pour les 3 périodes
        ##########################################@
        Y1 = np.zeros(shape=(decalagemaxi-decalagemini, len(labelsparam)))
        for period in range(nbperiodes):

                try:
                    Y1[decalage, :] = TAB_param_model[decalage][indexplace][period][:]
                except IndexError:
                    Y1[decalage, :] = 0.
                

            titre    = placefull + ' - SEIR1R2 parameters evolution for ' + labelsperiod[period]
            filename = prefFig   + '_ParamEvol_Period' + str(period) + '_abcd.png'
            plotData(TAB_decalage_corrige, Y1[:, :-1], titre, filename, labelsparam[:-1])
            titre    = placefull + ' - SEIR1R2 parameters evolution for ' + labelsperiod[period]
            filename = prefFig   + '_ParamEvol_Period' + str(period) + '_R0.png'
            plotData(TAB_decalage_corrige, Y1[:, -1].reshape(decalagemaxi-decalagemini, 1), titre, filename, [labelsparam[-1]])

        # plot pour les paramètres
        ##########################################@
        Y2 = np.zeros(shape=(decalagemaxi-decalagemini, len(labelsperiod)))
        for param in range(len(labelsparam)):

            for decalage in range(decalagemaxi-decalagemini):
                for period in range(len(labelsperiod)):
                    try:
                        Y2[decalage, period] = TAB_param_model[decalage][indexplace][period][param]
                    except IndexError:
                        Y2[decalage, period] = 0.
            titre    = placefull + ' - SEIR1R2 periods evolution for param ' + labelsparam[param]
            filename = prefFig   + '_PeriodEvol_Param' + labelsparam[param] + '.png'
            plotData(TAB_decalage_corrige, Y2, titre, filename, labelsperiod)

        # plot de l'EQM
        ##########################################
        fig = plt.figure(facecolor='w', figsize=figsize)
        ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
        
        Y3 = np.zeros(shape=(len(X)))
        for k in range(len(Y3)):
            try:
                Y3[k] = TAB_ListeEQM[k][indexplace]
            except IndexError:
                Y3[k] = 0.
        ax.plot(TAB_decalage_corrige, Y3, alpha=1.0, lw=2, label='EQM')

        ax.set_xlabel('Time shift (days)')
        ax.yaxis.set_tick_params(length=0)
        ax.xaxis.set_tick_params(length=0)
        ax.grid(b=True, which='major', c='w', lw=1, ls='-')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        legend = ax.legend()
        legend.get_frame().set_alpha(0.5)
        for spine in ('top', 'right', 'bottom', 'left'):
            ax.spines[spine].set_visible(False)

        plt.xlim([TAB_decalage_corrige[0], TAB_decalage_corrige[-1]])
        #plt.ylim([0, 1.0])

        # ajout d'un text d'annotation
        plt.title(placefull + ' - SEIR1R2, EQM on the number of detected cases' )
        plt.savefig(prefFig + '_EQM_Deriv.png', dpi=dpi)
        plt.close()


def plotData(X, Y, titre, filename, labels):
    
    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    
    for i in range(len(labels)):
        ax.plot(X, Y[:, i], alpha=1.0, lw=2, label=labels[i])

    ax.set_xlabel('Time shift (days)')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=1, ls='-')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    plt.xlim([X[0], X[-1]])
    if np.max(np.max(Y)) < 1.:
        plt.ylim([0, 1.])

    # ajout d'un text d'annotation
    plt.title(titre)
    plt.savefig(filename, dpi=dpi)
    plt.close()


if __name__ == '__main__':
    main()


