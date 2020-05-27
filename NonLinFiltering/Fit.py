#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt
from datetime       import datetime, timedelta

from common         import getDates, Plot, addDaystoStrDate, get_WE_indice, drawAnnotation
from common         import getLowerDateFromString, getNbDaysBetweenDateFromString
from ProcessSEIR1R2 import fit

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

def main():
    verbose       = 0
    plot          = 0

    # places     = 'France,Spain,Italy,United_Kingdom,Germany,Belgium'
    #places     = 'France,Spain,Italy,Germany'
    places      = 'France'
    #places     = 'United_Kingdom'
    #places     = 'France,69,75,01'
    #places     = 'France,93'
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

    # fit avec une seule période
    ##################################
    nbperiodes = 1
    decalage   = 0
    surplus    = 0
    modelSEIR1R2_allinone, ListeTextParamPlace_allinone, _, data_deriv_allinone, moldelR1_deriv_allinone = \
        fitplace(places, listplaces, nbperiodes, decalage, surplus, verbose, plot, 'All-In-One')

    for indexplaces in range(len(listplaces)):
        for texte in ListeTextParamPlace_allinone[indexplaces]:
            print(texte)
    

    # fit avec 3 périodes + décalage
    ##################################
    nbperiodes = -1
    decalage   = 5
    surplus    = 0
    modelSEIR1R2_piecewise, ListeTextParamPlace_piecewise, listepd, data_deriv_piecewise, moldelR1_deriv_piecewise = \
        fitplace(places, listplaces, nbperiodes, decalage, surplus, verbose, plot, 'Piecewise')

    for indexplaces in range(len(listplaces)):
        for texte in ListeTextParamPlace_piecewise[indexplaces]:
            print(texte)

    # Plot the two stratégies in one figure
    ##################################
    for indexplaces in range(len(listplaces)):
    
        place    = listplaces[indexplaces]

        # schéma 1, moins joli que le suivant
        # filename = './figures/DiffR1_BothFitbis_' + place + '_Shift' + str(decalage) + '.png'
        # title    = place + ' - Shift=' + str(decalage) + ' day(s)'
        # PlotAll(data_deriv_allinone[indexplaces], moldelR1_deriv_allinone[indexplaces], moldelR1_deriv_piecewise[indexplaces], title, filename)

        # Preparation plot pandas
        listheader = list(listepd[indexplaces])
        # on ajoute les dérivées numériques des cas et des morts
        listepd[indexplaces]['dcases']  = listepd[indexplaces][listheader[0]].diff()
        listepd[indexplaces]['ddeaths'] = listepd[indexplaces][listheader[1]].diff()
        longueur = len(listepd[indexplaces].loc[:, ('dcases')])

        if FrDatabase==True:
            DatesString = getDates('France', verbose)
        else:
            DatesString = getDates(place, verbose)

        # on ajoute les 4 dérivées numériques des cas et des morts
        # listepd[indexplaces].loc[:, ('dd_allinone')] = data_deriv_allinone     [indexplaces, :]
        listepd[indexplaces].loc[:, ('md_allinone')] = moldelR1_deriv_allinone [indexplaces, 0:longueur]
        #listepd[indexplaces].loc[:, ('dd_piecewise')]= data_deriv_piecewise    [indexplaces, :]
        listepd[indexplaces].loc[:, ('md_piecewise')]= moldelR1_deriv_piecewise[indexplaces, 0:longueur]
        #print('tail=', listepd[indexplaces].tail())

        filename = './figures/DiffR1_BothFit_' + place + '_Shift' + str(decalage) + '.png'
        title    = place + ' - Shift=' + str(decalage) + ' day(s)'
        PlotAllPandas(listepd[indexplaces], title, filename, y=['dcases', 'md_allinone', 'md_piecewise'], Dates=DatesString)

    # Plot the two SIER1R2
    ##################################

    for indexplaces in range(len(listplaces)):
    
        place    = listplaces[indexplaces]

        # Preparation plot pandas
        listheader = list(listepd[indexplaces])
        # on ajoute les dérivées numériques des cas et des morts
        longueur = len(listepd[indexplaces].loc[:, (listheader[0])])
        # print('longueur=', longueur)

        if FrDatabase==True:
            DatesString = getDates('France', verbose)
        else:
            DatesString = getDates(place, verbose)

        # on ajoute les 4 dérivées numériques des cas et des morts
        # listepd[indexplaces].loc[:, ('dd_allinone')] = data_deriv_allinone     [indexplaces, :]
        #print(len(modelSEIR1R2_piecewise [indexplaces, :, 0]))
        listepd[indexplaces].loc[:, ('Sp')]  = modelSEIR1R2_piecewise [indexplaces, :, 0]
        listepd[indexplaces].loc[:, ('Ep')]  = modelSEIR1R2_piecewise [indexplaces, :, 1]
        listepd[indexplaces].loc[:, ('Ip')]  = modelSEIR1R2_piecewise [indexplaces, :, 2]
        listepd[indexplaces].loc[:, ('R1p')] = modelSEIR1R2_piecewise [indexplaces, :, 3]
        listepd[indexplaces].loc[:, ('R2p')] = modelSEIR1R2_piecewise [indexplaces, :, 4]
        # print(listepd[indexplaces].head())
        # input('pause')

        title    = place + ' - Shift=' + str(decalage) + ' day(s)'
        filename = './figures/EIR1_piecewise_' + place + '_Shift' + str(decalage) + '.png'
        y=['Ep', 'Ip', 'R1p']
        PlotSEIR1R2Pandas(listepd[indexplaces], title, filename, listeString=y, Dates=DatesString)

        title    = place + ' - Shift=' + str(decalage) + ' day(s)'
        filename = './figures/R1_piecewise_' + place + '_Shift' + str(decalage) + '.png'
        y=['R1p']
        PlotSEIR1R2Pandas(listepd[indexplaces], title, filename, listeString=y, Dates=DatesString)


def PlotSEIR1R2Pandas(pd, titre, filenameFig, listeString, Dates=None):
    
    if len(listeString)==0 or listeString is None: pass

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # Dessin des courbes théoriques
    pd.plot(ax=ax, y=listeString, title=titre)
    
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


def PlotAllPandas(pd, titre, filenameFig, y, Dates=None):

    if len(y)==0 or y is None: pass

    # on supprime la première ligne de donnée (car ce sont des dérivées)
    pd = pd.iloc[1:]

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # Dessin des courbes
    styles     = ['gx-','b-','g-']
    linewidths = [0.5, 1.5, 1.5]
    alphas     = [1.0, 0.55, 0.7]
    labels     = [r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial R^1(t)}{\partial t}$ - All in one', r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise']
    for col, style, lw, label, alpha in zip(y, styles, linewidths, labels, alphas):
        pd[col].plot(title=titre, style=style, lw=lw, ax=ax, alpha=alpha, label=label)
 
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
    legend = ax.legend().get_frame().set_alpha(0.5)
    plt.legend(fontsize=9)
    #plt.tight_layout()
    plt.savefig(filenameFig, dpi=dpi)
    plt.close()


def PlotAll(data_deriv_allinone, moldelR1_deriv_allinone, moldelR1_deriv_piecewise, title, filename):

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    time = np.linspace(1, len(data_deriv_allinone[1:]), len(data_deriv_allinone[1:]))
    ax.plot(time, data_deriv_allinone     [1:], label=r'$\frac{\partial R^1(n)}{\partial n}$',              color='green', alpha=1.0,  lw=0.5, marker='x')
    ax.plot(time, moldelR1_deriv_allinone [1:], label=r'$\frac{\partial R^1(t)}{\partial t}$ - All-in-one', color='blue',  alpha=0.55, lw=1.5)
    ax.plot(time, moldelR1_deriv_piecewise[1:], label=r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise',  color='green', alpha=0.7,  lw=1.5)

    ax.set_xlabel('Time (days)')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=1, ls='-')
    
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    plt.title(title)
    plt.savefig(filename, dpi=dpi)
    plt.close()


def PlotFit(data_deriv, model_deriv, title, ch, filename):

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    time = np.linspace(1, len(data_deriv[1:]), len(data_deriv[1:]))
    ax.plot(time, data_deriv[1:],  label=r'$\frac{\partial R^1(n)}{\partial n}$ - ' + ch, color='green', alpha=1.0, lw=0.5, marker='x')
    ax.plot(time, model_deriv[1:], label=r'$\frac{\partial R^1(t)}{\partial t}$ - ' + ch, color='green', alpha=1.0, lw=1.5)

    ax.set_xlabel('Time (days)')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=1, ls='-')
    
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    plt.title(title)
    plt.savefig(filename, dpi=dpi)
    plt.close()


def fitplace(places, listplaces, nbperiodes, decalage, surplus, verbose, plot, ch):

    modelSEIR1R2, ListeTextParamPlace, listepd, data_deriv, moldelR1_deriv = fit([places, nbperiodes, decalage, surplus, verbose, plot])
    
    # Plot
    for indexplaces in range(len(listplaces)):
        place    = listplaces[indexplaces]
        filename = './figures/DiffR1_' + place + '_' + ch + '_Shift' + str(decalage) + '.png'
        title    = place + ' - Shift=' + str(decalage) + ' day(s)'
        #dataLength = listepd[indexplaces].shape[0]
        PlotFit(data_deriv[indexplaces], moldelR1_deriv[indexplaces], title, ch, filename)

    return modelSEIR1R2, ListeTextParamPlace, listepd, data_deriv, moldelR1_deriv


if __name__ == '__main__':
    main()


