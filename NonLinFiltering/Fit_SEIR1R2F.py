#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from datetime          import datetime, timedelta
  
from common            import getDates, addDaystoStrDate, get_WE_indice, drawAnnotation
from common            import getLowerDateFromString, getNbDaysBetweenDateFromString, getRepertoire
from ProcessSEIR1R2F   import fit 
from SolveEDO_SEIR1R2F import SolveEDO_SEIR1R2F

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

def main(argv):
    
    # Le decalalge ets traité comme argument pour faciliter les traitements systématiques
    if len(argv)>1:
        decalage3periodes = int(sys.argv[1])
    else:
        decalage3periodes = 10

    verbose    = 0
    UKF_filt   = False #False #True

    #places     = 'France,Spain,Italy,United_Kingdom,South_Korea'
    # places     = 'Germany'#,South_Korea'
    places     = 'France'
    # places     = 'France,69,75'
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
    decalage   = 2
    
    model_allinone, ListeTextParamPlace_allinone, liste_pd_allinone, data_deriv_allinone, modelR1_deriv_allinone, decalage_corrige, _, _ = \
            fit([places, nbperiodes, decalage, UKF_filt, verbose, 0])

    ListeTestPlace = []
    for indexplace in range(len(listplaces)):
        texteplace = ''
        for texte in ListeTextParamPlace_allinone[indexplace]:
            texteplace += '\n' + texte
            print(texte.replace('$', '').replace('\\', ''))
        ListeTestPlace.append(texteplace)

    PlotPlace(data_deriv_allinone, modelR1_deriv_allinone, listplaces, decalage_corrige, UKF_filt, FrDatabase, 'All-In-One', ListeTestPlace)
 
    # fit avec 3 périodes + décalage
    ##################################
    nbperiodes = -1

    model_piecewise, ListeTextParamPlace_piecewise, liste_pd_piecewise, data_deriv_piecewise, modelR1_deriv_piecewise, decalage_corrige, _, _ = \
            fit([places, nbperiodes, decalage3periodes, UKF_filt, verbose, 0])
    
    ListeTestPlace = []
    for indexplace in range(len(listplaces)):
        texteplace = ''
        for texte in ListeTextParamPlace_piecewise[indexplace]:
            texteplace += '\n' + texte
            print(texte.replace('$', '').replace('\\', ''))
        ListeTestPlace.append(texteplace)

    PlotPlace(data_deriv_piecewise, modelR1_deriv_piecewise, listplaces, decalage_corrige, UKF_filt, FrDatabase, 'Piecewise', ListeTestPlace)


    # Plot the two strategies in one figure
    ##################################
    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
        else:
            placefull   = place

        # Repertoire des figures
        repertoire = getRepertoire(UKF_filt, './figures/SEIR1R2F_UKFilt/'+placefull, './figures/SEIR1R2F/' + placefull)
        prefFig    = repertoire + '/' + placefull
        
        # Preparation plot pandas
        listheader = list(liste_pd_piecewise[indexplace])
        # on ajoute les dérivées numériques des cas et des morts
        liste_pd_piecewise[indexplace]['dcases']  = liste_pd_piecewise[indexplace][listheader[0]].diff()
        liste_pd_piecewise[indexplace]['ddeaths'] = liste_pd_piecewise[indexplace][listheader[1]].diff()
        longueur = len(liste_pd_piecewise[indexplace].loc[:, ('dcases')])

        if FrDatabase==True:
            DatesString = getDates('France', verbose)
        else:
            DatesString = getDates(place, verbose)

        # on ajoute les dérivées numériques des cas et des morts
        liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise')] = modelR1_deriv_piecewise[indexplace, 0:longueur, 0]
        liste_pd_piecewise[indexplace].loc[:, ('md_piecewise')] = modelR1_deriv_piecewise[indexplace, 0:longueur, 1]
        
        # filename = prefFig + '_Diff_BothFit_Shift' + str(decalage_corrige) + '.png'
        # title    = place + ' - Shift=' + str(decalage_corrige) + ' day(s)'
        # PlotAllPandas(liste_pd_piecewise[indexplace], title, filename, y=['dcases', 'mc_piecewise', 'ddeaths', 'md_piecewise'], Dates=DatesString, textannotation=ListeTestPlace[indexplace])

    # Plot the two SEIR1R2F
    ##################################

    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
        else:
            placefull   = place

        # Repertoire des figures
        repertoire = getRepertoire(UKF_filt, './figures/SEIR1R2F_UKFilt/'+placefull, './figures/SEIR1R2F/' + placefull)
        prefFig = repertoire + '/' + placefull

        # Preparation plot pandas
        listheader = list(liste_pd_piecewise[indexplace])
        longueur = len(liste_pd_piecewise[indexplace].loc[:, (listheader[0])])
        # print('longueur=', longueur)

        if FrDatabase==True:
            DatesString = getDates('France', verbose)
        else:
            DatesString = getDates(place, verbose)

        liste_pd_piecewise[indexplace].loc[:, ('Sp')]  = model_piecewise [indexplace, :, 0]
        liste_pd_piecewise[indexplace].loc[:, ('Ep')]  = model_piecewise [indexplace, :, 1]
        liste_pd_piecewise[indexplace].loc[:, ('Ip')]  = model_piecewise [indexplace, :, 2]
        liste_pd_piecewise[indexplace].loc[:, ('R1p')] = model_piecewise [indexplace, :, 3]
        liste_pd_piecewise[indexplace].loc[:, ('R2p')] = model_piecewise [indexplace, :, 4]
        liste_pd_piecewise[indexplace].loc[:, ('Fp')]  = model_piecewise [indexplace, :, 5]
        # print(liste_pd_piecewise[indexplace].head())
        # input('pause')

        title    = place + ' - Shift=' + str(decalage_corrige) + ' day(s)'
        listePlot=['Ep', 'Ip', 'R1p', 'Fp']
        filename = prefFig +'_' + ''.join(map(str, listePlot)) + '_piecewise_Shift' + str(decalage_corrige) + '.png'
        PlotSEIR1R2Pandas(liste_pd_piecewise[indexplace], title, filename, y=listePlot, Dates=DatesString, textannotation=ListeTestPlace[indexplace])

        title    = place + ' - Shift=' + str(decalage_corrige) + ' day(s)'
        listePlot=['R1p', 'Fp']
        filename = prefFig +'_' + ''.join(map(str, listePlot)) + '_piecewise_Shift' + str(decalage_corrige) + '.png'
        PlotSEIR1R2Pandas(liste_pd_piecewise[indexplace], title, filename, y=listePlot, Dates=DatesString, textannotation=ListeTestPlace[indexplace])


def PlotSEIR1R2Pandas(pd, titre, filenameFig, y, Dates=None, textannotation=''):
    
    if len(y)==0 or y is None: pass

    # juste pour récupérer des couleurs sohérentes
    solveur = SolveEDO_SEIR1R2F(N=1, dt=1, verbose=0)
    mod = solveur.modele
    listeColor = []
    for p in y:
        if p=='Sp' : listeColor.append(mod.getColor(0))
        if p=='Ep' : listeColor.append(mod.getColor(1))
        if p=='Ip' : listeColor.append(mod.getColor(2))
        if p=='R1p': listeColor.append(mod.getColor(3))
        if p=='R2p': listeColor.append(mod.getColor(4))
        if p=='Fp' : listeColor.append(mod.getColor(5))

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # Dessin des courbes théoriques
    pd.plot(ax=ax, y=y, title=titre, color=listeColor)
    
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

    # ajout d'un text d'annotation
    if textannotation != '':
        bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
        ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/3.), fontsize=6, bbox=bbox, ha="left", va="center") 

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
    plt.tight_layout(rect=(0, 0, 1., 0.95))
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(filenameFig, dpi=dpi)
    plt.close()


def PlotAllPandas(pd, titre, filenameFig, y, Dates=None, textannotation=''):

    if len(y)==0 or y is None: pass

    # on supprime la première ligne de donnée (car ce sont des dérivées)
    pd = pd.iloc[1:]

    # juste pour récupérer la bonne coukeur
    solveur = SolveEDO_SEIR1R2F(N=1, dt=1, verbose=0)
    mod     = solveur.modele
    colori  = [mod.getColor(3), mod.getColor(5)]

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # Dessin des courbes
    linestyles = ['-', '-', '-', '-']#, '-']
    markers    = ['x', '', '+', '']#, '']
    colors     = [colori[0], colori[0], colori[1], colori[1]]#, 'blue']
    linewidths = [0.5, 1.5, 0.5, 1.5]#, 1.5]
    alphas     = [1.0, 0.7, 1.0, 0.7]#, 0.55]
    labels     = [r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise', \
                    r'$\frac{\partial F(n)}{\partial n}$', r'$\frac{\partial F(t)}{\partial t}$ - Piecewise'] #, r'$\frac{\partial R^1(t)}{\partial t}$ - All in one']
    for col, ls, lw, l, a, c, m in zip(y, linestyles, linewidths, labels, alphas, colors, markers):
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

    # ajout d'un text d'annotation
    if textannotation != '':
        bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
        ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/3.), fontsize=6, bbox=bbox, ha="left", va="center") 

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
    plt.legend(fontsize=8)
    plt.tight_layout(rect=(0, 0, 1., 0.95))
    plt.savefig(filenameFig, dpi=dpi)
    plt.close()


def PlotPlace(data_deriv, model_derive, listplaces, decalage_corrige, UKF_filt, FrDatabase, ch, listetextannotation):

    # Plot
    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
        else:
            placefull   = place

        # Repertoire des figures
        repertoire = getRepertoire(UKF_filt, './figures/SEIR1R2F_UKFilt/'+placefull, './figures/SEIR1R2F/' + placefull)
        prefFig = repertoire + '/' + placefull

        filename = prefFig + '_Diff_' + ch + '_Shift' + str(decalage_corrige) + '.png'
        title    = place + ' - Shift=' + str(decalage_corrige) + ' day(s)'
        PlotFit(data_deriv[indexplace], model_derive[indexplace], title, filename, ch, listetextannotation[indexplace])


def PlotFit(data_deriv, model_deriv, title, filename, ch, textannotation=''):

    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # juste pour récupérer la bonne coukeur
    solveur = SolveEDO_SEIR1R2F(N=1, dt=1, verbose=0)
    mod     = solveur.modele
    colori  = [mod.getColor(3), mod.getColor(5)]

    time = np.linspace(1, len(data_deriv[1:]), len(data_deriv[1:]))
    ax.plot(time, data_deriv [1:, 0], label=r'$\frac{\partial R^1(n)}{\partial n}$ - ' + ch, color=colori[0], alpha=1.0, lw=0.5, marker='x')
    ax.plot(time, model_deriv[1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - ' + ch, color=colori[0], alpha=1.0, lw=1.5)
    ax.plot(time, data_deriv [1:, 1], label=r'$\frac{\partial F(n)}{\partial n}$ - '   + ch, color=colori[1], alpha=1.0, lw=0.5, marker='+')
    ax.plot(time, model_deriv[1:, 1], label=r'$\frac{\partial F(t)}{\partial t}$ - '   + ch, color=colori[1], alpha=1.0, lw=1.5)

    # ajout d'un text d'annotation
    if textannotation != '':
        bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
        ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=6, bbox=bbox, ha="left", va="center") 

    ax.set_xlabel('Time (days)')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=1, ls='-')
    
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    plt.tight_layout(rect=(0, 0, 1., 0.95))
    plt.title(title)
    plt.savefig(filename, dpi=dpi)
    plt.close()

if __name__ == '__main__':
    main(sys.argv)
