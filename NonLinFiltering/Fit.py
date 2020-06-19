#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt

from datetime          import datetime, timedelta
  
from common            import getDates, addDaystoStrDate, get_WE_indice, drawAnnotation
from common            import getLowerDateFromString, getNbDaysBetweenDateFromString, getRepertoire
from SolveEDO_SEIR1R2  import SolveEDO_SEIR1R2
from SolveEDO_SEIR1R2D import SolveEDO_SEIR1R2D
from ProcessSEIR1R2D   import fit as fitProcessSEIR1R2D
from ProcessSEIR1R2    import fit as fitProcessSEIR1R2

dpi     = 120    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

def main(sysargv):

    """
        :Example:

        For countries (European database)
        >> python3 Fit.py 
        >> python3 Fit.py France SEIR1R2 13 0 1 1
        >> python3 Fit.py Italy,Spain SEIR1R2D 13 1 1 1 # Italy and Spain, with UKF filtering

        For French Region (French database)
        >> python3 Fit.py France,69    SEIR1R2  13 0 1 1 # Dpt 69 (Rhône)
        >> python3 Fit.py France,69,01 SEIR1R2D 13 1 1 1 # Dpt 69 (Rhône) + Dpt 01 (Ain) with UKF filtering
        
        argv[1] : Country (or list separeted by ','), or 'France' followed by a list of dpts. Default: France 
        argv[2] : EDO model (SEIR1R2 or SEIR1R2D)                                             Default: SEIR2R2         
        argv[3] : Shift of periods (in days).                                                 Default: 0
        argv[4] : UKF filtering of data (0/1).                                                Default: 0
        argv[5] : Verbose level (debug: 3, ..., almost mute: 0).                              Default: 1
        argv[6] : Plot graphique (0/1).                                                       Default: 1
    """

    # Interpetation of arguments - reparation
    ######################################################@

    print('Command line : ', sysargv, flush=True)
    if len(sysargv) > 7:
        print('  CAUTION : bad number of arguments - see help')
        exit(1)

    # Default value for parameters
    places       = 'France'
    listplaces   = list(places.split(','))
    modeleString = 'SEIR1R2'
    decalage3P   = 10
    UKF_filt, UKF_filt01 = False, 0
    verbose      = 1
    plot         = True

    # Parameters from argv
    if len(sysargv)>1: places, listplaces = sysargv[1], list(sysargv[1].split(','))
    FrDatabase = False
    if listplaces[0]=='France' and len(listplaces)>1:
        try:
            int(listplaces[1]) #If this is a number, then it is a french dpt
        except Exception as e:
            FrDatabase=False
        else:
            FrDatabase = True
            listplaces = listplaces[1:]

    if len(sysargv)>2: modeleString  = sysargv[2]
    if len(sysargv)>3: decalage3P    = int(sysargv[3])
    if len(sysargv)>4 and int(sysargv[4])==1: UKF_filt, UKF_filt01 = True, 1
    if len(sysargv)>5: verbose       = int(sysargv[5])
    if len(sysargv)>6 and int(sysargv[6])==0: plot     = False

    # le modèle à traiter (SEIR1R2 or SEIR1R2D)
    if modeleString == 'SEIR1R2':
        fit = fitProcessSEIR1R2
    elif modeleString == 'SEIR1R2D':
        fit = fitProcessSEIR1R2D
    else:
        print('Wrong EDO model, only SEIR1R2 or SEIR1R2D available!')
        exit(1)


    # fit avec une seule période
    ##################################
    # nbperiodes = 1 # ne peut pas etre changé
    # decalage1P = 0 # ne peut pas etre changé
    
    # model_allinone, ListeTextParamPlace_allinone, liste_pd_allinone, data_deriv_allinone, model_deriv_allinone, _, _ = \
    #         fit([places, nbperiodes, decalage1P, UKF_filt01, 0, 0])

    # ListeTestPlace = []
    # for indexplace in range(len(listplaces)):
    #     texteplace = ''
    #     for texte in ListeTextParamPlace_allinone[indexplace]:
    #         texteplace += '\n' + texte
    #     ListeTestPlace.append(texteplace)

    # PlotPlace(modeleString, data_deriv_allinone, model_deriv_allinone, listplaces, decalage1P, UKF_filt, FrDatabase, 'All-In-One', ListeTestPlace)


    # fit avec 3 périodes + décalage
    ##################################
    nbperiodes = -1

    model_piecewise, ListeTextParamPlace_piecewise, liste_pd_piecewise, data_deriv_piecewise, model_deriv_piecewise, _, _ = \
            fit([places, nbperiodes, decalage3P, UKF_filt01, 0, 0])
    
    ListeTestPlace = []
    for indexplace in range(len(listplaces)):
        texteplace = ''
        for texte in ListeTextParamPlace_piecewise[indexplace]:
            texteplace += '\n' + texte
        ListeTestPlace.append(texteplace)

    # PlotPlace(modeleString, data_deriv_piecewise, model_deriv_piecewise, listplaces, decalage3P, UKF_filt, FrDatabase, 'Piecewise', ListeTestPlace)


    # Plot the two strategies in one figure
    ##################################
    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
        else:
            placefull   = place

        # Repertoire des figures
        repertoire = getRepertoire(UKF_filt, './figures/'+modeleString+'_UKFilt/'+placefull, './figures/'+modeleString+'/'+placefull)
        prefFig    = repertoire+'/'

        # Preparation plot pandas
        listheader = list(liste_pd_piecewise[indexplace])
        
        # on ajoute les dérivées numériques des cas et des morts
        liste_pd_piecewise[indexplace]['dcases']           = liste_pd_piecewise[indexplace][listheader[0]].diff()
        liste_pd_piecewise[indexplace]['ddeaths']          = liste_pd_piecewise[indexplace][listheader[1]].diff()
        liste_pd_piecewise[indexplace]['dcasesplusdeaths'] = liste_pd_piecewise[indexplace][listheader[2]].diff()
        longueur = len(liste_pd_piecewise[indexplace].loc[:, ('dcases')])

        if FrDatabase==True:
            DatesString = getDates('France', verbose)
        else:
            DatesString = getDates(place, verbose)

        # on ajoute les dérivées numériques des cas et des morts
        liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise')]          = model_deriv_piecewise[indexplace, 0:longueur, 0]
        liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise_residual')] = liste_pd_piecewise[indexplace].loc[:, ('mc_piecewise')]-liste_pd_piecewise[indexplace].loc[:, ('dcasesplusdeaths')]
        if modeleString=='SEIR1R2D':
            liste_pd_piecewise[indexplace].loc[:, ('md_piecewise')] = model_deriv_piecewise[indexplace, 0:longueur, 1]
            liste_pd_piecewise[indexplace].loc[:, ('md_piecewise_residual')] = liste_pd_piecewise[indexplace].loc[:, ('md_piecewise')]-liste_pd_piecewise[indexplace].loc[:, ('ddeaths')]

        # Dessin des dérivées
        filename = prefFig + 'Diff_FitPiecewise_Shift' + str(decalage3P) + '.png'
        title    = place + ' - Shift=' + str(decalage3P) + ' day(s)'
        if modeleString=='SEIR1R2':
            listPlots = ['dcasesplusdeaths', 'mc_piecewise']
        if modeleString=='SEIR1R2D':
            listPlots = ['dcases', 'mc_piecewise', 'ddeaths', 'md_piecewise']
        PlotFitPiecewise(modeleString, liste_pd_piecewise[indexplace], title, filename, y=listPlots, Dates=DatesString, textannotation=ListeTestPlace[indexplace])

        # Dessin des résidus
        filename = prefFig + 'Diff_FitPiecewiseResiduals_Shift' + str(decalage3P) + '.png'
        title    = place + ' - Shift=' + str(decalage3P) + ' day(s)'
        if modeleString=='SEIR1R2':
            listPlots = ['mc_piecewise_residual']
        if modeleString=='SEIR1R2D':
            listPlots = ['mc_piecewise_residual', 'md_piecewise_residual']
        PlotFitPiecewiseResidual(modeleString, liste_pd_piecewise[indexplace], title, filename, y=listPlots, Dates=DatesString, textannotation=ListeTestPlace[indexplace])

        input('temp')

    # Plot the two SEIR1R2
    ##################################

    for indexplace, place in enumerate(listplaces):

        # Get the full name of the place to process, and the special dates corresponding to the place
        if FrDatabase == True: 
            placefull   = 'France-' + place
        else:
            placefull   = place

        # Repertoire des figures
        repertoire = getRepertoire(UKF_filt, './figures/' + modeleString + '_UKFilt/'+placefull, './figures/' + modeleString + '/' + placefull)
        prefFig = repertoire + '/'

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
        if modeleString == 'SEIR1R2D':
            liste_pd_piecewise[indexplace].loc[:, ('Fp')]  = model_piecewise [indexplace, :, 5]

        title    = place + ' - Shift=' + str(decalage3P) + ' day(s)'
        listePlot=['Ep', 'Ip', 'R1p']
        if modeleString=='SEIR1R2D':
            listePlot.append('Fp')
        filename = prefFig + ''.join(map(str, listePlot)) + '_piecewise_Shift' + str(decalage3P) + '.png'
        PlotModel(modeleString, liste_pd_piecewise[indexplace], title, filename, y=listePlot, Dates=DatesString, textannotation=ListeTestPlace[indexplace])

        title    = place + ' - Shift=' + str(decalage3P) + ' day(s)'
        listePlot=['R1p']
        if modeleString=='SEIR1R2D':
            listePlot.append('Fp')
        filename = prefFig + ''.join(map(str, listePlot)) + '_piecewise_Shift' + str(decalage3P) + '.png'
        PlotModel(modeleString, liste_pd_piecewise[indexplace], title, filename, y=listePlot, Dates=DatesString, textannotation=ListeTestPlace[indexplace])


def PlotModel(modeleString, pd, titre, filenameFig, y, Dates=None, textannotation=''):
    
    if len(y)==0 or y is None: pass

    # juste pour récupérer des couleurs sohérentes
    if modeleString=='SEIR1R2':
        solveur = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
    else:
        solveur = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)

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

def PlotFitPiecewise(modeleString, pd, titre, filenameFig, y, Dates=None, textannotation=''):

    if len(y)==0 or y is None: pass

    # on supprime la première ligne de données (car ce sont des dérivées)
    pd = pd.iloc[1:]
    
    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # Dessin des courbes
    if modeleString=='SEIR1R2':
        solveur    = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
        colori     = solveur.modele.getColor(3)
        linestyles = ['-', '-']#, '-']
        markers    = ['x', '']#, '']
        colors     = [colori, colori]#, 'blue']
        linewidths = [0.5, 1.5]#, 1.5]
        alphas     = [1.0, 0.7]#, 0.55]
        labels     = [r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise']
    else:
        solveur    = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)
        colori     = [solveur.modele.getColor(3), solveur.modele.getColor(5)]
        linestyles = ['-', '-', '-', '-']#, '-']
        markers    = ['x', '', '+', '']#, '']
        colors     = [colori[0], colori[0], colori[1], colori[1]]#, 'blue']
        linewidths = [0.5, 1.5, 0.5, 1.5]#, 1.5]
        alphas     = [1.0, 0.7, 1.0, 0.7]#, 0.55]
        labels     = [r'$\frac{\partial R^1(n)}{\partial n}$', r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise', \
                      r'$\frac{\partial F(n)}{\partial n}$',   r'$\frac{\partial F(t)}{\partial t}$ - Piecewise']
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


def PlotFitPiecewiseResidual(modeleString, pd, titre, filenameFig, y, Dates=None, textannotation=''):

    if len(y)==0 or y is None: pass

    # on supprime la première ligne de données (car ce sont des dérivées)
    pd = pd.iloc[1:]
    
    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    # Dessin des courbes
    if modeleString=='SEIR1R2':
        solveur    = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
        colori     = solveur.modele.getColor(3)
        linestyles = ['-']
        markers    = ['x']
        colors     = [colori]
        linewidths = [0.5]
        alphas     = [1.0]
        labels     = [r'$\frac{\partial R^1(n)}{\partial n}-\frac{\partial R^1(t)}{\partial t}$ - Piecewise']
    else:
        solveur    = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)
        colori     = [solveur.modele.getColor(3), solveur.modele.getColor(5)]
        linestyles = ['-', '-']
        markers    = ['x', '+']
        colors     = [colori[0], colori[1]]
        linewidths = [0.5, 0.5]
        alphas     = [1.0, 1.0]
        labels     = [r'$\frac{\partial R^1(n)}{\partial n}-\frac{\partial R^1(t)}{\partial t}$ - Piecewise', \
                      r'$\frac{\partial F(n)}{\partial n}-\frac{\partial F(t)}{\partial t}$ - Piecewise']
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

# def PlotAll(data_deriv_allinone, model_deriv_allinone, model_deriv_piecewise, title, filename):

#     fig = plt.figure(facecolor='w', figsize=figsize)
#     ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#     # juste pour récupérer la bonne coukeur
#     solveur = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
#     mod   = solveur.modele
#     color = mod.getColor(3)

#     time = np.linspace(1, len(data_deriv_allinone[1:, 0]), len(data_deriv_allinone[1:, 0]))
#     ax.plot(time, data_deriv_allinone     [1:, 0], label=r'$\frac{\partial R^1(n)}{\partial n}$',              color=color,  alpha=1.0,  lw=0.5, marker='x')
#     ax.plot(time, model_deriv_allinone [1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - All-in-one', color='blue', alpha=0.55, lw=1.5)
#     ax.plot(time, model_deriv_piecewise[1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - Piecewise',  color=color,  alpha=0.7,  lw=1.5)

#     ax.set_xlabel('Time (days)')
#     ax.yaxis.set_tick_params(length=0)
#     ax.xaxis.set_tick_params(length=0)
#     ax.grid(b=True, which='major', c='w', lw=1, ls='-')
    
#     legend = ax.legend()
#     legend.get_frame().set_alpha(0.5)
#     for spine in ('top', 'right', 'bottom', 'left'):
#         ax.spines[spine].set_visible(False)

#     plt.title(title)
#     plt.savefig(filename, dpi=dpi)
#     plt.close()


# def PlotPlace(modeleString, data_deriv, model_derive, listplaces, decalage, UKF_filt, FrDatabase, ch, listetextannotation):

#     # Plot
#     for indexplace, place in enumerate(listplaces):

#         # Get the full name of the place to process, and the special dates corresponding to the place
#         if FrDatabase == True: 
#             placefull   = 'France-' + place
#         else:
#             placefull   = place

#         # Repertoire des figures
#         repertoire = getRepertoire(UKF_filt, './figures/' + modeleString + '_UKFilt/'+placefull, './figures/' + modeleString + '/' + placefull)
#         prefFig = repertoire + '/'

#         filename = prefFig + 'Diff_' + ch + '_Shift' + str(decalage) + '.png'
#         title    = place + ' - Shift=' + str(decalage) + ' day(s)'
#         PlotFit(modeleString, data_deriv[indexplace], model_derive[indexplace], title, filename, ch, listetextannotation[indexplace])


# def PlotFit(modeleString, data_deriv, model_deriv, title, filename, ch, textannotation=''):

#     fig = plt.figure(facecolor='w', figsize=figsize)
#     ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#     time = np.linspace(1, len(data_deriv[1:]), len(data_deriv[1:]))

#     if modeleString=='SEIR1R2':
#         solveur = SolveEDO_SEIR1R2(N=1, dt=1, verbose=0)
#         colori  = [solveur.modele.getColor(3)]
#     elif modeleString=='SEIR1R2D':
#         solveur = SolveEDO_SEIR1R2D(N=1, dt=1, verbose=0)
#         colori  = [solveur.modele.getColor(3), solveur.modele.getColor(5)]
#     else:
#         print('PlotFit Impossible')
#         exit(1)

#     # les plots
#     ax.plot(time, data_deriv [1:, 0], label=r'$\frac{\partial R^1(n)}{\partial n}$ - ' + ch, color=colori[0], alpha=1.0, lw=0.5, marker='x')
#     ax.plot(time, model_deriv[1:, 0], label=r'$\frac{\partial R^1(t)}{\partial t}$ - ' + ch, color=colori[0], alpha=1.0, lw=1.5)
#     if modeleString=='SEIR1R2D':
#         ax.plot(time, data_deriv [1:, 1], label=r'$\frac{\partial D(n)}{\partial n}$ - '   + ch, color=colori[1], alpha=1.0, lw=0.5, marker='+')
#         ax.plot(time, model_deriv[1:, 1], label=r'$\frac{\partial D(t)}{\partial t}$ - '   + ch, color=colori[1], alpha=1.0, lw=1.5)

#     # ajout d'un text d'annotation
#     if textannotation != '':
#         bbox=dict(boxstyle='round4, pad=.3', alpha=0.5, color='xkcd:light turquoise') #, facecolor='0.8', edgecolor='xkcd:light turquoise', lw=0.5)
#         ax.annotate(textannotation, xy=(ax.get_xlim()[1], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])/2.6), fontsize=6, bbox=bbox, ha="left", va="center") 

#     ax.set_xlabel('Time (days)')
#     ax.yaxis.set_tick_params(length=0)
#     ax.xaxis.set_tick_params(length=0)
#     ax.grid(b=True, which='major', c='w', lw=1, ls='-')
    
#     legend = ax.legend()
#     legend.get_frame().set_alpha(0.5)
#     for spine in ('top', 'right', 'bottom', 'left'):
#         ax.spines[spine].set_visible(False)

#     plt.tight_layout(rect=(0, 0, 1., 0.95))
#     plt.title(title)
#     plt.savefig(filename, dpi=dpi)
#     plt.close()



if __name__ == '__main__':
    main(sys.argv)
