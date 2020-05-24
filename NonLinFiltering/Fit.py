#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy             as np
import matplotlib.pyplot as plt
from datetime       import datetime, timedelta

from common         import readDataEurope, getDates, Plot, addDaystoStrDate
from common         import getLowerDateFromString, getNbDaysBetweenDateFromString
from ProcessSEIR1R2 import fit

dpi     = 150    # plot resolution of saved figures
figsize = (8, 4) # figure's size (width, height)

def main():
    verbose       = 0
    plot          = 0

    #countries     = 'France,Spain,Italy,United_Kingdom,Germany,Belgium'
    countries     = 'France,Spain,Belgium'
    listcountries = list(countries.split(','))

    # fit avec une seule période
    ##################################
    nbperiodes = 1
    decalage   = 0
    surplus    = 0

    pd_exerpt, data_deriv_0, moldelR1_deriv_0 = fitcountry(countries, nbperiodes, decalage, surplus, verbose, plot, '1period')

    # fit avec 3 périodes + décalage
    ##################################
    nbperiodes = 3
    decalage   = 10
    surplus    = 0

    d_exerpt, data_deriv_1, moldelR1_deriv_1 = fitcountry(countries, nbperiodes, decalage, surplus, verbose, plot, '3periods')


    # Plot all in one figure
    for indexcountry in range(len(listcountries)):
        country  = listcountries[indexcountry]
        filename = './figures/DiffR1_BothFit_' + country + '_Shift' + str(decalage) + '.png'
        title    = country + ' - Shift=' + str(decalage) + ' day(s)'
        PlotAll(data_deriv_0[indexcountry], moldelR1_deriv_0[indexcountry], data_deriv_1[indexcountry], moldelR1_deriv_1[indexcountry], title, filename)


def PlotAll(data_deriv_0, model_deriv_0, data_deriv_1, model_deriv_1, title, filename):
    fig = plt.figure(facecolor='w', figsize=figsize)
    ax  = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

    time = np.linspace(1, len(data_deriv_0[1:]), len(data_deriv_0[1:]))
    ax.plot(time, data_deriv_0[1:],  label='Instant cases - 1period',                          color='green', alpha=1.0, lw=0.5, marker='x')
    ax.plot(time, model_deriv_0[1:], label=r'$\frac{\partial R^1(t)}{\partial t}$ - 1period',  color='blue', alpha=0.4, lw=1.3)
    # ax.plot(time, data_deriv_1[1:],  label='Instant cases - 3periods',                         color='blue',   alpha=1.0, lw=0.5, marker='x')
    ax.plot(time, model_deriv_1[1:], label=r'$\frac{\partial R^1(t)}{\partial t}$ - 3periods', color='green',   alpha=0.7, lw=1.3)

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
    ax.plot(time, data_deriv[1:],  label='Instant cases - '  + ch,                        color='green', alpha=1.0, lw=0.5, marker='x')
    ax.plot(time, model_deriv[1:], label=r'$\frac{\partial R^1(t)}{\partial t}$ - ' + ch, color='green', alpha=1.0, lw=2.0)

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

def fitcountry(countries, nbperiodes, decalage, surplus, verbose, plot, ch):

    pd_exerpt, data_deriv, moldelR1_deriv = fit([countries, nbperiodes, decalage, surplus, verbose, plot])
    dataLength = pd_exerpt.shape[0]
    # print('dataLength=', dataLength)
    # print('shape data_deriv=', np.shape(data_deriv))

    # Plot
    listcountries = list(countries.split(','))
    for indexcountry in range(len(listcountries)):
        country  = listcountries[indexcountry]
        filename = './figures/DiffR1_' + country + '_' + ch + '_Shift' + str(decalage) + '.png'
        title    = country + ' - Shift=' + str(decalage) + ' day(s)'
        PlotFit(data_deriv[indexcountry], moldelR1_deriv[indexcountry], title, ch, filename)

    return pd_exerpt, data_deriv, moldelR1_deriv


if __name__ == '__main__':
    main()


