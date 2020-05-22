#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from datetime  import datetime, timedelta
from common    import readDataEurope, readDataGouvFr, getDates, PlotData, addDaystoStrDate

# constante
fileLocalCopy = True         # if we upload the file from the url (to get latest results) or from a local copy file

def main():
    """
        Program to perform UKF filtering on Covid Data.
 
        :Example:

        >> python3 PlotDataCovid.py 
        >> python3 PlotDataCovid.py Italy
        >> python3 PlotDataCovid.py Italy 1 
        >> python3 PlotDataCovid.py France,Germany 1 # Shortcut for processing the two countries successively
       
        argv[1] : Country (or list separeted by ','), or Continent, or 'World'. Default: France 
        argv[2] : Verbose level (debug: 3, ..., almost mute: 0). Default: 1
    """

    print('Command line : ', sys.argv, flush=True)
    if len(sys.argv) > 3:
        print('  CAUTION : bad number of arguments - see help')
        exit(1)

    # Default value for parameters
    listcountries = ['France']
    verbose       = 1
    
    # Parameters from argv
    if len(sys.argv)>1: listcountries = list(sys.argv[1].split(','))
    if len(sys.argv)>2: verbose       = int(sys.argv[2])

    dt = 1
    readStartDateStr = "2020-03-01" 
    readStopDateStr  = None
    
    for country in listcountries:

        print('PROCESSING of ', country, ' in ', listcountries)
        prefFig = './figures/Data_' + country

        # These are the date of confinement and deconfinement + other. 
        # See function getDates on how to add or delete dates to put the focus on
        Dates = getDates(country, verbose)
        
        # Lecture des données et copy of the observation
        #############################################################################
        pd_exerpt, HeadData, N, readStartDateStr, readStopDateStr = readDataEurope(country, readStartDateStr, readStopDateStr, fileLocalCopy, verbose=0)
        readStartDate = datetime.strptime(readStartDateStr, "%Y-%m-%d")
        readStopDate  = datetime.strptime(readStopDateStr,  "%Y-%m-%d")
        dataLength = pd_exerpt.shape[0]
        if verbose>1:
            print(readStartDateStr, readStopDateStr)
            print(readStartDate, readStopDate)
            #input('pause')

        print('tail=', pd_exerpt.tail())

        PlotData(pd_exerpt, titre=country, filenameFig=prefFig+'.png', y=HeadData, color='green', Dates=Dates)



if __name__ == '__main__':
    main()
