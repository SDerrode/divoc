# Projet divoc

Collection of some Jupyter Notebooks and Python3 scripts to test non-linear filtering techniques (in particular UKF) to process Covid Data collected on the web. 

The programs are based on the Kalman-like filters implemented by Roger R Labbe Jr. in the [FilterPy](https://filterpy.readthedocs.io/en/latest/index.html) python module. 

## Notebooks repository

  - *readPlotDataGouvFr.ipynb*: Data reading from the website [data.gouv.fr](data.gouv.fr) (more specifically [these ones](https://static.data.gouv.fr/resources/donnees-hospitalieres-relatives-a-lepidemie-de-covid-19/20200327-154414/metadonnees-donnees-hospitalieres-covid19.csv)) and production of some figures for verification.

  - *readPlotDataEurope.ipynb*:  Data reading from the website [European Center for Disease Prevention and Control](https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide) and production of some figures for verification.

Start reading online now by clicking the binder or Azure badge below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SDerrode/divoc/master)

[![Binder](http://mybinder.org/badge.svg)](https://beta.mybinder.org/v2/gh/rlabbe/Kalman-and-Bayesian-Filters-in-Python/master)
<a href="https://notebooks.azure.com/import/gh/rlabbe/Kalman-and-Bayesian-Filters-in-Python"><img src="https://notebooks.azure.com/launch.png" /></a>


## NL-Filtering

  - *UKF_SEIR.py*: Main program to test UKF on SEIR-type models. Only Bacaer model implemented

  - *SEIR_Bacaer.py*: Class containing the description of the model described in the paper of Nicolas Bacaer (IRD). See file for exact referencing. 

  - *Covid_SpecialDates*: Class to set-up and manage special Covid dates for a country (confinement date(s), deconfinement date(s)
, others dates. Mainly used to annotate plots.

  - *common.py*: Collection of functions for reading data on the web and for plotting

