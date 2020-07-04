# Projet divoc

Collection of Python3 scripts to process Covid Data collected on the web. 

<!--The programs are based on the Kalman-like filters implemented by Roger R Labbe Jr. in the [FilterPy](https://filterpy.readthedocs.io/en/latest/index.html) python module. -->


Here is a result of using the script _Fit.py_ when France is selected (see text below).

![Image fit France](./France_DiffR1_BothFit.png "Fit result for France")

#TO REWRITE ENTERELY!!!



## NonLinFiltering repository

The explanations are given for the SEIR1R2 model (ie files which filename is terminating by SEIR1R2). I will add other models (for example the model extension called SEIRR1R2) which usae is more or less the same.

### Presentations of python files

- The files _common.py_, _Covid\_SpecialDates.py_ are for internal use:
   
    + *common.py*: Collection of functions for reading data on the web (or in data repository), for plotting, and many others.
    + *Covid_SpecialDates*: Class to set-up and manage special Covid dates for a country (confinement date(s), deconfinement date(s), others dates).

- The file _SEIR1R2.py_ is the implementation of the pandemic model used by N. Bacaert in his paper (see file for exact referencing). It is used in connection with _SolveEDO\_SEIR1R2.py_.
      
- The file _SolveEDO\_SEIR1R2.py_ contains a class used by programs to solve the EDO SEIR1R2, but it is also runnable to simulate the behavior of a SEIR1R2 model (see usage below).

- The files _PlotDataCovid.py_, _ProcessSEIR1R2.py_, _Fit\_SEIR1R2.py_ and _PlotTS\_SEIR1R2.py_ are runnable scripts (see usage below).

### Usages

1. The first script to test is the _PlotDataCovid.py_ one, which generates plots of cases and deaths for whatever country, or for whatever department of France. Typical usages are:

> python3 PlotDataCovid.py Germany
> 
> python3 PlotDataCovid.py France,Spain,Italy,United_Kingdom,Germany,Belgium
> 
> python3 PlotDataCovid.py France,69,75


See what kind of plots is generated in the _figures_ repository.


2. Then you will be interested in running the main program in _SolveEDO\_SEIR1R2.py_:

> python3 SolveEDO_SEIR1R2.py

Modify the parameters in the main program (at the end of the file) to get different simulations. Here is an example of the plot generated (see _figures_ repo):

![SEIR1R2 simulation](./SEIR1R2model_01234.png "SEIR1R2 simulation")

3. The next step is to run the script _ProcessSEIR1R2.py_. This is the main script, with (nearly) all the parameters editable by command arguments:

> python3 ProcessSEIR1R2.py France,69 1 0 0 1 1
> 
> python3 ProcessSEIR1R2.py Italy,Spain 2 10 0 1 1

First example: the program deals with the French department 69 ('Auvergne-Rhone-Alpes'), and try to find the parameters that best fit all the data (all the data is only one period). 
    
Second example: You can also deals with severall periods, where a shift of 10 days with respect to the the lock-down date of the considered countries.

4. Eventually, you can run these scripts for processing several countries or French counties at one

> python3 Fit_SEIR1R2.py

The main goal of this program is to call the previous program to generate plots that merge all the periods. 

You have to modify manually the parameters at the beginning of the script (change the countries or the department, change the time shift from 10 to 5 for example).
See what happens in the _figures_ repository.
    
> python3 PlotTS_SEIR1R2.py

This last call runs _ProcessSEIR1R2.py_ severall times, varying the time-shift (the range of values must be set by hand in the _ProcessSEIR1R2.py_ file) and plot parameters (a, b, c, f, R_0) evolution according to time-shift. The program also saves the EQM between the derivative of R1 data and the derivative of model's R1, to get a possible idea on the "optimal" time-shift.

## Contact

stephane dot derrode at ec-lyon dot fr

