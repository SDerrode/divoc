
This file describes how to regenerate all the figures in the preprint entitled: **Piecewise estimation of R0 by a simple SEIR model. Application to COVID-19 in French regions and departments.** If I have time ;-), I plan to include all the commands in a Jupiter script. 

The commands exposed below give an idea on how to play with the scripts. Please refer to the first lines of the scripts to have more information on the possible arguments. 

To run the scripts, you will need a number of python library, including among the less popular ones : datetime, dateutil, contextily, geopandas, warnings, requests, zipfile, io, pathlib, sklearn, filterpy, lmfit, lmfit

In the preprint, we only deal with the SEIR1R2 model. I have started to work on a variation of this model called SEIR1R2D, which include dead in another way, but the model has not been tested very much, so be careful with this model!


**Figure 1:**
    
    >> python Fit.py FRANCE,MetropoleD+ 0 SEIR1R2 18 0 1 1

Files are generated in directory *./figures/SEIR1R2/France/sexe_0_delay_18/*

Remark: Here *FRANCE,MetropoleD+* means that we add the data for all the departements in the database, to get results for France. You can also try *python Fit.py FRANCE,R53,D06,R84+ 0 SEIR1R2 18 0 1 1* to get the result for the department 06, for all the departments in the region 53 and for the region 84 (ie the sum of data for all the departments in the region 84).

**Figure 2:**

    >> python PlotR0Est_TS.py

Files are generated in directory *./figures/SEIR1R2/R0Est_TS/sexe_0_delay_18/*

Remark: The parameters (e.g. the start date of analysis, sex...) have to be modified in the python script directly (I plan to include a command line arguments).

Remark: It can take several minutes to get the figures.

**Figure 3:**

    >> python PlotTimeShift.py FRANCE,MetropoleD+ 0 SEIR1R2 0 14,23 1 1

Files are generated in directory *./figures/France/sexe_0_delay_14_23/*

Remark: It can take several minutes to get the figures.


**Figure 4:** 4 France maps at region scale
    
    >> python PlotTimeShift.py FRANCE,MetropoleR+ 0 SEIR1R2 0 18,19 0 1
    >> python MapFranceR0.py REG R0Moyen_18_19.csv

Remark: The first command generates a csv file in directory *./figures/SEIR1R2/TimeShift/*. This one used by the second call to plot the 4 maps in the same directory.

Remark: The first call can take up to 3/4 minutes.


**Figure 5:** 4 maps at department scale
    
    >> python PlotTimeShift.py FRANCE,MetropoleD 0 SEIR1R2 0 18,19 0 1
    >> python MapFranceR0.py DPT R0Moyen_18_19.csv

Remark: The first command generates a csv file in directory *./figures/SEIR1R2/TimeShift*. This one used by the second call to plot the 4 maps in the same directory. Be aware that this second call will delete previous files, so don't forget to rename the previous file if you want to keep all the maps for regions and departments.

Remark: The first call can take up to 20 minutes.
