
**Figure 1:**
    
    >> python Fit.py FRANCE,MetropoleD+ 0 SEIR1R2  13 0 1 1

File: *./figures/SEIR1R2/France-MetropoleD+/sexe_0_shift_13/Fit_13_Diff_Piecewise.png*

**Figure 2:**

    >> python PlotR0Est_TS.py

File: *./figures/SEIR1R2/R0Est_TS_sexe_0_shift_13/ROEst_TS.csv*

Remark: A few parameters (e.g. the start date of analysis, sex...) can be modified in the python script.

Remark: It can take several minutes to get the fig.


**Figure 3:** 4 maps at region scale
    
    >> python PlotTimeShift.py FRANCE,MetropoleR+ 0 SEIR1R2  0 13,14 0 1
    >> python MapFranceR0.py REG R0Moyen_13_14.csv

Remark: The first command generates a file in repo *./figures/SEIR1R2/TimeShift*. This one used by the second call to plot the 4 maps in the same repo.

Remark: It can take up to 2/3 minutes to get the fig.


**Figure 4:** 4 maps at department scale
    
    >> python PlotTimeShift.py FRANCE,MetropoleD 0 SEIR1R2  0 13,14 0 1
    >> python MapFranceR0.py DPT R0Moyen_13_14.csv

Remark: The first command generates a file in repo *./figures/SEIR1R2/TimeShift*. This one used by the second call to plot the 4 maps in the same repo.

Remark: It can take up to 20/25 minutes to get the fig.


**Figure 5:** similar to Figure 1

    >> python Fit.py FRANCE,D06 0 SEIR1R2  13 0 1 1

File: *./figures/SEIR1R2/France-D06/sexe_0_shift_13/Fit_13_Diff_Piecewise.png*
