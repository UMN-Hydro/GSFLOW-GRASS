# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 15:26:24 2017

@author: gcng
"""

# plot_gsflow_csv.m
# 
# List of StatVarNames: see Table 12 of GSFLOW manual and 
#  create_table_gsflowcsv.m

import sys
import platform
import numpy as np
from matplotlib import pyplot as plt
import settings_test
import pandas as pd
import datetime as dt
import matplotlib.dates as mdates
import gsflow_csv_table as gvar  # all variable names, units, and descriptions

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'
 
# *** SPECIFY PLOT-VARIABLES HERE
 
# -- variables to plot (see list in gsflow_csv_table.py):
PlotVar = []
PlotVar.append('basinppt')
PlotVar.append('basinactet')

# -- data files
gsflow_csv_fil = settings_test.PRMSoutput_dir + slashstr + 'gsflow.csv'  # gsflow time series data

# -- save figure to this file
savefigfile = 'fig.png'

#%%


descr = []
unit_prev = []
for jj in range(len(PlotVar)):
    for ii in range(len(gvar.varname)):
        if gvar.varname[ii] == PlotVar[jj]:
            unit = gvar.unit[ii]
            descr.append(gvar.descr[ii])
            if (jj > 0) and (unit != unit_prev):
                    print "Error! Plot variables do not have same units.  Exiting..."
                    sys.exit()
            unit_prev = unit
            break
    

## Read in data

# header{1,NVars}: variable name
# data{NVars}: all data 
data = pd.read_csv(gsflow_csv_fil);   # 
    
dateList = [dt.datetime.strptime(date, '%m/%d/%Y').date() for date in data['Date']]


# - plot data
plt.figure()
#plt.close("all")
for ii in range(len(PlotVar)):
    plt.plot(dateList, data[PlotVar[ii]])
#    plt.plot(range(len(data)), data[PlotVar[ii]])
    plt.show
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
#    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()    
    plt.ylabel(unit)
    
plt.legend(PlotVar)
plt.title(settings_test.PROJ_CODE)

plt.savefig(savefigfile, dpi = 300)


print '\n-------------------'
print 'Plotting: '
for ii in range(len(PlotVar)):
    print '- ' + PlotVar[ii] + ' (' + unit + '): ' + descr[ii]
print '-------------------\n'


