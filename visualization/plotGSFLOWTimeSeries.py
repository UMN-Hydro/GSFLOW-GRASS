#! /usr/bin/env python

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
import pandas as pd
import datetime as dt
import matplotlib.dates as mdates

if platform.system() == 'Linux' or platform.system() == 'Darwin':
    slashstr = '/'
else:
    slashstr = '\\'

# add path containing readSettings.py
sys.path.append('..' + slashstr + 'Run')

# Read in user-specified settings
from readSettings import Settings
# Set input file
if len(sys.argv) < 2:
    settings_input_file = 'settings.ini'
    print 'Using default input file: ' + settings_input_file
else:
    settings_input_file = sys.argv[1]
    print 'Using specified input file: ' + settings_input_file
Settings = Settings(settings_input_file)

if Settings.GSFLOW_ver == '1.2.0':
#    sys.exit('Need to update Time Series script to work for GSFLOW 1.2.2 outputs! Exiting...')
    import GSFLOWcsvTable as gvar  # all variable names, units, and descriptions
elif Settings.GSFLOW_ver == '1.2.2':
    import GSFLOWcsvTable_1p2p2 as gvar  # all variable names, units, and descriptions



#%% *** SET THE FOLLOWING *****************************************************
 
# *** enter variables to plot 
 # (see list in gsflow_csv_table.py for GSFLOW 1.2.0 output and 
 #gsflow_csv_table_1p2p2.py for GSFLOW 1.2.2 output and ):
PlotVar = []
if Settings.GSFLOW_ver == '1.2.0':
    PlotVar.append('basinstrmflow')
    PlotVar.append('uzf_recharge')
    PlotVar.append('basinsroff')
    PlotVar.append('gwflow2strms')
    #PlotVar.append('basininterflow')
    #PlotVar.append('streambed_loss')
#    PlotVar.append('basinpervet')
    
    PlotVar.append('basinppt')
    #PlotVar.append('basinactet')

elif Settings.GSFLOW_ver == '1.2.2': 
    PlotVar.append('StreamOut_Q')
    PlotVar.append('RechargeUnsat2Sat_Q')
    PlotVar.append('HortSroff2Stream_Q')
    PlotVar.append('StreamExchng2Sat_Q')
    #PlotVar.append('Interflow2Stream_Q')
    
    PlotVar.append('Precip_Q')
    #PlotVar.append('CapET_Q')


# *** save figure to this file
savefigfile = 'fig.png'

#%% *** CHANGE FILE NAMES IF NEEDED *******************************************
# (default is to use entries from Settings File) 
gsflow_csv_fil = Settings.PRMSoutput_dir + slashstr + 'gsflow.csv'  # gsflow time series data
plot_title = Settings.PROJ_CODE

#  *** CHANGE BASIN AREA AS NEEDED (SEE prms.out FILE FOR AREA IN [acres]) ***
# (default is to use entries from Settings File) 
HRUfil = Settings.GISinput_dir + slashstr + 'HRUs.txt'
HRUdata = pd.read_csv(HRUfil)
HRUarea = HRUdata['hru_area']  # [acre]
basin_area = sum(HRUarea) * 4046.85642  # acre -> m2


#%%

# make sure basinppt and basinacet are listed last, because of twin y-axis
PlotVar0 = PlotVar[:]
ctr = 1
ctr_end = len(PlotVar)
for ii in range(len(PlotVar)):
    if (PlotVar[ii] == 'basinppt') or (PlotVar[ii] == 'basinactet') or (PlotVar[ii] == 'Precip_Q') or (PlotVar[ii] == 'CapET_Q'):
        PlotVar0[ctr_end-1] = PlotVar[ii]
        ctr_end = ctr_end - 1
    else:
        PlotVar0[ctr-1] = PlotVar[ii]
        ctr = ctr + 1
PlotVar = PlotVar0[:]

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
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

# convert volume units [m^3] to length [mm]
if unit[0:2] == 'm^3':
    conv = 1. / basin_area * 1000.
    unit0 = 'mm' 
    if len(unit) > 3:
        unit0 + unit[3:]
    unit = unit0[:]
    

#plt.close("all")
for ii in range(len(PlotVar)):
    
    if (PlotVar[ii] == 'basinppt') or (PlotVar[ii] == 'basinactet') or (PlotVar[ii] == 'Precip_Q') or (PlotVar[ii] == 'CapET_Q'):
        ln0 = ax2.plot(dateList, data[PlotVar[ii]], '--k')
        ax2.set_ylabel(PlotVar[ii] + ' ' + unit)
    else:
        ln0 = ax1.plot(dateList, data[PlotVar[ii]])            
        ax1.set_ylabel(unit)

    if ii == 0:
        lns = ln0
    else:
        lns = lns + ln0
     
#    plt.plot(range(len(data)), data[PlotVar[ii]])
    plt.show
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
#    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()    

#lns = lns1+lns2+lns3
#labs = [l.get_label() for l in lns]
ax1.legend(lns, PlotVar, loc=0)
    
#plt.legend(PlotVar)
plt.title(plot_title)

plt.show()

plt.savefig(savefigfile, dpi = 300)


print '\n-------------------'
print 'Plotting: '
for ii in range(len(PlotVar)):
    print '- ' + PlotVar[ii] + ' (' + unit + '): ' + descr[ii]
print '-------------------\n'


