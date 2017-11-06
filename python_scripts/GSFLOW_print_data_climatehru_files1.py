# Based on: GSFLOW_print_data_climatehru_files1.m

# GSFLOW_print_climate_hru_files1.m
# 11/24/16
#
# based on: PRMS_print_climate_hru_files2.m

# Creates climate data inputs files to run PRMS for climate_hru mode using 
# by uniformly applying data from a single weather station:
#   Specify precip, tmax, tmin, and solrad for each HRU for each
#   simulation day.  Can have start day earlier and end day later than
#   simulation period.
# 

import pandas as pd
import numpy as np
import settings_test
import platform

if platform.system() == 'Linux':
    slashstr = '/'
else:
    slashstr = '\\'



# *** CUSTOMIZE TO YOUR COMPUTER! *****************************************

# directory for PRMS input files (include slash ('/') at end) - cbh files
# outputed to here
#PRMSinput_dir = GSFLOW_DIR + '/inputs/PRMS/'
PRMSinput_dir = settings_test.PRMSinput_dir

# PRMSinput_dir = '/home/gcng/workspace/ProjectFiles/AndesWaterResources/GSFLOW/ChimTest'

# GIS data (for nhru)
in_GISdata_dir  = settings_test.GISinput_dir
HRUfil = in_GISdata_dir + slashstr + 'HRUs_tmp.txt'

# *************************************************************************

# Project-specific entries ->

# -- Number of HRU's over which to apply data uniformly
# nhru should be generated dynamically (from GIS data)
HRUdata = pd.read_csv(HRUfil)
nhru = HRUdata.shape[0]

## -- Set any description of data here (example: location, units)
descr_str = settings_test.PROJ_CODE

# -- Read in daily values from 1 station for:
#   - precip (check 'precip_units')
#   - tmax (check 'temp_units')
#   - tmin (check 'temp_units') 
#   [- swrad [langleys for F temp_units] ]-currently unavailable at Chim
#   - ymdhms_v

in_datafil = settings_test.climate_data_file
Data = pd.read_csv(in_datafil, skiprows=5, delim_whitespace=True, header=None)
Data_columns = ['year', 'month', 'day', 'hour', 'min', 'sec', 'tmax', 'tmin', 'precip']

# daily data
#tmax must be F
#tmin must be F
#precip must be inch

if Data.shape[1] == 5:
    Data_columns += ['swrad']

Data.columns = Data_columns

if 'swrad' in Data.columns:
    Data['swrad'] = Data['swrad'] * 1.e+6 / 41480.
    # MJ/m2 -> Langeley (1 Langely = 41840 J/m2)

## ------------------------------------------------------------------------
# Generally, do no change below here

# -- These files will be generated (names must match those in Control file!)
#    (generally don't change this)
empty_datafil = PRMSinput_dir + slashstr + 'empty.day'
precip_datafil = PRMSinput_dir + slashstr + 'precip.day'
tmax_datafil = PRMSinput_dir + slashstr + 'tmax.day'
tmin_datafil = PRMSinput_dir + slashstr + 'tmin.day'
hum_datafil = PRMSinput_dir + slashstr + 'humidity.day' # needed for Penman-Monteith
solrad_datafil = PRMSinput_dir + slashstr + 'swrad.day'

# - how many of above met variables to process (0 thru N)
N = 3

# - Write to data variable files
print_fmt1 = ["%4d"] + ["%2d"] * 5
print_fmt2 = ["%4d"] + ["%2d"] * 5 + ["%6.2f"] * nhru

 
try:
    swrad
except:
    pass
else:
    N = 5

for ii in range(0, N+1):
    if ii==0:
        outdatafil = empty_datafil
        data = []
        label = ['precip 0', 'tmax 0', 'tmin 0']
    elif ii==1:
        outdatafil = precip_datafil
        data = Data['precip']
        label = ['precip {}'.format(nhru)]
    elif ii==2: 
        outdatafil = tmax_datafil
        data = Data['tmax']
        label = ['tmaxf {}'.format(nhru)]
    elif ii==3: 
        outdatafil = tmin_datafil
        data = Data['tmin']
        label = ['tminf {}'.format(nhru)]
    elif ii==4: 
        outdatafil = hum_datafil
        data = Data['tmax']/100.
        label = ['humidity_hru {}'.format(nhru)]
    elif ii==5: 
        outdatafil = solrad_datafil
        data = Data['swrad']
        label = ['swrad {}'.format(nhru)]

    fid = open(outdatafil, "w+")

    header = "\n".join([descr_str] + label) + "\n##########\n"
    fid.write(header)

    if len(data) > 0:
        data0 = pd.concat([Data[['year', 'month', 'day', 'hour', 'min', 'sec']],
                          pd.DataFrame(np.dstack([data] * nhru)[0])], axis=1)
        np.savetxt(fid, data0, fmt=print_fmt2, delimiter=" ")

    else:
        data0 = Data[['year', 'month', 'day', 'hour', 'min', 'sec']]
        np.savetxt(fid, data0, fmt=print_fmt1, delimiter=" ")

    fid.close()

