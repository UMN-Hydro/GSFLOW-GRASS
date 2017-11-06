import numpy as np
from matplotlib import pyplot as plt
import sys
# GRASS
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import miscellaneous as m
from grass.pygrass.gis import region
from grass.pygrass import vector # Change to "v"?
from grass.script import vector_db_select
from grass.pygrass.vector import Vector, VectorTopo
from grass.pygrass.raster import RasterRow
from grass.pygrass import utils
from grass import script as gscript

def get_columns_in_order(vect, cols, nodata_value=-999):
  colNames = np.array(gscript.vector_db_select(vect, layer=1)['columns'])
  colValues = np.array(gscript.vector_db_select(vect, layer=1)['values'].values())
  outlist = []
  for col in cols:
    newcol = colValues[:,colNames == col].squeeze()
    # If column does not exist, populate with nodata value
    # Strangely, has shape (nrows, 0)
    if np.prod(newcol.shape) == 0:
      newcol = (-999*np.ones(colValues.shape[0])).astype(int).astype(str)
    outlist.append(newcol)
  return outlist
  
# Reaches
columns_in_order = ['KRCH', 'IRCH', 'JRCH', 'ISEG', 'IREACH', 'RCHLEN', 'STRTOP', 'SLOPE', 'STRTHICK', 'STRHC1', 'THTS', 'THTI', 'EPS', 'UHC']
outcols = get_columns_in_order('reaches', columns_in_order)
outarray = np.array(outcols).transpose()
outtable = np.vstack((columns_in_order, outarray))
np.savetxt('reach_data.txt', outtable, fmt='%s', delimiter=',')

# Segments
columns_in_order = ['NSEG', 'ICALC', 'OUTSEG', 'IUPSEG', 'IPRIOR', 'NSTRPTS', 'FLOW', 'RUNOFF', 'ETSW', 'PPTSW', 'ROUGHCH', 'ROUGHBK', 'CDPTH', 'FDPTH', 'AWDTH', 'BWDTH']
outcols = get_columns_in_order('segments', columns_in_order)
outarray = np.array(outcols).transpose()
outtable = np.vstack((columns_in_order, outarray))
np.savetxt('segment_data_4A_INFORMATION.txt', outtable, fmt='%s', delimiter=',')

columns_in_order = ['HCOND1', 'THICKM1', 'ELEVUP', 'WIDTH1', 'DEPTH1', 'THTS1', 'THTI1', 'EPS1', 'UHC1']
outcols = get_columns_in_order('segments', columns_in_order)
outarray = np.array(outcols).transpose()
outtable = np.vstack((columns_in_order, outarray))
np.savetxt('segment_data_4B_UPSTREAM.txt', outtable, fmt='%s', delimiter=',')

columns_in_order = ['HCOND2', 'THICKM2', 'ELEVDN', 'WIDTH2', 'DEPTH2', 'THTS2', 'THTI2', 'EPS2', 'UHC2']
outcols = get_columns_in_order('segments', columns_in_order)
outarray = np.array(outcols).transpose()
outtable = np.vstack((columns_in_order, outarray))
np.savetxt('segment_data_4C_DOWNSTREAM.txt', outtable, fmt='%s', delimiter=',')

