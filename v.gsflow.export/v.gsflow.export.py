#!/usr/bin/env python
############################################################################
#
# MODULE:       v.gsflow.export
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Export database tables and pour point for GSFLOW input and control files
#
# COPYRIGHT:    (c) 2016-2017 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  uses inputs from the v.gsflow series
 
# More information
# Started December 2016

#%module
#% description: Export databse tables and pour point for GSFLOW input and control files
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: GSFLOW
#%end

#%option G_OPT_V_INPUT
#%  key: hru_input
#%  label: Sub-basin hydrologic response units
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_INPUT
#%  key: grid_input
#%  label: MODFLOW grid
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_INPUT
#%  key: gravres_input
#%  label: Union of MODFLOW grid and HRUs
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_INPUT
#%  key: pour_point_input
#%  label: Pour point at the outlet of the basin
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_OUTPUT
#%  key: hru_output
#%  label: PRMS-style HRU output table for GSFLOW input, no file ext
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_OUTPUT
#%  key: grid_output
#%  label: MODFLOW output table for GSFLOW input, no file ext
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_OUTPUT
#%  key: gravres_output
#%  label: Gravity Reservoir output table for GSFLOW input, no file ext
#%  required: yes
#%  guidependency: layer,column
#%end

#%option G_OPT_V_OUTPUT
#%  key: pour_point_output
#%  label: Pour point coordinates for GSFLOW input, no file ext
#%  required: yes
#%  guidependency: layer,column
#%end

##################
# IMPORT MODULES #
##################

# PYTHON
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

###############
# MAIN MODULE #
###############

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

def main():
    """
    Build gravity reservoirs in GSFLOW: combines MODFLOW grid and HRU sub-basins
    These define the PRMS soil zone that connects to MODFLOW cells
    """

    ##################
    # OPTION PARSING #
    ##################

    # Parser
    options, flags = gscript.parser()

    # Input
    reaches = options['reaches_input']
    segments = options['segments_input']
    gravity_reservoirs = options['gravres_input']
    pour_point = options['pour_point_input']
    
    # Output
    out_reaches = options['reaches_output']
    out_segments = options['segments_output']
    out_gravity_reservoirs = options['gravres_output']
    out_pour_point = options['pour_point_output']
    
    ##############
    # PROCESSING #
    ##############

    # Reaches
    ##########
    columns_in_order = ['KRCH', 'IRCH', 'JRCH', 'ISEG', 'IREACH', 'RCHLEN', 'STRTOP', 'SLOPE', 'STRTHICK', 'STRHC1', 'THTS', 'THTI', 'EPS', 'UHC']
    outcols = get_columns_in_order(reaches, columns_in_order)
    outarray = np.array(outcols).transpose()
    outtable = np.vstack((columns_in_order, outarray))
    np.savetxt(out_reaches+'.txt', outtable, fmt='%s', delimiter=',')

    # Segments
    ###########
    columns_in_order = ['NSEG', 'ICALC', 'OUTSEG', 'IUPSEG', 'IPRIOR', 'NSTRPTS', 'FLOW', 'RUNOFF', 'ETSW', 'PPTSW', 'ROUGHCH', 'ROUGHBK', 'CDPTH', 'FDPTH', 'AWDTH', 'BWDTH']
    outcols = get_columns_in_order(segments, columns_in_order)
    outarray = np.array(outcols).transpose()
    outtable = np.vstack((columns_in_order, outarray))
    np.savetxt(out_segments+'_A_INFORMATION.txt', outtable, fmt='%s', delimiter=',')

    columns_in_order = ['HCOND1', 'THICKM1', 'ELEVUP', 'WIDTH1', 'DEPTH1', 'THTS1', 'THTI1', 'EPS1', 'UHC1']
    outcols = get_columns_in_order(segments, columns_in_order)
    outarray = np.array(outcols).transpose()
    outtable = np.vstack((columns_in_order, outarray))
    np.savetxt(out_segments+'_segment_data_4B_UPSTREAM.txt', outtable, fmt='%s', delimiter=',')

    columns_in_order = ['HCOND2', 'THICKM2', 'ELEVDN', 'WIDTH2', 'DEPTH2', 'THTS2', 'THTI2', 'EPS2', 'UHC2']
    outcols = get_columns_in_order(segments, columns_in_order)
    outarray = np.array(outcols).transpose()
    outtable = np.vstack((columns_in_order, outarray))
    np.savetxt(out_segments+'_4C_DOWNSTREAM.txt', outtable, fmt='%s', delimiter=',')

    # Gravity reservoirs
    #####################
    columns_in_order = ['gvr_hru_id', 'gvr_hru_pct', 'gvr_cell_id', 'gvr_cell_pct']
    outcols = get_columns_in_order(gravity_reservoirs, columns_in_order)
    outarray = np.array(outcols).transpose()
    outtable = np.vstack((columns_in_order, outarray))
    np.savetxt(out_gravity_reservoirs+'.txt', outtable, fmt='%s', delimiter=',')

    # Pour Point
    #####################
    _y, _x = np.squeeze(gscript.db_select(sql='SELECT row,col FROM '+pour_point))
    outstr = 'discharge_pt: row_i '+_y+' col_i '+_x
    outfile = file(out_pour_point+'.txt', 'w')
    outfile.write(outstr)
    outfile.close()
    
if __name__ == "__main__":
    main()

