import os
import numpy as np
import time
from build_ini import BuildINI

bi = BuildINI()
bi.DEM_input = ''
bi.proj_name = 'Shullcas'
bi.start_date = '2013-08-26'
bi.end_date = '2016-09-29'
bi.init_start_date = '2013-08-26'
bi.climate_data_file = '/home/awickert/GSFLOW/Shullcas/Huaytapallana_MetData_metric.txt'
#bi.threshold_drainage_area_meters2 = '9000000'
#bi.MODFLOW_grid_resolution_meters = '500'
bi.outlet_point_x = '482452'
bi.outlet_point_y = '8672978'
bi.NLAY = '1'
bi.DZ = '200'
bi.fl_create_hydcond = '0'
#bi.hydcond = 'hydcond_test.txt'
bi.finf = '0.0015'

inifile = os.getcwd() + '/test.ini'

drainage_thresholds = np.linspace(1E6, 17E6, 5)
MODFLOW_grid_sizes = np.linspace(100, 900, 5)

#drainage_thresholds = np.linspace(9E6, 10E6, 1)
#MODFLOW_grid_sizes = np.linspace(500, 900, 1)

outlist = []
for da_thresh in drainage_thresholds:
    for modflow_gs in MODFLOW_grid_sizes:
        bi.threshold_drainage_area_meters2 = str(da_thresh)
        bi.MODFLOW_grid_resolution_meters = str(modflow_gs)
        bi.proj_name = 'gridTestShullcas_' + 'D' + str(da_thresh) + '_M' + str(modflow_gs)
        bi.writeINI(inifile)
        
        print modflow_gs, da_thresh
        
        t1 = time.time()
        os.system('python /home/awickert/models/GSFLOW-GRASS/domain_builder/buildDomainGRASS.py ' + inifile)
        t2 = time.time()
        dt_domain = t2-t1
        
        t1 = time.time()
        os.system('sh /home/awickert/models/GSFLOW-GRASS/Run/goGSFLOW.sh ' + inifile)
        t2 = time.time()
        dt_input_run = t2-t1
        
        outlist.append([bi.MODFLOW_grid_resolution_meters, bi.threshold_drainage_area_meters2, dt_domain, dt_input_run, bi.proj_name])
        
np.savetxt('outlist_'+bi.proj_name.split('_')[0]+'.txt', outlist, '%s')


# Organize output
#outlist = np.genfromtxt('outlist.txt', dtype=str)
outlist = np.array(outlist)
runnames = outlist[:,-1]
params = outlist[:,:2].astype(float)
runtimes = outlist[:,2:4].astype(float)

from matplotlib import pyplot as plt

MODFLOW_grid_sizes__mod = np.hstack(( MODFLOW_grid_sizes[:4], MODFLOW_grid_sizes[-1] ))

# Plot time to create array
t_creation = np.flipud(runtimes[:,0].reshape(5,5)) # drainage area x, modflow y
t_creation__mod = np.flipud(np.hstack(( t_creation[:,:4], np.array([t_creation[:,-1]]).transpose() )))
plt.figure()
#plt.imshow(t_creation)
#plt.colorbar()
ax = plt.subplot()
plt.contourf(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, t_creation__mod, [100, 150, 200, 250, 300, 350, 400], alpha=0.5)
cs = plt.contour(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, t_creation__mod, [125, 150, 200, 250, 300, 350], colors='k')
ax.set_xscale('log')
ax.set_yscale('log')
plt.clabel(cs, fmt='%d s', inline_spacing=5)
plt.xlabel('MODFLOW grid cell area [km$^2$]', fontsize=16)
plt.ylabel('Drainage threshold [km$^2$]', fontsize=16)
plt.tight_layout()
plt.savefig('t_creation_'+bi.proj_name.split('_')[0]+'.svg')
plt.show()

# Plot time to run array
t_run = np.flipud(runtimes[:,1].reshape(5,5)) # drainage area x, modflow y
t_run__mod = np.flipud(np.hstack(( t_run[:,:4], np.array([t_run[:,-1]]).transpose() )))
plt.figure()
#plt.imshow(t_run)
#plt.colorbar()
ax = plt.subplot()
plt.contourf(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, t_run__mod, [0, 50, 100, 150, 200, 250, 300], alpha=0.5)
cs = plt.contour(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, t_run__mod, [10, 50, 100, 150, 200, 250], colors='k')
ax.set_xscale('log')
ax.set_yscale('log')
plt.clabel(cs, fmt='%d s', inline_spacing=5)
plt.xlabel('MODFLOW grid cell area [km$^2$]', fontsize=16)
plt.ylabel('Drainage threshold [km$^2$]', fontsize=16)
plt.tight_layout()
plt.savefig('t_run_'+bi.proj_name.split('_')[0]+'.svg')
plt.show()

# Load all
Qarray = []
for runname in runnames:
    try:
        Q = np.genfromtxt('/home/awickert/GSFLOW2018/'+runname+'/outputs/PRMS_GSFLOW/gsflow.csv', delimiter=',',skip_header=1)[:,6]
    except:
        Q = np.nan*np.zeros(2004)
    Qarray.append(Q)

Qarray = np.array(Qarray)

# mm/day to m3/s
Qarray_m3s = Qarray / (24. * 3600.)

rmse = []
for Qrow in Qarray_m3s:
    _rmse = np.mean( ( (Qarray_m3s[0] - Qrow)**2 )**.5 )
    rmse.append(_rmse)

rmse = np.array(rmse).reshape(5, 5)
rmse = np.flipud(rmse)

ddt = np.mean(np.diff(drainage_thresholds))/2.
ddm = np.mean(np.diff(MODFLOW_grid_sizes))/2.

MODFLOW_grid_sizes__mod = np.hstack(( MODFLOW_grid_sizes[:4], MODFLOW_grid_sizes[-1] ))
rmse__mod = np.flipud(np.hstack(( rmse[:,:4], np.array([rmse[:,-1]]).transpose() )))

"""
plt.figure()
plt.imshow(rmse/np.mean(Qarray_m3s[0])*100., aspect=10000**2./100E6, extent=[((MODFLOW_grid_sizes[0]-ddm)/1000.)**2., ((MODFLOW_grid_sizes[-1]+ddm)/1000.)**2, (drainage_thresholds[0]-ddt)/1000.**2, (drainage_thresholds[-1]+ddt)/1000.**2])
plt.colorbar()
cs = plt.contour(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, rmse__mod/np.mean(Qarray_m3s[0])*100., [10, 15, 20, 25, 30, 35], colors='k')
plt.clabel(cs, fmt='%r %%')
plt.xlabel('MODFLOW grid cell area [km$^2$]', fontsize=16)
plt.ylabel('HRU area [km$^2$]', fontsize=16)
"""

plt.figure()
ax = plt.subplot()
plt.contourf(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, rmse__mod/np.mean(Qarray_m3s[0])*100., np.linspace(0, 130, 14), cmap='plasma', alpha=.5)
cs = plt.contour(MODFLOW_grid_sizes__mod**2/1E6, drainage_thresholds/1E6, rmse__mod/np.mean(Qarray_m3s[0])*100., [5, 10, 20, 30, 40, 60, 80, 100], colors='k')
ax.set_xscale('log')
ax.set_yscale('log')
plt.clabel(cs, fmt='%d%%')
plt.xlabel('MODFLOW grid cell area [km$^2$]', fontsize=16)
plt.ylabel('Drainage threshold [km$^2$]', fontsize=16)
plt.tight_layout()
plt.savefig('rmse_'+bi.proj_name.split('_')[0]+'.svg')
plt.show()

