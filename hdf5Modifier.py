import numpy as np 
import matplotlib.pyplot as pl
import h5py as hdf5
import random
import os

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data


simDataName = os.path.join('DDC_data','LabData_20201211.h5')
outputFile = os.path.join('DDC_data','IceOn_supershort_20201211.h5')

h5Data = HDF5data(simDataName,            
        t_path=('Timeseries/scales' + '/time'),
        z_path= ('Timeseries/scales' + '/z'),
        var_path=('Timeseries/Temperature' + '/T_Chain'),
        # format='LabData',
        )

(t_lab,z_lab,rho_lab) = h5Data.loadData()

iceOnTime = 22350   #timeIndex

t_lab = t_lab[iceOnTime:iceOnTime+1500]
rho_lab = rho_lab[:,iceOnTime:iceOnTime+1500]

with hdf5.File(outputFile,'a') as h5:
    
    ## Base Case 
    # save_data = {'t': np.array(t_Array)[::(len(t_Array)//512)],
    #         'z':z,
    #         'rho':np.array(rho_Array).T[:,::(len(t_Array)//512)],}
    save_data = {'t': t_lab,
            'z':z_lab,
            'rho':rho_lab
            }

    for key in save_data.keys():
        h5[key] = save_data[key]