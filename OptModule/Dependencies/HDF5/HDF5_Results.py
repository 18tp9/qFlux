### Save Data to HDF5
### HDF5 class for processing and plotting post-simulation data

from cmath import isnan
import os
import h5py as hdf
import numpy as np

from ..Tools.Utilities import tools
from ...Parameters import shared_parameters


class HDF5data():
    def __init__(self,filename,
        t_path=(r'/LabData/t'),
        z_path= (r'/LabData/z'),
        var_path=(r'rho'),
        **kwargs):
        
        self.filename   = filename
        self.var_path   = var_path
        self.t_path     = t_path
        self.z_path     = z_path
        self.tArrPath   = (r'/Timeseries/scales/time')
        self.zArrPath   = (r'/Timeseries/scales/z')
        self.RawZ_path  = (r'/Timeseries/scales/z_raw')
        self.RawT_path  = (r'/Timeseries/Temperature/T_raw')
        self.groupnames = []
        self.precision  = 1e-2
        self.sigdigs    = 2
        self.finalIter  = (('iter_' + str(shared_parameters.num_iter)))

        
    def loadMeasured(self,**kwargs):
        with hdf.File(self.filename,'r') as file:
            t       = file[self.t_path][...]
            z       = file[self.z_path][...]
            Data    = file[self.finalIter][self.var_path][...]

        return self.__clean_data(t,z,Data)


    def __clean_data(self,t_lab,z_lab,Data):
        '''
        Truncate and filter
        '''
        params       = shared_parameters()

        t_clean      = tools.truncate(t_lab,shared_parameters.trunc_lowerB,shared_parameters.trunc_upperB)
        z_clean      = params.grid.z
        Data   = tools.truncate(Data,shared_parameters.trunc_lowerB,shared_parameters.trunc_upperB)
        labData      = tools.grid_change(z_lab,Data,z_clean)
        

        ### Refine 
        t_clean = (t_clean - t_clean[0])
        ### if the lab data is more fine than the simulation, dump data  
        dt_lab = np.amin(np.diff(t_clean))
        if (dt_lab/shared_parameters.dt0)>1:
            t_clean[::np.ceil(dt_lab/shared_parameters.dt0)]
            labData[:,::np.ceil(dt_lab/shared_parameters.dt0)]


        return (t_clean,z_clean,labData)


    def save_data(self,filename,**kwargs):
        '''
        Add the data to the file of interest. 
        -- kwargs is a dictionary of what you want to include in the file 
        '''
        with hdf.File(filename,'a') as file:
            file.require_group('LabData')
            for key in kwargs.keys():
                if not key in file['LabData']:
                    file['LabData'][key] = kwargs[key]
        
    
        


    
    


'''
testData = np.ones([5,5])

test.new_SimData('Test1','RhoData',testData)
test.get_groups()
print (test.pullMatrix('Test1','RhoData'))
#test.loadMeasured(test.dataPath,dataset='20201211')
index = (test.measured_index(dataset='20201211',index = 500))
print (index)'''