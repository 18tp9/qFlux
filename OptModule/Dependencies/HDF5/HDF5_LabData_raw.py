### Save Data to HDF5

from cmath import isnan
import os
import h5py as hdf
import numpy as np

from ..Tools.Utilities import tools
from ...Parameters import shared_parameters


class HDF5data():
    def __init__(self,filename ):
        
        self.filename   = filename
        self.tArrPath   = (r'/Timeseries/scales/time')
        self.zArrPath  = (r'/Timeseries/scales/z_raw')
        self.dataPath  = (r'/Timeseries/Temperature/T_raw')
        self.groupnames = []
        self.precision  = 1e-2
        self.sigdigs    = 2

        self.params = shared_parameters()
        
    def loadMeasured(self,str=None):
        with hdf.File(self.filename,'r') as file:
            t       = file[self.tArrPath][...]
            z       = file[self.zArrPath][...]
            Data    = file[self.dataPath][...]

        return self.__clean_data(t,z,Data)

    def __clean_data(self,t,z,Data):
        '''
        Truncate and filter
        '''
        t_clean      = tools.truncate(t,shared_parameters.trunc_lowerB,shared_parameters.trunc_upperB)
        z_clean      = z
        data_clean   = tools.truncate(Data,shared_parameters.trunc_lowerB,shared_parameters.trunc_upperB)
        return (t_clean,z_clean,data_clean)



    def save_data(self,filename,**kwargs):
        '''
        Add the data to the file of interest. 
        -- kwargs is a dictionary of what you want to include in the file 
        '''
        with hdf.File(filename,'a') as file:
            file.create_group('LabData')
            for key in kwargs.keys():
                file['LabData'][key] = kwargs[key]
        
    
        


    
    


'''
testData = np.ones([5,5])

test.new_SimData('Test1','RhoData',testData)
test.get_groups()
print (test.pullMatrix('Test1','RhoData'))
#test.loadMeasured(test.dataPath,dataset='20201211')
index = (test.measured_index(dataset='20201211',index = 500))
print (index)'''