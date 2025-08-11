### Save Data to HDF5

from cmath import isnan
import os
from re import S
import h5py as hdf
import numpy as np

# from ..Tools.Utilities import tools
# from ..Parameters import shared_parameters


class HDF5data():
    def __init__(self,filename,t_path='sim_time',z_path='z',var_path=None,
        iteration=-1,relative_path=True,**kwargs):
        
        self.iteration  = iteration
        self.filename   = filename

        self.groupnames = []
        ### Right now z and time are the same for each iteration
        self.t_path          = f'{self.group}sim_time'  if relative_path else t_path
        self.z_path          = f'{self.group}z'         if relative_path else z_path

        self.revative_path   = relative_path
        self._var_path        = f'{self.group}{var_path}'  if relative_path else var_path

    @property
    def group(self):
        return f'/iter_{str(self.iteration)}/'   

    @property
    def var_path(self):
        return  self._var_path
    
    @var_path.setter
    def var_path(self,value):
        self._var_path = f'{self.group}{value}'  if self.relative_path else value


    def loadMeasured(self,**kwargs):
        '''
        Load the data. Itation -1 is the default where no data has been collected
        -- in this case, return nothing
        '''
        
        if self.iteration<0:
            return (None,None,None)


        with hdf.File(self.filename,'r') as file:
            # file.visit(lambda x : print(x))
            t       = file[self.t_path][...]
            z       = file[self.z_path][...]
            Data    = file[self.var_path][...] #if self.var_path in file.keys() else None

        return (t,z,Data)

    

    def save_data(self,filename=None,advance=True,**kwargs):
        '''
        Add the data to the file of interest. 
        -- kwargs is a dictionary of what you want to include in the file 
        '''
        if advance:
            self.iteration += 1;
        
        if filename is None: filename = self.filename
        with hdf.File(filename,'a') as file:
            file.require_group(self.group)
            for key in kwargs.keys():
                file[self.group][key] = kwargs[key]
        
    
    


'''
testData = np.ones([5,5])

test.new_SimData('Test1','RhoData',testData)
test.get_groups()
print (test.pullMatrix('Test1','RhoData'))
#test.loadMeasured(test.dataPath,dataset='20201211')
index = (test.measured_index(dataset='20201211',index = 500))
print (index)'''