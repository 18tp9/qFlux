### Save Data to HDF5

from cmath import isnan
import os
import h5py as hdf
import numpy as np

from ..Tools.Utilities import tools


class HDF5data():
    def __init__(self,filename,subset,t_path='t',z_path='z',var_path=None,relative_path=True,**kwargs):
        
        self.subset     = subset
        self.filename   = filename

        self.groupnames = []
        ### Right now z and time are the same for each iteration
        self.t_path          =  t_path
        self.z_path          =  z_path

        self.relative_path   = relative_path
        self._var_path        =  var_path

    @property
    def group(self):
        return self.subset   

    @property
    def var_path(self):
        return  self._var_path
    
    @var_path.setter
    def var_path(self,value):
        self._var_path = f'{self.subset}{value}'  if self.relative_path else value


    def loadMeasured(self,**kwargs):
        '''
        Load the data. Itation -1 is the default where no data has been collected
        -- in this case, return nothing
        '''
        

        with hdf.File(self.filename,'r') as file:

            # file.visit(lambda x : print(x))
            t       = file[self.subset][self.t_path][...]
            z       = file[self.subset][self.z_path][...]
            Data    = file[self.subset][self.var_path][...] #if self.var_path in file.keys() else None

        return (t,z,Data)

    

    def save_data(self,filename=None,advance=True,**kwargs):
        '''
        Add the data to the file of interest. 
        -- kwargs is a dictionary of what you want to include in the file 
        '''

        
        if filename is None: filename = self.filename
        with hdf.File(filename,'a') as file:
            file.require_group(self.group)
            for key in kwargs.keys():
                file[self.group][key] = kwargs[key]