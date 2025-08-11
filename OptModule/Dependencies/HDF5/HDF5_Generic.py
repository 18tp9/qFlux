import os
import numpy as np 
import h5py as hdf

from .HDF5_LabData import HDF5data as HDF5data_lab
from .HDF5_LabData_raw import HDF5data as HDF5data_lab_raw
from .HDF5_SimData import HDF5data as HDF5data_sim
from .HDF5_FakeData import HDF5data as HDF5data_Fake
from .HDF5_Results import HDF5data  as HDF5data_results
# from OptModule.Tools.Utilities import tools
# from ..Parameters import shared_parameters



class HDF5data():


    def __init__(self,filename,format='LabData',**kwargs):
        '''
        One point access to the data 
        '''
        self.filename = filename
        self.format = format 

        if self.format == "LabData":
            self.h5 = HDF5data_lab(self.filename,**kwargs)
            
        elif self.format == "LabData_raw":
            self.h5 = HDF5data_lab_raw(self.filename,**kwargs)
            
        elif self.format == "SimData":
            self.h5 = HDF5data_sim(self.filename,**kwargs)

        elif self.format == 'FakeData':
            self.h5 = HDF5data_Fake(self.filename,**kwargs)

        elif self.format == 'default':
            self.h5 = HDF5data_results(self.filename,**kwargs)

        else:
            raise ValueError("Format Does not Exist")

    def loadData(self):
        '''
        Read data from a specific format 
        -- var_name is only important for the sim data 
        '''
        (t,z,Data) = self.h5.loadMeasured();
        return (t,z,Data)


        
    def save_data(self,filename,**kwargs):
        '''
        save the data to a given file
        '''
        self.h5.save_data(filename,**kwargs)
        


if __name__ == "__main__":

    import matplotlib.pyplot as pl

    filename = "..\\..\\InputData\\LabData_20201211.h5"
    h5 = HDF5data(filename, format="LabData")
    
    (Lab_t,Lab_z,LabData) = h5.loadData()
    pl.contourf(Lab_t,Lab_z,LabData)




