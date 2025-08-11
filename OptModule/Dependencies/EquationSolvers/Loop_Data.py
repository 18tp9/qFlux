### file for indexing/error calcs functions/plotting etc


import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
# from scipy import array,arange,exp
from matplotlib.widgets import Slider
import matplotlib.pyplot as pl
import h5py
import os

from scipy.interpolate import interp1d
from ..HDF5.HDF5_Generic import HDF5data
from scipy.signal import savgol_filter


class DataManager():

    def __init__(self,filename=None,var_name=None,format="SimData",
        t_sim=None,z=None,data=None,NO_LOAD=False,
        filter_data=False,
        **kwargs):

        self.index = 0 

        ### Load the data from the filename 

        if not (filename is None):
            self.h5 = HDF5data(filename, format=format,**kwargs)
        
            self.var_name = self.h5.h5.var_path
            if not NO_LOAD:
                (t,self.z,self.data) = self.h5.loadData()
                if filter_data:
                    self.data = savgol_filter(self.data, 5, 3,axis=0,mode='mirror')
                self.t_sim = t
            
        else:
            self.set_data(t_sim,z,data,var_name)




    def update_variable(self,var_name):
        ### Update the variable of interest
        self.var_name = var_name
        self.index=0
        self.h5.h5.var_path = var_name
        (_,_,self.data) = self.h5.loadData(var_name)


    def set_data(self,t_sim,z,data,var_name=None):
        '''
        Hard code the input data 
        '''
        self.var_name = var_name
        self.t_sim,self.z = t_sim,z
        self.data = data
        self.index=0

    def reset(self):
        # Reset the index
        self.index = 0 

    def __add__(self,data2):
        '''
        Add two datasets together
        -- Currently assuming that they are the same length 
        '''
        return DataManager(data=self.data  + data2.data,t_sim = self.t_sim,z=self.z,var_name=self.var_name)


    def __mul__(self,value):
        '''
        Multiply by a float 
        '''
        if isinstance(value,float) or isinstance(value,int):
            return DataManager(data=value*self.data ,t_sim = self.t_sim,z=self.z,var_name=self.var_name)
        else:
            raise ValueError()
        
    def __rmul__(self,other):
        ## reverse
        return self.__mul__(other)


    def get_Data(self,t,reverse=False):
        '''
        Return the data at the requested time t 
        -- t is relative to start of the lab data 
        -- assumed that the units of t and lab time are identical 
         ---- reverse if the index can decrease 
        '''

        closest = np.argmin(abs(t-self.t_sim))

        return self.data[:,closest]



    def change_grid(self,data,newZ):
        ###takes in data set and interpolates values in space to match given grid

        interpolate = interp1d(self.z,data,axis=0,kind='cubic')
        z_above  = newZ[(newZ>self.z[-1])]
        z_within = newZ[(newZ<=self.z[-1])* (newZ>=self.z[0])]
        z_below  = newZ[(newZ<self.z[0] )]
        new_data = interpolate(z_within)

        if data.ndim==1:
            if len(z_below)>0: new_data = np.append(data[0]*np.ones(len(z_below)),new_data,axis=0) 
            if len(z_above)>0: new_data = np.append(new_data,data[-1]*np.ones(len(z_above)),axis=0) 
            
        elif data.ndim==2:
            if len(z_below)>0: new_data = np.append(np.tile(data[0,:],(len(z_below),1)),new_data,axis=0) 
            if len(z_above)>0: new_data = np.append(new_data,np.tile(data[-1,:],(len(z_above),1)),axis=0)
            
        else:
            raise IndexError('Note the right number of dimensions')
        
        return new_data


    def save_data(self,filename,**kwargs):
        '''
        save the data 
        '''

        self.h5.save_data(filename,**kwargs)