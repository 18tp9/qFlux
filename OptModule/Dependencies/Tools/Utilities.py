### file for indexing/error calcs functions/plotting etc


import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
# from scipy import array,arange,exp
from matplotlib.widgets import Slider
import matplotlib.pyplot as pl
import h5py
# from Graphics.VideoGen import generate_video
# from Graphics.VideoGen import make_gif
import os
from ...Parameters import shared_parameters
# from .Dusolt_data import *

class tools():
    def __init__(self) -> None:
         self.precision = 1
         self.sigdigs   = 2
    
    @staticmethod
    def data_index(data,index):
            '''
            index into measured/lab data
            '''
            #data = group[dataset][...]
            #self.curr_rho = data[:,index]
            return data[:,index]


    # @staticmethod 
    # def get_hData(dataset,t,func,precision = shared_parameters.dt0, sigdigs = 2):
    #         '''
    #         t is index -- to be fixed
    #         '''
    #         #index = np.argmax(self.time==t)
    #         ### get data at time t
    #         #data = self.data_index(int(t/self.precision))


    #         '''
    #         Jason --- There is a bug here, It t= 0.001 , the try statement could be fine but the data will be from some other time  
    #         '''
    #         # try:
    #         #     ### get data at time t
    #         if (t/precision) % 1 == 0:
    #             data = func(dataset,int(t/precision))

    #         # except:
    #         else:
    #             ###if no data point at time t, round
    #             roundedT    = round(t,sigdigs)
    #             index1      = int(roundedT/precision)

    #             ###make sure we rounded down
    #             if roundedT > t:
    #                 index1 = index1 - 1
    #             ###add 1 to get second index and average the two values
    #             index2 = index1 + 1
    #             if index1 < 10**(sigdigs):          #make sure index 2 is not out of bounds
    #                 data   = (func(dataset,index1) + func(dataset,index2))/2.0
    #             else:
    #                 data   = (func(dataset,(index1 - 1)))
    #         return data

    @staticmethod
    def sortProfile(rho):
        
        rho[::-1] = np.sort(rho)

        return rho

    @staticmethod
    def error_calc(rho,rho_lab,z):
            '''
            *** TO FILL 
            '''
            frho = (rho-rho_lab)**2
            err1 = np.trapz(frho,z)
            return err1

    @staticmethod
    def error_pointwise(rho,rho_lab,n):
            '''
            *** TO FILL 
            '''
            frho = sum((rho-rho_lab)**2)
            frho = frho/n
            return frho
    
    @staticmethod 
    def truncate(data,lowerB,upperB):
        '''
        selects time window for labData
        '''
        # newData = []
        # for i in range(upperB-lowerB+1):
        #     index = lowerB + i
        #     newData.append(data[:,index])
        
        if data.ndim==2:
            newData = data[:,lowerB:upperB+1]
            return newData
        elif data.ndim==1:
            newData = data[lowerB:upperB+1]
            return newData
        else:
            raise IndexError

        # return np.asarray(newData).T

    @staticmethod
    def trunc_vert(data,num_removed):
        '''
        removes points from lowest depth of labData
        '''
        newData = []

        if data.ndim == 1:
            return data[(num_removed-1):]
        
        else:

            for i in range(data.shape[0]-num_removed):
                index = (data.shape[0] -1) - i
                newData.append(data[index,:])
            
            return np.asarray(newData).T

    @staticmethod
    def grid_change(z,data,newGrid):
        ###takes in data set and interpolates values in space to match given grid
        #newGrid = np.repeat(newGrid.z,int(self.Data.shape[1]),axis=0)
        interpolate = interp1d(z,data,axis=0,kind='cubic',fill_value='extrapolate')

        newRho = interpolate(newGrid)
        # newRho = np.flipud(newRho)

        '''for i in range(data.shape[1]):
            spatial_data = data[:,i]


        newRho = griddata(z,data,newGrid,method='linear')'''
        return newRho
    
    @staticmethod
    def subsample(z,data,newZ):
            
            raw_location = [np.argmin(np.abs(z - value)) for value in newZ]

            newData = []
            for i in range(len(newZ)):
                zL0 = newZ[i]
                z_close = z[raw_location[i]]
                
                if (zL0 < z_close) and ((i+1) < (len(raw_location))): #closest point is less than labZ

                    zp1 = z[raw_location[i+1]]  #opposite side of closest point

                    #weights for interpolation
                    w1 = abs(z_close - zL0)/abs(zL0-zp1)
                    w2 = abs(z_close - zp1)/abs(zL0-zp1)

                    interpolated = (data[raw_location[i]]*w2 + data[raw_location[i+1]]*w1)
                    newData.append(interpolated)


                elif (zL0 > z_close):     # closest point is less than labZ

                    zm1 = z[raw_location[i-1]]

                    w1 = abs(z_close - zL0)/abs(zL0-zm1)
                    w2 = abs(z_close - zm1)/abs(zL0-zm1)

                    interpolated = (data[raw_location[i]]*w1 + data[raw_location[i-1]]*w2)
                    newData.append(interpolated)

                
                else:   #just append closest point
                    newData.append(data[raw_location[i]])

            return np.asarray(newData)

        # newData = np.zeros(newZ.shape)

        # ### Weighted interpolation
        # for i,zL0 in enumerate(newZ):

        #     arg_upper = np.argmax(z>zL0) if np.argmax(z>zL0)>0 else 0
        #     arg = np.argmax(z>zL0)

 
        #     if arg_upper >0:
        #         z1        = z[arg_upper] 
        #         z0        = z[arg_upper-1]
        #         w1 = abs(z1 - zL0)/abs(z0-z1) 
        #         w0 = abs(zL0 - z0)/abs(z0-z1)
        #         newData[i] = data[arg_upper-1]*w1 + data[arg_upper]*w0
        #     else:
        #         newData[i] = data[0]
            

    @staticmethod
    def normalizeMax(temps):
        max1 = temps[0,0]
        min1 = temps[-1,0]

        for i in range(temps.shape[1]):
             
             maxVal = temps[0,i]
             minVal = temps[-1,i]

             temps[:,i] = min1*((temps[:,i] - minVal)/(maxVal-minVal)) +(max1 - min1)

        return temps
        

    @staticmethod
    def revert_grid(z,data,newGrid,extrap_func):
        interpolate = interp1d(z,data,axis=0,kind='cubic',fill_value='extrapolate')
        extrapolate = extrap_func(interpolate)
        newRho = extrapolate()
        newRho = np.flipud(newRho)
        return newRho

    @staticmethod
    def neg_filter(data):
        for i in range(data.shape[0]):
            for k in range(data.shape[1]):
                if data[i,k] < 0:
                    data[i,k] = 0
        return data
    
    # @staticmethod
    # def get_dusoltData(outpath,grid):

    #     hypso = spline(grid)


    #     with hdf5.File(outpath,'w') as h5:

    #     ### Putting time in terms of days to not worry about datetime
    #         save_data = {'t': range(1,365),
    #                         'z':depths.astype('float64'),
    #                         'A':hypso.astype('float64'),
    #                         'T':temps.astype('float64')
    #                         }
    #         for key in save_data.keys():
    #             h5[key] = save_data[key]

    @staticmethod
    def __Q_LH_calc(windspeed,surfaceTemp,airTemp,relHum):
         ### Latent heat flux calc
         Le     = 2.453e6   #J/kg
         Cl     = 1.3e-3    
         p_a    = 1.2       #kg/m3
         U_ten  = windspeed
         T_s    = surfaceTemp
         T_a    = airTemp
         R_h    = relHum
         P      = 1013.25   #hpa

         ea = (R_h/100)*np.exp(2.303*((7.5*T_a)/(T_a + 237.3) + 0.7858))
         es = np.exp(2.3026*((7.5*T_s)/(T_s + 237.3) + 0.7858))

         Q_LH = (0.622/P)*(Le*Cl*p_a*U_ten)(ea - es)
         
         return Q_LH
    
    @staticmethod
    def __Q_SH_calc(windspeed,surfaceTemp,airTemp):
         ### Sensible heat flux calc
         C_s    = 1.3e-3
         p_a    = 1.2       #kg/m3
         C_p    = 1003      #J/kg*C
         U_ten  = windspeed
         T_s    = surfaceTemp
         T_a    = airTemp

         Q_SH   = C_s*p_a*C_p*U_ten*(T_a - T_s)

         return Q_SH
    
    @classmethod
    def surfaceFlux(self,shortwave,longwave,windspeed,surfaceTemp,airTemp,relHum):
         '''
         returns sum of surface heat fluxes
         '''
         SH = self.__Q_SH_calc(windspeed,surfaceTemp,airTemp)
         LH = self.__Q_LH_calc(windspeed,surfaceTemp,airTemp,relHum)

         H = shortwave + longwave + SH + LH

         return H
    
         
        
    @staticmethod
    def adjTimeConv(timeVector):
        return np.fliplr(timeVector)
        
            
        




    