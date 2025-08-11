### reading dusolt excel data

import numpy as np
import matplotlib.pyplot as pl
import plotly.express as px
import pandas as pd
import os
from datetime import datetime
# from OptModule.Dependencies.EquationSolvers.DiffusionSolver_base.Operators import CHEB
import scipy.interpolate as interp
import h5py as hdf5
import scipy.signal as sig


class Dusolt(object):
    z_grid = [1,5,9,13,19,25,31,36]
    z_grid = np.asarray(z_grid).T
    

    wq1     = pd.read_csv('Dusolt\LakeModelData\WQ_1.csv')
    times   = wq1['time']
    wq1T    = np.asarray(wq1['temp'])

    wq5     = pd.read_csv('Dusolt\LakeModelData\WQ_5.csv')
    wq5T    = np.asarray(wq5['temp'])

    wq9     = pd.read_csv('Dusolt\LakeModelData\WQ_9.csv')
    wq9T    = np.asarray(wq9['temp'])

    wq13    = pd.read_csv('Dusolt\LakeModelData\WQ_13.csv')
    wq13T   = np.asarray(wq13['temp'])

    wq19    = pd.read_csv('Dusolt\LakeModelData\WQ_19.csv')
    wq19T   = np.asarray(wq19['temp'])

    wq25    = pd.read_csv('Dusolt\LakeModelData\WQ_25.csv')
    wq25T   = np.asarray(wq25['temp'])

    wq31    = pd.read_csv('Dusolt\LakeModelData\WQ_31.csv')
    wq31T   = np.asarray(wq31['temp'])

    wq36    = pd.read_csv('Dusolt\LakeModelData\WQ_36.csv')
    wq36T   = np.asarray(wq36['temp'])

    startTime = 0
    endTime = 365

    Dusolttime = []
    for timestamp in times:
        t = datetime.strptime(timestamp,'%Y-%m-%d %H:%M:%S')
        Dusolttime.append(t.date())

    Dusolttime = np.asarray(Dusolttime)

    # temps = np.zeros([9,365])
    temps = np.matrix([wq1T,wq5T,wq9T,wq13T,wq19T,wq25T,wq31T,wq36T])[:,startTime:endTime] # extra row of surface temperature

    surfT = np.array(np.concatenate(temps[-1,:]))[0]#[startTime:endTime]
    # pl.contourf(time,z_grid,temps)
    # ticks = [time[len(time)//5],time[2*len(time)//5],time[3*len(time)//5],time[4*len(time)//5]]
    # pl.xticks(ticks)
    # pl.xlabel('Date')
    # pl.ylabel('z (m)')
    # pl.colorbar()

    # pl.show()

    hypsography   = pd.read_csv('Dusolt\LakeModelData\hypsography.csv')
    depths  = hypsography['depths (m)']
    areas   = hypsography['areas (m2)']
    areas = np.asarray(areas)

    # NULL,chebGrid = CHEB.cheb(len(depths)-1,xmin=0,xmax=36)

    spline = interp.interp1d(depths,areas,'cubic',fill_value='extrapolate')

    # Decay coefficient for beer's law
    nu = 0.2

    #Cloud cover 
    CC = 0.2

    ### Time in day number
    timeNum = np.linspace(startTime,endTime,endTime-startTime)

    # pl.plot(chebGrid,spline(chebGrid))
    # pl.show()



    metData = pd.read_csv('Dusolt\LakeModelData\Harpmet_new.csv')
    Sw      = metData['ShortWave'][startTime:endTime]
    Lw      = metData['LongWave'][startTime:endTime]
    T_a     = metData['AirTemp'][startTime:endTime]     ### csv has 366 data points Jan1 for both yrs
    Rh      = metData['RelHum'][startTime:endTime]
    U_ten   = metData['WindSpeed'][startTime:endTime]

    bottom  = wq36T[startTime:endTime]
    surface  = wq1T[startTime:endTime]

    @staticmethod
    def __Q_LH_calc(windspeed,surfTemp,airTemp,relHum):
         ### Latent heat flux calc
         Le     = 2.453e6   #J/kg
         Cl     = 1.3e-3    
         p_a    = 1.2       #kg/m3
         U_ten  = windspeed
         T_s    = surfTemp
         T_a    = airTemp
         R_h    = relHum
         p      = 1013.25   #hpa

         ea = (R_h/100)*np.exp(2.303*((7.5*T_a)/(T_a + 237.3) + 0.7858))
         es = np.exp(2.3026*((7.5*T_s)/(T_s + 237.3) + 0.7858))

         Q_LH = (0.622/p)*(Le*Cl*p_a*U_ten)*(ea - es)
         
         return Q_LH
    
    @staticmethod
    def __Q_SH_calc(windspeed,surfTemp,airTemp):
         ### Sensible heat flux calc
         C_s    = 1.3e-3
         p_a    = 1.2       #kg/m3
         C_p    = 1003      #J/kg*C
         U_ten  = windspeed
         T_s    = surfTemp
         T_a    = airTemp

         Q_SH   = C_s*p_a*C_p*U_ten*(T_a - T_s)

         return Q_SH
    
    @staticmethod
    def __netSW(shortwave,depth,t,nu):
         #apply beer's law to shortwave

         r_a = 0.08 + 0.02*np.sin((2*np.pi*t)/365 + np.pi/2)

         shortwave = shortwave*(1-r_a)

         SW = shortwave*(1-np.exp(-nu*depth))

         return SW
    
    @staticmethod
    def __LW(longwave,surfTemp,airTemp,CC):

        eps_a   = (9.37e-6)*((273+airTemp)**2)
        eps_w   = 0.96
        sigma   = 5.6697e-8
        r_a     = 0.03


        # Hw = - eps_w*sigma*(abs(surfTemp)+273)**4
        LW = (1-r_a)*(1+0.17*(CC**2))*eps_a*sigma*((airTemp+273)**4) 
        # LW_nett = (1-r_a)*longwave - eps_w*sigma*(abs(surfTemp)+273)**4

        # LW = (1-r_a)*LW_nett

        return LW
    @staticmethod
    def __infraRad(surfTemp):
        eps_w   = 0.96
        sigma   = 5.6697e-8

        Hw = - eps_w*sigma*(abs(surfTemp)+273)**4

        return Hw
    
    @classmethod
    def __surfaceFlux(self):
         '''
         returns sum of surface heat fluxes
         '''
        #  SH = self.__Q_SH_calc(self.U_ten[t],surfaceTemp,self.T_a[t])
        #  LH = self.__Q_LH_calc(self.U_ten[t],surfaceTemp,self.T_a[t],self.Rh[t])
        #  ## Jason -- comment this out for SW 
        #  SW = self.__netSW(self.Sw[t],max(self.depths),t,self.nu)
        #  LW = self.__LW(self.Lw[t],surfaceTemp,self.T_a[t],self.CC)
         
         SH = self.__Q_SH_calc(self.U_ten,self.surfT,self.T_a)
         LH = self.__Q_LH_calc(self.U_ten,self.surfT,self.T_a,self.Rh)
         ## Jason -- comment this out for SW 
         SW = self.__netSW(self.Sw,max(self.depths),self.timeNum,self.nu)
         LW = self.__LW(self.Lw,self.surfT,self.T_a,self.CC)
         
         HW = self.__infraRad(self.surfT)

        #  H =  LW + SW + SH + LH
         H =  LW + SH + LH + HW # Jason

         return H
    
    @classmethod
    def getSW_I(self,nu):
        shortWave_I = []
        for i in range(self.startTime,self.endTime):
            SW = self.__netSW(self.Sw[i],self.z_grid,i,nu)
            # print(SW.shape)
            shortWave_I.append(SW)

        shortWave_I = np.asarray(shortWave_I).T

        return shortWave_I
    
    # @classmethod
    # def flux_components(self,):
    #      '''
    #      returns sum of surface heat fluxes
    #      '''
    #      SH = self.__Q_SH_calc(self.U_ten,self.surfT,self.T_a)
    #      LH = self.__Q_LH_calc(self.U_ten,self.surfT,self.T_a,self.Rh)
    #      ## Jason -- comment this out for SW 
    #      SW = self.__netSW(self.Sw,max(self.depths),self.timeNum[self.startTime:self.endTime],self.nu)
    #      LW = self.__LW(self.Lw,self.surfT,self.T_a,self.CC)
    #      HW = self.__infraRad(self.surfT)

        #  return [SW,LW,SH,LH,HW,(LW+SH+LH+HW)]

    # with hdf5.File('InputData/DusoltTemps.h5','w') as h5:

    #     ### Putting time in terms of days to not worry about datetime
    #     save_data = {'t': range(365),
    #                     'z':z_grid.astype('float64'),
    #                     'T':temps.astype('float64')#[:,150,250],
    #                     }
    #     for key in save_data.keys():
    #         h5[key] = save_data[key]
    

    @classmethod
    def heatFlux(self):
        ### returns H
        H = self.__surfaceFlux()

        H = sig.medfilt(H,27)

        H = H*(3600*24)     #converted to J/d*m^2

        return H
    
    @classmethod
    def bottTemp(self):
        # returns bottom temp for all time
        return self.bottom
    
    @classmethod
    def surfTemp(self):
        #returns surf temp for all time
        return self.surface
