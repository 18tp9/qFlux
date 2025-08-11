###shared parameters
'''
class to share parameters such as grid, total time, data size, etc.
will use operator to create grid in this class and then pass it to phi, simloop etc.
that way if we change grid size, time - don't have to do it in each file
'''

import os
import numpy as np
from .Dependencies.EquationSolvers.DiffusionSolver_base.Operators import Operator
# from logisticsManager import admin
import datetime

    ####################################
#####
##### Uncomment for Lab Data 
'''
Parameters for the lab data stuff
'''
# Input Info 
labDataName        = os.path.join('..','..','InputData','Erie2.h5')
kappa0             = 1e-4#1e-5
gamma              = 1e-2

timestamp          = str(datetime.datetime.now())[0:-10]
timestamp          = timestamp.replace(" ","_")
timestamp          = timestamp.replace(":","-")


output_path        = os.path.join(os.getcwd(),'Analyses')#,'LakeErieData')

simDataName        = os.path.join(output_path,(timestamp + '.h5'))
outfileName        = simDataName
DataSet            = ''


# inputData = admin.inData(os.path.join('ErieData','Erie_day203.h5'))

# inputData.z -= np.min(inputData.z)

# zMax = np.max(inputData.z)
# zMin = np.min(inputData.z)
# Lz = zMax

# inputData.t_sim -= np.min(inputData.t_sim)
# inputData.t_sim = np.round(inputData.t_sim,2)

# fTime = np.max(inputData.t_sim)
# nt = len(inputData.t_sim)
# dt0 = np.diff(inputData.t_sim)
#params.construct_grid(64+1,(zMax-zMin),gridMin=zMin)

#time parameters

dt0             = 130.5#130.5
nt              = 51#1001#1001
fTime           = 6525#2037


Nz              = 64+1
Lz              = 10.5#5.5
gridMin         = 0
zType           = 'Cosine'
offset          = None#Lz/2
Bc1             = 0
Bc2             = 0

rho0_offset     = 0 
rho0_width      = 0.3


num_iter        = 300

class shared_parameters():

    ####################################
    #####
    ##### Uncomment for Lab Data 
    '''
    Parameters for the lab data stuff
    '''
    # Input Info 

    kappa0             = 7e-2#1e-5

    timestamp          = str(datetime.datetime.now())[0:-10]
    timestamp          = timestamp.replace(" ","_")
    timestamp          = timestamp.replace(":","-")


    output_path        = os.path.join(os.getcwd(),'Analyses')#,'LakeErieData')

    simDataName        = os.path.join(output_path,(timestamp + '.h5'))
    outfileName        = simDataName
    DataSet            = ''

    
    # inputData = admin.inData(os.path.join('ErieData','Erie_day203.h5'))
    
    # inputData.z -= np.min(inputData.z)

    # zMax = np.max(inputData.z)
    # zMin = np.min(inputData.z)
    # Lz = zMax
    
    # inputData.t_sim -= np.min(inputData.t_sim)
    # inputData.t_sim = np.round(inputData.t_sim,2)

    # fTime = np.max(inputData.t_sim)
    # nt = len(inputData.t_sim)
    # dt0 = np.diff(inputData.t_sim)
    #params.construct_grid(64+1,(zMax-zMin),gridMin=zMin)
    
    #time parameters

    dt0             = 28.7#130.5
    nt              = 70#1001#1001
    fTime           = 1980#6525#2037


    Nz              = 64+1
    Lz              = 15.5#5.5
    gridMin         = 0
    zType           = 'Cosine'
    offset          = None#Lz/2
    Bc1             = 0
    Bc2             = 0

    rho0_offset     = 0 
    rho0_width      = 0.3


    num_iter        = 1000


    # ####################################
    # #####
    # ##### Uncomment for Fake Data 
    # '''
    # Parameters for fake data 
    # '''
    # kappa0             = 1e-2
    # rawZ           = np.array([-3.        , -2.22580645, -1.4516129 , -0.67741935,  0.09677419,
    #     0.87096774,  1.64516129,  2.41935484,  3.        ])


    # # output_path        = os.path.join(os.getcwd(),'Analyses')
    # # simDataName        = os.path.join(output_path,'Constant.h5')
    # # DataSet            = 'Constant'
     
    
    # dt0             = 0.01
    # fTime           = 10
    

    # Nz              = 32
    # Lz              = 6
    # zType           = 'Cosine'
    # offset          = 3
    # Bc1             = 0
    # Bc2             = 0

    # rho0_offset     = 0 
    # rho0_width      = 0.3

    # num_iter        = 20
    # #############################################

    def __init__(self):

        #time parameters
        self.datasize           = int(self.fTime/self.dt0) +1

        #grid parameters
        self.construct_grid(Nz = self.Nz, Lz = self.Lz)

        
        
    def construct_grid(self,Nz,Lz,gridMin=0):
        self.Nz     = Nz
        self.Lz     = Lz
        self.grid   = Operator(Lz,Nz,zType=self.zType,
            offset=self.offset,
            BC1=self.Bc1,BC2=self.Bc2,gridMin = gridMin)
        self.grid.z = self.grid.z 
        
    def update_BC2(self,surfaceTemp):

        self.grid.BC2 = surfaceTemp

    def edit_simLength(self,fTime,dt0 = 1e-2):
        self.fTime  = fTime
        self.dt0    = dt0
