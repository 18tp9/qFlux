### "Front end" manager
### macros for setup - inputs and parameters etc

import numpy as np
import matplotlib.pyplot as pl
import os
from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager

class admin():
    """logistics manager"""
    def __init__(self, arg):
        pass
    
    @classmethod
    def inData(self,filePath,format='LabData',t_path = 't',z_path = 'z',var_path='T'):
            inputData = DataManager(filePath,
            format=format,
            t_path=t_path,
            z_path= z_path,
            var_path=var_path
            )
            return inputData
        
    @classmethod
    def restartF(self,filePath,iteration,format='SimData',t_path = 't',z_path = 'z',var_path='F'):
            fData = DataManager(filePath,
            format=format,
            iteration = iteration,
            t_path=t_path,
            z_path= z_path,
            var_path=var_path,
            relative_path=True,
            filter_data=False
            )
            return fData
        
    @classmethod
    def interpolateZmat(self,inputData,zGrid):
        
        ''' takes data on 2d zgrid and interpolates/normalizes'''
        
        from OptModule.Dependencies.Tools.Utilities import tools
        interpolated = []

        for i in range(len(inputData.t_sim)):
                sampled = tools.subsample(inputData.z[:,i],inputData.data[:,i],zGrid)
                interpolated.append(sampled)
        
        interpolated = np.asarray(interpolated).T
        inputData.data = interpolated
        inputData.z = zGrid

        inputDataNorm = tools.normalizeMax(inputData.data)
        inputData.data = inputDataNorm[:,:-2]
        inputData.t_sim = inputData.t_sim[:-2]
        
        return inputData
    
    @classmethod
    def dirichletBC(self,inputData):
        '''return arrays of top and bottom points
        for Dirichlet BC'''
        
        bottom = inputData.data[-1,:]
        top = inputData.data[0,:]
        
        return bottom,top
        
    
    
    
    
    
    

    