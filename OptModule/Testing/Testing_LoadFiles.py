import numpy as np 
import os 
from OptModule.Dependencies.EquationSolvers.DiffusionSolver import Diffusion1D

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data


from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from OptModule.Dependencies.Tools.Utilities import tools
from OptModule.Parameters import shared_parameters

import matplotlib.pyplot as pl 





if __name__ == "__main__":

    '''
    Basic load and plot 
    '''
    filepath = os.path.join("C:\\Users\\jason\\OneDrive - Queen's University\\Research\\Projects\\Adjoint Loop Rewrite\\InputData")
    filename = os.path.join(filepath,"LabData_20201211.h5")
    h5 = HDF5data(filename, format="LabData")

    (t,z,d) = h5.loadData()

    pl.contourf(t,z,d)

    
    '''
    Filter data and reset 
    '''
    params       = shared_parameters()

    labData = DataManager(params.labDataName,
        format='LabData',
        var_path=(r'/Timeseries/Temperature/T_Chain'),
        t_path=(r'/Timeseries/scales/time'),
        z_path= (r'/Timeseries/scales/z'),
    )
    

    pl.figure()
    pl.contourf(labData.t_sim,labData.z,labData.data)