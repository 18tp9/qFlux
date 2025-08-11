import numpy as np 
# import h5py as h5 
import matplotlib.pyplot as pl 
import seaborn as sns
import os

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data
from OptModule.Dependencies.Tools.Error_Analysis import error_analysis
from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from OptModule.Dependencies.Tools.Utilities import tools


def getMedianStep(simDataName):

    # simDataName = os.path.join('Analyses',simDataName)

    labData = DataManager(simDataName,
                format="SimData",
                iteration=0,
                t_path= 'LabData/t',
                z_path=  'LabData/z',
                var_path= 'LabData/T',
                relative_path=False,
        )

    I_list = []
    eps_list = []
    gammList = []
    print("Getting median step from optimized case...")
    for iteration in range(0,2000):
            #### Simulated Data 
            try:
                ErrorData = HDF5data(simDataName,  
                    format="SimData",
                    iteration=iteration  ,
                    var_path='epsilon'       
                    )

                (t_sim,z_sim,eps) = ErrorData.loadData()


            except:
                eps = np.nan

            try:
                ErrorData = HDF5data(simDataName,  
                    format="SimData",
                    iteration=iteration  ,
                    var_path='gamma'       
                    )

                (t_sim,z_sim,gamma) = ErrorData.loadData()


            except:
                gamma = np.nan
                
            if eps != 0:
                I_list.append(iteration); eps_list.append(eps); gammList.append(gamma) 

    eps_list = np.asarray(eps_list)
    gammList = np.asarray(gammList)

    return np.nanmedian(eps_list*gammList)