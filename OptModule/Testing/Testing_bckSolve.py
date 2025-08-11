import numpy as np 
import os 

from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from OptModule.Dependencies.EquationSolvers.Loop_Backwards import phi
from OptModule.Parameters import shared_parameters


import matplotlib.pyplot as pl 

if __name__ == "__main__":

    #Lab Data
    filepath = os.path.join("C:\\Users\\jason\\OneDrive - Queen's University\\Research\\Projects\\Adjoint Loop Rewrite\\InputData")
    filename = os.path.join(filepath,"LabData_20201211.h5")
    LabData = DataManager(filename, format="LabData")


    #Parameters
    params       = shared_parameters()
    SimData = DataManager(params.simDataName,var_name='rho',format="SimData",iteration=0)


    #Create a forward solver 
    adj_solver = phi(params.grid,SimData,LabData)
    adj_solver.adj_loop()
    SimData.save_data(params.simDataName,advance=False,**adj_solver.get_output_values()) 
    