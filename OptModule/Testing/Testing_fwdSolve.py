import numpy as np 
import os 

from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from OptModule.Dependencies.EquationSolvers.Loop_Forward import Sim_Loop
from OptModule.Parameters import shared_parameters


import matplotlib.pyplot as pl 

if __name__ == "__main__":

    #Lab Data
    filepath = os.path.join("C:\\Users\\jason\\OneDrive - Queen's University\\Research\\Projects\\Adjoint Loop Rewrite\\InputData")
    filename = os.path.join(filepath,"LabData_20201211.h5")
    LabData = DataManager(filename, format="LabData")
    SimData = DataManager(filename,format="SimData",NO_LOAD=True)

    #Parameters
    params       = shared_parameters()

    #Create a forward solver 
    fwd_solver      = Sim_Loop(params.grid,LabData,phiData=None,plot=False) 
    error           = fwd_solver.get_error(0)

    LabData.save_data(params.simDataName,t=LabData.t_sim,z=LabData.z,T=LabData.data)
    SimData.save_data(params.simDataName,**fwd_solver.get_output_values())
    # SimData.save_data(params.simDataName,**fwd_solver.get_output_values())
    