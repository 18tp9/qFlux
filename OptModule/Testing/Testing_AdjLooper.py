
import matplotlib.pyplot as pl
import numpy as np 
import os 
import time

from OptModule.Dependencies.Minimization.Methods import minimization
from OptModule.Dependencies.Tools.Error_Analysis import error_analysis 
from OptModule.Dependencies.Tools.Utilities import tools

from OptModule.Dependencies.EquationSolvers.Loop_Forward import Sim_Loop
from OptModule.Dependencies.EquationSolvers.Loop_Backwards import phi
from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager

from OptModule.Parameters import shared_parameters
from OptModule.Adjoint_Looper import Adjoint_Solver




if __name__ == "__main__":
    ###loop conditions
    params       = shared_parameters()
    num_iter        = params.num_iter

    labData = DataManager(params.labDataName,
        format='LabData',
        t_path=(r'/Timeseries/scales/time'),
        z_path= (r'/Timeseries/scales/z'),
        var_path=(r'/Timeseries/Temperature/T_Chain'),
    )

    Adj_Sol = Adjoint_Solver(
        labData,params.simDataName,params.grid,
        plot=False,
        )

    for i in range(0,num_iter):
        
        ###run through adjoint loop and save phi values
        print ("Iterations Remaining: " + str(num_iter - i))
        Adj_Sol.full_adj_loop()
        
        err_str = f'{time.time() - Adj_Sol.startTime}, \
                    {Adj_Sol.iteration:03d},\
                    {Adj_Sol.epsilon:0.03e},\
                    {Adj_Sol.fwd_error:0.03e},\
            '
        error_analysis.write_error('Error.txt',err_str)

        executionTime = (time.time() - Adj_Sol.startTime)
        print('Execution time in seconds: ' + str(executionTime))

