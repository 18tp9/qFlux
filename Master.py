'''vFlux is a 1D scalar flux parameterization tool that 
uses an adjoint method to determine the turbulent heat 
flux present in observed data. see Pendergast et al. 2025'''

import matplotlib.pyplot as pl
import numpy as np 
import os 
import time

from OptModule.Dependencies.Tools.Error_Analysis import error_analysis 
from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from logisticsManager import admin
from OptModule.Parameters import shared_parameters  #class
from OptModule import Parameters as params          #module
from OptModule.Adjoint_Looper import Adjoint_Solver
from FindMedianStep import getMedianStep

import sys

def main():

    
    inputDataName = sys.argv[1]
    minStyle = sys.argv[2]
    
    if minStyle == 'Saddles' or minStyle == 'GSS':
        outPut = minStyle + '_' + inputDataName
        fraction = None
    else:
        fraction = sys.argv[3]
        outPut = str(fraction) + minStyle + '_' + inputDataName


    return inputDataName,minStyle,fraction,outPut

# raise ValueError

if __name__ == "__main__":

    ### key parameters class
    paramClass       = shared_parameters()

    ### if using previous forcing
    restart =False

    ### if variable z array (SCAMP cast etc.) z is interpolated
    varZ= False

    ### load input data
    inputDataFolder = 'FakeData'
    inputDataName = 'FakeData_SinForcing.h5'

    minStyle='Saddles'

    params.kappa0 = 1.0e-2
    params.gamma = 1.0


    ## for batch files
    # inputDataName,minStyle,fraction,outPut = main()
    # if minStyle == 'GSS':
    #     fraction = 1.0

    ###work around _ comment out later
    # fraction = None

    # if fraction != None:
    #     fName = 'Saddles_' + inputDataName
        
    #     medianStep = getMedianStep(fName)


    #     params.gamma = np.float64(fraction)*np.float64(medianStep)

    # ### Hardcode step size
    # params.gamma = 0.001
    


    
    inputData = admin.inData(os.path.join(inputDataFolder,inputDataName),var_path='SineForcing/rho',t_path='SineForcing/t',z_path='SineForcing/z')
    # inputData = admin.inData(os.path.join(inputDataFolder,inputDataName),var_path='Constant/rho',t_path='Constant/t',z_path='Constant/z')
    # inputData = admin.inData(os.path.join(inputDataFolder,inputDataName),var_path='Timeseries/Temperature/T_Chain',t_path='Timeseries/scales/time',z_path='Timeseries/scales/z')
    # inputData = admin.inData(inputDataName,var_path='rho',t_path='t',z_path='z')
    # inputData = admin.inData(inputDataName,var_path='T',t_path='t',z_path='z')




    params.simDataName = minStyle+'_kPhi1e-2_long_'+inputDataName
    
    # # for batch runs
    # if minStyle == 'SimpleStep':
    #     params.simDataName = fraction + '_' + minStyle+'_'+inputDataName




    inputData.z -= np.min(inputData.z)

    zMax = np.max(inputData.z)
    zMin = np.min(inputData.z)
    params.Lz = zMax
    params.Nz = inputData.z.shape[0]

    paramClass.zType = 'Cheb'
    params.zType = 'Cheb'

    paramClass.construct_grid(params.Nz,(zMax-zMin))
    params.grid = paramClass.grid


    if varZ:
        inputData = admin.interpolateZmat(inputData,params.grid.z)

    inputData.data = np.flipud(inputData.data)

    ### match data sizes for t and z
    inputData.t_sim -= np.min(inputData.t_sim)
    inputData.t_sim = np.round(inputData.t_sim,4)
    params.fTime = np.max(inputData.t_sim)
    paramClass.fTime=np.max(inputData.t_sim)
    params.nt = len(inputData.t_sim)
    params.dt0 = np.diff(inputData.t_sim)

    # inputData.z =params.grid.z

    ### set number of iterations
    params.num_iter = 500

    ### empty forcing for first iteration
    fData = DataManager(t_sim=np.arange(params.nt),z=params.grid.z,data = np.zeros([params.Nz,params.nt]),var_name='T')

    if restart:
        ''' Provide forcing from previous test file to restart
        give path and iteration to restartF'''

        fData = admin.restartF('Analyses//perFlux_time_0423//PABLO_tanh_+_Per._Flux.h5',88)
        # fData.data = fData.data + 0.1*np.random.rand(*[params.Nz,params.nt])





    num_iter        = params.num_iter
    
    
        

    fData.t_sim = inputData.t_sim
    fData.z = params.grid.z
    # params.grid.z = inputData.z



    Adj_Sol = Adjoint_Solver(
        inputData,params.simDataName,params.grid,Fm1_data=fData,
        plot=False, min_method=minStyle
        )

    Adj_Sol.smoothing = False


    error_analysis.write_error('Error.txt','RunTime , Iteration , Epsilon , Error' )

    pl.ioff()

    ### Loop
    for i in range(num_iter+1):

        Adj_Sol.plot = False
        
        ###run through adjoint loop and save phi values
        print ("Iterations Remaining: " + str(num_iter - i))
        Adj_Sol.full_adj_loop()

        executionTime = (time.time() - Adj_Sol.startTime)
        
        err_str = f'{time.time() - Adj_Sol.startTime}, \
                    {Adj_Sol.iteration:03d},\
                    {Adj_Sol.epsilon:0.03e},\
                    {Adj_Sol.fwd_error:0.03e},\
                '
        # error_analysis.write_error('Error.txt',err_str)

        
        print('Execution time in seconds: ' + str(executionTime))
        
        if Adj_Sol.fwd_error < 1e-12:
            break

    with open('ExecutionTime.txt','a') as f:
        f.write(params.simDataName + '\n')
        f.write(str(executionTime)+'\n')

        f.close()


