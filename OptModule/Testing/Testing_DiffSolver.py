import numpy as np 
import os 
from OptModule.Dependencies.EquationSolvers.DiffusionSolver import Diffusion1D
from OptModule.Dependencies.EquationSolvers.DiffusionSolver_base.Operators import Operator

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data

import matplotlib.pyplot as pl 

if __name__ == "__main__":

    Nz     = 32
    Lz     = 10
    grid   = Operator(Lz,Nz,zType='Cheb',offset=None,BC1=1,BC2=0)



    rho_solver = Diffusion1D(0*grid.z,grid,fTime=1)
    for (t,rho) in rho_solver.loop():
        pl.plot(rho,grid.z)
        pl.pause(0.1)