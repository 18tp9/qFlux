# from ast import operator
from typing import Callable
import numpy as np 
import matplotlib.pyplot as pl
import os 

# from scipy.special import erf

from .DiffusionSolver_base.Operators import Operator
from .DiffusionSolver_base.Base_DIffusionSolver import timestepper    
from ..Tools.Utilities import tools


class Diffusion1D(timestepper):
    


    def __init__(self,rho0,operator,bottomT = None,surfaceT = None,flux = None,hypso=None,kappa0 = 1e-3, dt = 0.01, fTime=10,Nt = 10,explicit_ForceFunc=None) -> None:
        '''
        *** Collect all of the eps data in __init__
        '''
        self.kappa0 = kappa0
        self.fTime  = fTime
        self.dt0    = dt
        self.Nt     = Nt
        self.flux = flux
        self.bottom = bottomT
        self.surface = surfaceT

        ## Grid Parameters
        self.operator   = operator
        
        self.z      = self.operator.z
        self.dz     = self.operator.Lz/(self.operator.Nz-1)
        self.rho0   = rho0
        self.hypso  = hypso


        ##### Initialize
        if explicit_ForceFunc is None: explicit_ForceFunc = lambda t: 0*rho0    ## No forcing by default
        self.initialize_forcing(explicit_ForceFunc = explicit_ForceFunc)
        


    def initialize_forcing(self,explicit_ForceFunc):
        '''
        reset rho data for next loop
        '''
        self.explicit_ForceFunc    = explicit_ForceFunc
        super(Diffusion1D,self).__init__(self.rho0,self.operator,self.kappa0) ## Jason 


    def get_explicit_forcing(self,t):
        #Compute the explicit terms on the RHS of solver
        #Input -- Constant `extra' diffusivity
        #Output -- -F to solve for dT/dt = kappa0*d^2/dz^2 T - F 
        # return -1.0*self.operator.get_grad(F0)
        return self.explicit_ForceFunc(t)

    

    def loop(self):
        '''
        Goes through forwards part of loop
        returns error per time step and time for total error calculation
        '''


        ### loop through rho calculation and error calculation
        for i in range(len(self.dt0)):
        # while self.t<self.fTime: 


            #####################
            # dT_dz = ((-1/self.kappa0)*self.flux[int(self.t/self.dt0)]*(1/4182)*(1/1000))
            # self.operator.BC2 = dT_dz
            # self.operator.BC1 = self.bottom[int(self.t/self.dt0)]
            # self.operator.BC2 = self.surface[int(self.t/self.dt0)]
            
            

            ## define rho
            yield (self.t,self.rho)
                        
            ### get forcing at current time given epsilon - get F'
            self.Fp     = self.get_explicit_forcing(self.t)
                
            ### Time Step size
            # dt   = min(self.dz**2/self.kappa0/2,self.dt0)

            dt = self.dt0[i]

            ###Take a step
            t = self.step(self.Fp,dt0=dt)


        # output the very last rho
        #sort profile
        # self.rho = tools.sortProfile(self.rho,self.z)
        yield(self.t,self.rho)



    # def _store_values(self): 
    #     '''
    #     ** clean up function 
    #     Store the values in arrays 
    #     '''
    #     self.time.append(self.t)
    #     self.rho_arr.append(self.rho)
    #     self.rho_grad.append(self.operator.get_grad(self.rho))


class Diffusion2D(timestepper):
    


    def __init__(self,rho0,operator,hypso=None,kappa0 = 1e-3, dt = 0.01, fTime=10,explicit_ForceFunc=None) -> None:
        '''
        *** Collect all of the eps data in __init__
        '''
        self.kappa0 = kappa0
        self.fTime  = fTime
        self.dt0    = dt

        ## Grid Parameters
        self.operator   = operator
        
        self.z      = self.operator.z
        self.dz     = self.operator.Lz/(self.operator.Nz-1)
        self.rho0   = rho0
        self.hypso  = hypso


        ##### Initialize
        if explicit_ForceFunc is None: explicit_ForceFunc = lambda t: 0*rho0    ## No forcing by default
        self.initialize_forcing(explicit_ForceFunc = explicit_ForceFunc)
        


    def initialize_forcing(self,explicit_ForceFunc):
        '''
        reset rho data for next loop
        '''
        self.explicit_ForceFunc    = explicit_ForceFunc
        super(Diffusion1D,self).__init__(self.rho0,self.operator,self.kappa0) ## Jason 


    def get_explicit_forcing(self,t):
        #Compute the explicit terms on the RHS of solver
        #Input -- Constant `extra' diffusivity
        #Output -- -F to solve for dT/dt = kappa0*d^2/dz^2 T - F 
        # return -1.0*self.operator.get_grad(F0)
        return self.explicit_ForceFunc(t)

    

    def loop(self):
        '''
        Goes through forwards part of loop
        returns error per time step and time for total error calculation
        '''

        ### loop through rho calculation and error calculation
        while self.t<self.fTime:  
            #####################
            ## define rho
            yield (self.t,self.rho)
                        
            ### get forcing at current time given epsilon - get F'
            self.Fp     = self.get_explicit_forcing(self.t)
                
            ### Time Step size
            # dt   = min(self.dz**2/self.kappa0/2,self.dt0)
            dt = self.dt0

            ###Take a step
            t = self.step(self.Fp,dt0=dt)


        # output the very last rho
        yield(self.t,self.rho)