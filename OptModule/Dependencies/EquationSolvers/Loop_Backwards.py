#class for phi to step backwards in time

from ctypes import _SimpleCData
import os
import numpy as np
import matplotlib.pyplot as pl


from .DiffusionSolver import Diffusion1D
from .DiffusionSolver_base.Operators import Operator
from ..Tools.Utilities import tools

# from ...Parameters import shared_parameters as params
from ... import Parameters as params


class phi():


    def __init__(self,grid,SimData,labData,hypso,flux=None,iteration = 1,plot=False):

        # Initialize 
        self.plot       = plot
        self.operator   = grid
        self.hypso      = hypso
        self.iteration  = iteration
        self.flux       = flux
        self.iteration  = iteration

        self.init_data(labData,SimData)      
        self.fTime = np.max(labData.t_sim)
        self.dt0 = np.diff(labData.t_sim)[::-1]

        self.init_adjSolver(1.0,self.dt0,self.fTime); 
        
        

    def init_data(self,labData,SimData):
        '''
        reset rho data for next loop
        '''
        self.SimData = SimData

        self.labData = labData


    ################################
    '''
    Initialize Diffusion Solver 
    '''
    ################################
    def init_adjSolver(self, kappa0, dt, fTime):
        '''
        reset for next loop
        --- kappa0,dt,dTime are scalars
        --- Note:: Diffusion 1D needs a function (forcingFunc) defined below 
        '''
        
        self.adj_time = []
        self.phi_arr  = []


        # rho     = self.SimData.get_Data(-1,reverse=True)
        # rho_m   = self.labData.get_Data(-1,reverse=True)
        # rho_L   = tools.grid_change(self.SimData.z,rho,self.labData.z)

        # phi0 = 2.0*(rho_L-rho_m)
        # phi0 = tools.grid_change(self.labData.z,phi0,self.SimData.z)

        # self.phi0 = phi0
        self.operator.BC1 = 0
        self.operator.BC2 = 0
        
        self.z              = self.operator.z
        self.phi0           = 0*self.z
        # self.phi0 = self.forcingFunc(0)
        self.phi_Solver     = Diffusion1D(self.phi0,self.operator,
                                kappa0=kappa0, dt = dt, fTime=fTime,flux=self.flux,
                                explicit_ForceFunc=self.forcingFunc,Nt=params.nt
                                ) ## Jason 

        

    def get_realtime(self,t_adj):
        '''
        converting t'(adjoint) to t(forward) for pulling rho values
        '''
        realtime = params.fTime - t_adj

        return realtime

        
    def forcingFunc(self,t_adj):
        '''
        Forcing function for the diffusion equation 
         --- needed : (rho - rho_m) to solve for phi
        ## Update to use regularized Gaussian at data points
        '''
        

        ###convert back to real(forwards) time
        
        t = self.get_realtime(t_adj)

        zL = self.SimData.z

        rho     = self.SimData.get_Data(t,reverse=True)
        rho_m   = self.labData.get_Data(t,reverse=True)
        
        if np.all(self.labData.z == self.z):
            ## If the lab data and the simulated data are on the same grid, just take the difference pointwise
            Forcing_Sum = 2.0*(rho-rho_m)
        else:

            Forcing_Sum = 0*zL
            dz_lab_vec = np.append(np.diff(self.labData.z),self.labData.z[-1] - self.labData.z[-2])
            
            for ii,z0 in enumerate(self.labData.z):#[:,t]):
                # regularizing coefficient -- smoothed over 1-ish grid cells
                delta = np.min([dz_lab_vec[ii],dz_lab_vec[ii-1]])#/5 if ii>0 else dz_lab_vec[ii]/5 ## 3 standard deviations is 99% reduction to the next grid point

                ### This is a Gaussian Function with a peak equal to the difference and whose magnitude decreases to 0.01 by the next grid cell
                forcingZ   = 2.0*(rho-rho_m)[ii]*np.exp(-(zL-z0)**2/2/delta**2)#/np.sqrt(2*np.pi*delta**2)

                Forcing_Sum = Forcing_Sum +  forcingZ


        self.iteration +=1

        return Forcing_Sum

    ################################
    ################################



    ################################
    '''
    Adjoint Loop 
    '''
    ################################
    def adj_loop(self):
        '''
        Goes through adjoint loop and updates arrays
        '''

        # # #set phi0 to error at last point in time
        # rho     = self.SimData.get_Data(-1,reverse=True)
        # rho_m   = self.labData.get_Data(-1,reverse=True)
        # # rho_L   = tools.grid_change(self.SimData.z,rho,self.labData.z)

        # phi0 = 2.0*(rho-rho_m)
        

        # self.phi0 = phi0


        for (t_adj,phi) in self.phi_Solver.loop():  
    


            ### store values for later saving
            
            self._store_values(t_adj,phi)


            ### Plot out the data for testing purposes.
            if (self.plot) and (self.phi_Solver.jj % 10 == 0):
                self._plotting_function()
            



    def _store_values(self,t_adj,phi): 
        '''
        ** clean up function 
        Store the values in arrays 
        '''
        self.adj_time.append(t_adj)
        self.phi_arr.append(phi)

    def get_output_values(self):
        '''
        Return the output data as a dictionary
        '''
        return {
            'adj_time':np.array(self.adj_time),
            # 'adj_z':self.z,
            'phi':(np.array(self.phi_arr).T),
            'real_time':self.get_realtime(np.array(self.adj_time)),
            'phi_grad_flipped':self.get_FunctionalGradient(),
        }


    def get_FunctionalGradient(self):
        '''
        ????
        Function to output the gradient of phi
        '''
        phi_arr = np.asarray(self.phi_arr).T
        phi_grad = self.operator.get_grad_array(phi_arr)
        return np.fliplr(phi_grad)


    def _plotting_function(self):
        '''
        ** Clean up function
        Plot the value of phi
        '''
        pl.figure(1); pl.clf()
        pl.plot(self.phi0,self.z,'k')
        pl.plot(self.rho/np.max(np.abs(self.rho)),self.z,label='T')
        
        pl.legend()
        pl.xlabel(r'$\phi$')
        pl.ylabel('z')
        pl.title(f'Step = {self.jj},Time = {self.t:0.4f}')
        pl.xlim([-1.1,1.1])
        pl.pause(0.1)


            

        


        


