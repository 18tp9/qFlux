### Testing the combination of getEpsilon and Phi
### Completing one full loop
import matplotlib.pyplot as pl
import numpy as np 
import os
import time

from .Dependencies.Minimization.Methods import minimization
from .Dependencies.Tools.Error_Analysis import error_analysis 

from .Dependencies.EquationSolvers.Loop_Forward import Sim_Loop
from .Dependencies.EquationSolvers.Loop_Backwards import phi
from .Dependencies.EquationSolvers.Loop_Data import DataManager
# from OptModule.Parameters import shared_parameters as params
from OptModule import Parameters as params



class Adjoint_Solver():

    def __init__(self,labData,simDataName,grid,Fm1_data = None,fluxArr = None,iteration=-1,
        plot=False,
        min_method="Backtracking",hypso = None):
        
        self.startTime    = time.time()
        # Create basic functionality 
        error_analysis._create_error_files('.')

        self.simDataName = simDataName

        # Grid
        self.grid = grid 
        self.plot = plot

        # Setup the forward Loop 
        self.iteration  = iteration
        self.labData    = labData
        self.simData    = self.load_simData('rho',iteration)
        self.fRestart   = Fm1_data
        self.filterF    = False
        self.smoothing  = True
        self.errorM1    = 100000
        self.errOrder   = 0
        
        ## Minimization parameters
        self.min_method = min_method

        ### Depth AVG parameters
        self.hypso          = hypso
        self.labDataName    = params.labDataName
        self.fluxArr = fluxArr



    def collect_labData(self):
        '''
        load in the LabData 
        '''        
        return DataManager(self.labDataName,format="LabData")


    def load_simData(self,var_name,iteration,filter_data=False):
        '''
        Load the sim data 
        '''
        return DataManager(self.simDataName,var_path=var_name,iteration=iteration,format="SimData",filter_data=filter_data)



    def full_adj_loop(self):
        '''
        Complete one full adjoint loop
        '''
        
        ## If first iteration, save the lab data and set the previous (non-existant) sim data to none
        if self.iteration == -1:
            self.labData.save_data(self.simDataName,t=self.labData.t_sim,z=self.labData.z,T=self.labData.data)
            simData = None
        #Otherwise perform an adjoint loop 
        else:
            simData = self.load_simData('rho',iteration=self.iteration)

            adj_solver = self.adj_iteration(simData,self.labData)

            ## Save the phi_grad for iteration
            adj_data   = adj_solver.get_output_values()
            t_adj_real = self.labData.t_sim #these are the same right now

            z          = self.labData.z #these are the same right now
            phi_grad   = adj_data['phi_grad_flipped']
            simData.set_data(t_adj_real,z,phi_grad) ### Hack to not have to make a new DataManager

        
        # Perform a forward loop 
        if not (self.iteration == -1): self.prev_error = self.fwd_error
        (self.fwd_error,self.fwd_solver) = self.fwd_iteration(self.labData,simData)
        
        ## Update the number of iterations
        self.iteration +=1

        

    def fwd_iteration(self,labData,phiData):
        ### do a forward loop

        if self.iteration == -1:
            ## NO forward data at all
            
            fwd_solver      = Sim_Loop(self.grid,labData,Fm1_data=self.fRestart,fluxArr=self.fluxArr,phiData=None,hypso=self.hypso) 
            self.epsilon    = 0 
            self.epsM1      = 0

            self.fwd_error  = fwd_solver.get_error(0)
            self.errorM1    = self.fwd_error
            
            # delta = np.mean(np.abs(labData.data - np.asarray(fwd_solver.rhoOpt).T))
            # print(delta)
            
            # params.gamma = delta/100
            # raise ValueError
            
            self.errOrder   = minimization.get_exponent(self.fwd_error)
            fwd_solver.storeOpt(fwd_solver.rho_arr,fwd_solver.F_arr)




        elif self.iteration == 0:
           
            ## load F from previous simulation or Fm1 = 0 

            if self.fRestart != None:
                simData_Fm1         = self.load_simData(var_name='F',iteration=self.iteration,filter_data=False)
            else:
                simData_Fm1 = phiData*0

            fwd_solver          = Sim_Loop(self.grid,labData,fluxArr=self.fluxArr,phiData=phiData,Fm1_data=simData_Fm1,hypso=self.hypso) 
            (self.epsilon,self.fwd_error) = self.minimize(fwd_solver.get_error, method=self.min_method,gradPhi=phiData)
            self.optF       = fwd_solver.optF

            

            # decreasing search range scaling if not converging
            # set counter to quit if stuck in while
            counter = 0
            while ((self.epsilon) ==0):
                params.gamma = params.gamma/7
                print('Decreasing step size for linear search...')
                (self.epsilon,self.fwd_error) = self.minimize(fwd_solver.get_error, method=self.min_method,gradPhi=phiData)
            
                counter+=1

                if counter == 100:
                    raise ValueError("NOT CONVERGING.")
                
            # increasing search range scaling if not converging
            # set counter to quit if stuck in while
            counter = 0
            if self.min_method == 'Saddles':
                while (int(self.epsilon) ==1):
                    params.gamma = 5*params.gamma
                    print('Increasing step size for linear search...')
                    (self.epsilon,self.fwd_error) = self.minimize(fwd_solver.get_error, method=self.min_method,gradPhi=phiData)

                    counter+=1

                    if counter == 100:
                        raise ValueError("NOT CONVERGING.")
            
            ## Plot first full iteration J profile
            # alphas,Jvals = minimization.Plotting(fwd_solver.get_error,params.gamma,self.iteration,init=1)
            # fwd_solver.alphas = np.asarray(alphas)
            # fwd_solver.Jvals = np.asarray(Jvals)

            

            if self.fwd_error<self.errorM1:
                fwd_solver.get_error(self.epsilon,printing=False)
                fwd_solver.time = time.time() - self.startTime
                fwd_solver.storeOpt(fwd_solver.rho_arr,fwd_solver.F_arr)
                
            

        else:
            ## full iteration      
            
            
            simData_Fm1         = self.load_simData(var_name='F',iteration=self.iteration,filter_data=self.filterF)

            fwd_solver          = Sim_Loop(self.grid,labData,fluxArr=self.fluxArr,
                                           phiData=phiData,Fm1_data=simData_Fm1,hypso=self.hypso) 
            fwd_solver.alphas = []
            fwd_solver.Jvals = []
            
            (self.epsilon,self.fwd_error) = self.minimize(fwd_solver.get_error, method=self.min_method,gradPhi=phiData)
            
            # # Plotting - save extensive functional view for select iterations
            # if self.iteration <100:
            #     alphas,Jvals = minimization.Plotting(fwd_solver.get_error,params.gamma,self.iteration)
            #     fwd_solver.alphas = np.asarray(alphas)
            #     fwd_solver.Jvals = np.asarray(Jvals)


            # #### HARD CODING METHOD CHANGE ###

            # if self.iteration>5:
            #     self.min_method = 'Backtracking'



            counter = 0
            
            while (int(self.epsilon) ==1) or (self.epsilon == 0):
                
                if int(self.epsilon)==1 and (self.min_method=='Saddles'):
                    params.gamma = 5*params.gamma
                    print('Increasing step size for linear search...')
                    (self.epsilon,self.fwd_error) = self.minimize(fwd_solver.get_error, method=self.min_method,gradPhi=phiData)
                
                    counter+=1
                elif int(self.epsilon)==1:
                    break
                if self.epsilon == 0:
                    if self.min_method == 'Saddles':
                        params.gamma = params.gamma/7
                        print('Decreasing step size for linear search...')
                        (self.epsilon,self.fwd_error) = self.minimize(fwd_solver.get_error, method=self.min_method,gradPhi=phiData)
                    
                    counter+=1
                if counter == 20:
                    raise ValueError("NOT CONVERGING.")

            

            if self.fwd_error<self.errorM1:
                fwd_solver.get_error(self.epsilon,printing=False)
                fwd_solver.time = time.time() - self.startTime
                fwd_solver.storeOpt(fwd_solver.rho_arr,fwd_solver.F_arr)
                fwd_solver.epsOpt = self.epsilon
                


        if self.plot:
            fwd_solver.plot=True
            fwd_solver.get_error(self.epsilon)

        

        self.simData.save_data(self.simDataName,**fwd_solver.get_output_values())
        return (self.fwd_error,fwd_solver)


    def adj_iteration(self,simData_rho,labData):
        ### Solve the adjoint equations
        
        if self.iteration < 0:
            adj_solver = phi(self.grid,simData_rho,labData,self.hypso,flux=self.fluxArr,iteration=-1)

        else:
            adj_solver = phi(self.grid,simData_rho,labData,self.hypso,flux=self.fluxArr)
            
        adj_solver.adj_loop()
        self.simData.save_data(self.simDataName,advance=False,**adj_solver.get_output_values()) 
        return adj_solver


    def minimize(self,func,method="Brent",gradPhi=None):
        '''
        select the minimization method
        '''
        if method=="Brent":
            self.ax = 0
            self.bx = 1
            self.cx =0.5
            return minimization.Brent(self.ax,self.cx,self.bx,func)

        elif method == 'GSS':
            self.ax = 0.0
            self.cx = 1.0
            return minimization.gss(func,self.ax,self.cx)

        elif method=='Backtracking':
            cx  = 1; 
            c1  = 0.6 ; # Wolfe Condition Coefficient
            nu = 0.6 ; # Contraction ratio

            # assert not (gradPhi is None),"phi not provided "
            tau = gradPhi.t_sim[-1] - gradPhi.t_sim[0]
            Lz  = gradPhi.z[-1]    - gradPhi.z[0]
            norm2_gradPhi = c1*np.sum(np.trapz(gradPhi.data**2,gradPhi.z,axis=0)/Lz)/tau
            return minimization.BackTrack(func,cx,nu,c1,norm2_gradPhi)
        
        elif method == 'SimpleStep':
            stepsize = 1.0
            return minimization.SimpleStep(func,stepsize)
        
        elif method == 'BigSteps':
            alpha = 0.8
            c1 = 0.6
            return minimization.BigSteps(func,alpha,c1)
        
        elif method == 'Saddles':
            alpha = 0.7
            c1 = 0.2
            return minimization.SaddleFinder(func,alpha,c1)



    def output_data(self):
        '''
        Get the current forward iteration data 
        -- t,z,rho,F
        '''
        return self.fwd_solver.get_output_values()


if __name__ == "__main__":
    ###loop conditions
    from .Parameters import shared_parameters

    params       = shared_parameters()
    num_iter        = params.num_iter

    Adj_Sol = Adjoint_Solver(
        params.labDataName,params.simDataName,params.grid,
        plot=True,
        )
    for i in range(0,num_iter):
        
        ###run through adjoint loop and save phi values
        print ("Iterations Remaining: " + str(num_iter - i))
        Adj_Sol.full_adj_loop()
        
        err_str = f'{time.time() - Adj_Sol.startTime:0.01f},' + \
                    f'{Adj_Sol.iteration:03i},' + \
                    f'{Adj_Sol.epsilon:0.03e},'+ \
                    f'{Adj_Sol.fwd_error:0.03e},'

        error_analysis.write_error('Error.txt',err_str)

        executionTime = (time.time() - params.startTime)
        print('Execution time in seconds: ' + str(executionTime))



    #tools.plot_comp(epsilon.rho_Solver.time,epsilon.z,labData,num_iter)
    #tools.plot_kappaE(epsilon.rho_Solver.time,epsilon.z,num_iter)


