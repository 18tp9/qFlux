
#simulation loop - runs through forward loop
import os
import numpy as np
import matplotlib.pyplot as pl
from datetime import datetime, timedelta


from .DiffusionSolver import Diffusion1D

from ..Tools.Error_Analysis import error_analysis as err_calc
# from ...Parameters import shared_parameters as params
from ... import Parameters as params
# from ..Tools.Dusolt_data import Dusolt
from ..Tools.Utilities import tools



class Sim_Loop():

    def __init__(self,grid,labData,fluxArr=None,phiData=None,hypso=None,Fm1_data=None,plot=False):
        '''
        *** Collect all of the eps data in __init__
        '''
        #file names as attributes makes interpolation alot easier
        self.errD_path      = 'fwdIteration_error.txt'
        self.errDeg         = 1

        self.plot = plot
        self.img  = False
        ## Initialize the data
        self.phiData = phiData if phiData is not None else 0*labData
        self.labData   = labData
        self.Fm1_data = Fm1_data if Fm1_data is not None else 0*labData
        self.hypso     = hypso
        self.fluxArr   = fluxArr
        self.epsOpt    = 0 
        self.time      = 0
        

        self.alphas = []
        self.Jvals = []

        self.evaluations = 0


        params.dt0 = np.diff(self.labData.t_sim)
        params.fTime = np.max(self.labData.t_sim)

        #### Diffusion solver. 
        self.operator       = grid
        rho0                = self.labData.data[:,0]
        self.init_fwdSolver(rho0,params.kappa0,params.dt0,params.fTime)       ### Define rho_Solver


        ## Iteration Parameters
        self.Fm1    = 0*self.z
        self.Fp     = 0*self.z
        self.rho_lab = []
        self.rhoOpt         = np.zeros(labData.data.shape).T
        self.optF           = np.zeros(labData.data.shape).T

        self.forward_iteration_number = -1
        self.errorMin                 = 1e15

        #for plotting GLM
        self.day1 = datetime(1983,1,2)
    

    def initialize(self):
        '''
        Reset the different data Readers
        '''
        # self.labData.reset()
        if self.phiData is not None: 
            self.phiData.reset()

        if self.Fm1_data is not None: 
            self.Fm1_data.reset()

        self.init_fwdSolver(self.rho0,params.kappa0,params.dt0,params.fTime)

    ################################
    '''
    Initialize Diffusion Solver 
    '''
    ################################
    def init_fwdSolver(self,rho0,kappa0,dt,fTime,explicit_ForceFunc=None):
        '''
        taking labData initial point for truncated series
        '''
        self.sim_time       = []
        self.rho_arr        = []
        
        self.F_arr          = [] 
        self.fluxStore      = []
        self.fSmoothed      = False

        self.rho0           = rho0
        self.z              = self.operator.z
        self.rho_Solver     = Diffusion1D(self.rho0,self.operator,
                                kappa0=kappa0,flux=self.fluxArr,
                                dt = dt, fTime=fTime,explicit_ForceFunc=explicit_ForceFunc,Nt = params.nt,
                                hypso = self.hypso
                                ) ## Jason 
        


    def get_newF(self,epsilon):
        '''
            Gets z X t matrix F value by pulling previous Fm1 and adding given epsilon X  z by t phi_grad
            '''

        if self.Fm1_data is not None:
            return  (self.Fm1_data + ((-1.0*epsilon*params.gamma)*self.phiData))
        
        return None if self.phiData is None else ((-1.0*epsilon*params.gamma)*self.phiData)


    def __explicit_ForceFunc(self,t):
        '''
        Forcing function for the diffusion equation 
        '''
        F0     = self.F.get_Data(t)     ### Get F at specific time
        
        ### depth AVG
        # F_A    = self.__update_F_A()
        
        return -1.0*self.operator.get_grad(F0)
        



    def init_forcing(self,epsilon):
        '''
        Initialize the diffusion solver to the appropriate forcing 
        '''
        self.F = self.get_newF(epsilon)
        if self.F is None:
            self.Forcing_function = lambda t : self.__update_F_1()
        else:
            self.Forcing_function = lambda t : self.__explicit_ForceFunc(t)

        self.rho_Solver.initialize_forcing(self.Forcing_function)

    def __update_F_1(self):
        '''
        returns area forcing
        '''
        return params.kappa0*(self.operator.get_grad(self.rho_Solver.rho))#*(self.operator.get_grad(self.hypso))

    ################################
    ################################


    ################################
    '''
    Error Calculations
    '''
    ################################
    def get_error(self,epsilon,printing=True):
        '''
        Total error calculation for a given kappa value
        -runs loop and saves data
        '''
        self.initialize()
        self.init_forcing(epsilon)

        ### Loop over the error 
        err_pointwise = []
        meanErr_pointwise = []
        err_int = []

        self.evaluations+=1

        for (t,rho) in self.rho_Solver.loop():

            
            # get surface temperature at time t and set BC2
            rho_lab = self.labData.get_Data(t,reverse=False)

            surfTemp = rho_lab[0]

            bottTemp = rho_lab[-1]
            
            # gradRho = self.operator.get_grad_array(rho_lab)
            # print(gradRho.shape)

            self.operator.BC2 = bottTemp 
            self.operator.BC1 = surfTemp
            # self.operator.BC2 = gradRho[-1]
            # print(self.operator.BC2)
            # self.operator.BC1 = gradRho[0]
            # print(self.operator.BC1)
            
            # rho = tools.sortProfile(rho)
            
            
            # if t<params.fTime:  #time as index - must be one less size
                
            self._store_values(t,rho)


            # Compute the error if LabData is at fixed points (not interpolated)

            err_pointwise.append(
                err_calc.error_at_fixed_points(rho,self.z,rho_lab,self.labData.z))
                # err_calc.error_at_fixed_points(rho,self.z,rho_lab,self.labData.z[:,t]))
            
            # meanErr_pointwise.append(np.mean(err_calc.error_at_fixed_points(rho,self.z,rho_lab,self.labData.z[:,t])))
            # Compute the integrate error 
            err_int.append(
                err_calc.error_integrated(rho,self.z,rho_lab,self.labData.z))
                # err_calc.error_integrated(rho,self.z,rho_lab,self.labData.z[:,t]))



            # ### Plot out the data for testing purposes.
            if (self.plot):# and (self.rho_Solver.jj %1 == 0):
                self._plotting_function(t,rho,rho_lab)

            
        ####### 
        ## Return the error 

        ### Saving Data to Arrays
        err_J           = (np.trapz(np.array(err_int),np.array(self.sim_time)))
        # err_RMSE        = np.sqrt(np.trapz(np.array(err_pointwise),np.array(self.sim_time))/(self.rho_Solver.fTime*params.Lz))
        err_RMSE        = np.sqrt(np.mean(np.asarray(err_pointwise)**2))
        # err_int         = (np.trapz(np.array(err_int),      np.array(self.sim_time))/self.rho_Solver.fTime)
        # meanErr         = sum(meanErr_pointwise)/self.rho_Solver.fTime

        

        # error_str = f'Epsilon = {epsilon:0.06e} ==> Error = {err_int:0.03e} ==> Degree Error = {np.sqrt(err_int):0.03e} ==> Pointwise error = {np.sqrt(err_pointwise):0.03e}'
        error_str = f'Epsilon = {epsilon:0.06e} ==> Functional J Error = {err_J:0.03e} ==> RMSE error = {err_RMSE:0.03e}'
        if printing:
            print(error_str)
            err_calc.write_error(self.errD_path,error_str)

        ### Store best rho, F within sub-iterations

        # if err_J<self.errorMin:
        #     self.errorMin = err_J
        #     self.errMinRMSE = err_RMSE
        #     self.epsilon  = epsilon
        #     # print(self.epsilon)
        #     self.__storeOpt()




        # return np.sqrt(err_pointwise)
        return err_J
    
    def storeOpt(self,rho,F):
        self.rhoOpt = rho
        self.optF   = F
    
    def __labGridInterp(self):

        z_sim = self.rho_Solver.z
        rho_sim_lab = np.zeros(self.labData.data.shape)

        for j in range(max(self.sim_time)):

            rho_sim = np.asarray(self.rho_arr).T[:,j]
            # rho_lab = self.labData.data[:,j]
            rawZ = self.labData.z
            # rawZ = self.labData.z[:,j]

            ### Weighted interpolation
            for i,zL0 in enumerate(rawZ):

                arg_upper = np.argmax(z_sim>zL0) if np.argmax(z_sim>zL0)>0 else -1
                z1        = z_sim[arg_upper] 
                z0        = z_sim[arg_upper-1] if arg_upper >0 else z_sim[0]
                
                w1 = (z1 - zL0)/(z1-z0) 
                w0 = (zL0 - z0)/(z1-z0) 
                rho_sim_lab[i,j] = rho_sim[arg_upper-1]*w1 + rho_sim[arg_upper]*w0

        return rho_sim_lab

    def _store_values(self,t,rho): 
        '''
        ** clean up function 
        Store the values in arrays 
        '''
        
        self.sim_time.append(t)
        self.rho_arr.append(rho)
        if self.F is None:
            self.F_arr.append(0)    
        else:
            self.F_arr.append(self.F.get_Data(t))
  

    def get_output_values(self):
        return {
            'rho':(np.array(self.rhoOpt).T),
            # 'rho':(np.array(self.rho_arr).T),

            # 'rho_labGrid':(self.__labGridInterp()),
            'z':self.z,
            'sim_time':np.array(self.sim_time),
            'F':(np.array(self.optF).T),
            # 'F':(np.array(self.F_arr).T),
            'kappa0': params.kappa0,
            'gamma': params.gamma,
            'epsilon': self.epsOpt,
            'tComp': self.time,
            'evals': self.evaluations,
            'alphas':self.alphas,
            'Jvals':self.Jvals
            # 'RMSE': self.errMinRMSE
            # 'H_d':np.asarray(self.fluxStore),
            # 'H_s':np.asarray(self.fluxStore)*(4182*1000*(1/(3600*24))),
            # 'rho_grad':np.asarray(self.operator.get_grad(np.asarray(self.rhoOpt).T))[0,:]

        }




    ################################
    ################################




    def _plotting_function(self,t,rho_sim,rho_lab):
        '''
        ** Cleanup function 
        Plotting function 
        '''
        if t%10==0:
            pl.ion()
            pl.figure(1); pl.clf()
            rho_sim = np.flipud(rho_sim)
            # pl.plot(self.rho0,self.z,'k')

            plotting_z = -1*(self.z - max(self.z))  #so 0 is at surface

            pl.plot(rho_sim,plotting_z,label='T')
            pl.plot(rho_lab,self.labData.z,label='Ref. T')
            # pl.plot(range(0,26),np.zeros(25),color='k')
            pl.ylim(-0.1,10.1)
            pl.gca().invert_yaxis()
            # pl.fill_betweenx(plotting_z,rho_sim,rho_lab,color='grey',alpha=0.2)
            pl.legend(loc = 'lower right')

            pl.xlabel(r'$^\circ$C')
            pl.ylabel('z')

            date = self.day1 + timedelta(days=t)
            
            # pl.title(date.strftime('%B %d'))
            pl.xlim([10,26])
            pl.show()
            pl.pause(0.01)

            if self.img == True:     #image saving
                #img_ind+=1
                path_name = (f'C:/Users/pende/Documents/Summer22/Introduction/rhoAnimation/{self.rho_Solver.jj:04}.png')
                pl.savefig(path_name)
        # pl.pause(0.001)
