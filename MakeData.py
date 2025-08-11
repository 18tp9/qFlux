'''
Script to make Fake Data for validation
'''

import numpy as np 
import matplotlib.pyplot as pl
import h5py as hdf5
import random
import os

from scipy.special import erf

from OptModule.Dependencies.EquationSolvers.DiffusionSolver_base.Base_DIffusionSolver import timestepper
from OptModule.Dependencies.EquationSolvers.DiffusionSolver_base.Operators import Operator
from scipy.interpolate import interp1d
from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data







class Diffusion1D(timestepper):

    def get_explicit_forcing(self,kz):
        '''
        Compute the explicit terms on the RHS of solver
        Input -- Constant `extra' diffusivity
        Output -- -F to solve for dT/dt = kappa0*d^2/dz^2 T - F 
        '''
        # return -1.0*F
        rho_z = self.op.get_grad(self.rho)
        F     = self.op.get_grad((kz)*rho_z)
        # F     = 0
        return -1.0*F
    
    def setBC(self):
        self.op.BC1 = 1
        self.op.BC2 = -1


# simDataName = os.path.join('Analyses','Erie2Best.h5')


# fluxData = HDF5data(simDataName,  
#         format="SimData",
#         iteration=305 ,
#         var_path='F'       
#         )
# (_,_,F) = fluxData.loadData()

# CompData = HDF5data(simDataName,            
#         t_path=('LabData' + '/t'),
#         z_path= ('LabData' + '/z'),
#         var_path=('LabData' + '/T'),
#         # format='LabData',
#         )

# (t_lab,z_lab,rho_lab) = CompData.loadData()

if __name__=="__main__":


    pl.ion(); pl.show()
    outputFile = 'FakeData/extraDiff2_ke-1.h5'
    
    #   Grid Data
    Nz,Lz   = 64+1,1
    kappa0 = 1e-2
        
    grid = Operator(Lz,Nz,zType='Cheb',offset=Lz/2)
    z    = grid.z
    dz   = Lz/Nz
    # dz = np.mean(np.diff(z))

    ## Initial Condition
    rho0 = -(erf((z-0)/0.1))
    # rho0 = np.linspace(1,-1,65)
    # rho0 = z*0
    # rho0 = (-z**2 + 0.25)*10
    

    # Time Stepping variables
    dt0 = 5e-4
    fTime = 10.0
    
    # tGrid = Operator(fTime,int(fTime//dt0),zType='Cheb')
    # tCheb = tGrid.z
    # dt = np.diff(tCheb)

    
    for testCase in [
        'Constant']:
    
        ## Solver Data
        # Kz_Test1  = 0.5*(1+np.tanh(z/0.1))*mystery_kappa
        # Kz_Test2  = 0.5*(1+np.tanh(z/0.1))*mystery_kappa * (t)
        if testCase == 'Constant':
            mystery_kappa = 2e-3
            forcingFunc = lambda t: 0*t +  0*z + mystery_kappa
        elif testCase == 'Var_z':
            mystery_kappa = 1e-5
            forcingFunc = lambda t: 0*t +  (1+np.tanh(z/0.1))*mystery_kappa
        elif testCase == 'Var_zt':
            mystery_kappa = 1e-5
            forcingFunc = lambda t: 0.5*t*(1+np.tanh(z/0.1))*mystery_kappa


        ### Construct the solevr 
        rho_Solver  = Diffusion1D(rho0,grid,kappa0)

        rho_Solver.setBC()



        ## Initial output
        print('Starting Forward Solve ')
        print(' Starting to time step ...\n' + 
            f'Grid Size Nz={Nz} \n'  \
            + '-------')


        ## Output Data 
        t_Array = []
        rho_Array = []
        fluxArr = []

        i=0       

        while rho_Solver.t<=fTime:
            i+=1
            ### get forcing at current time given epsilon - get F'
            Kz     = forcingFunc(rho_Solver.t)
            #Fp     = rho_Solver.get_explicit_forcing(Kz) + 50*np.diff(z,append=0.5)*np.sin(2*t)*(-z**2 + 0.25)
            Fp     = rho_Solver.get_explicit_forcing(Kz) #+ rho_Solver.op.get_grad(0.5*np.sqrt(np.pi)*erf((z)/0.25))*np.sin(2*t)
            Fp1    = rho_Solver.get_explicit_forcing(Kz) #+ (0.5*np.sqrt(np.pi)*erf((z)/0.25))*np.sin(2*t)
            fluxArr.append(Fp1)

            ### Time Step size
            dt   = min(dz**2/(kappa0+2e-3)/2,dt0)

            ###Take a step
            t = rho_Solver.step(Fp,dt0=dt)
            
            # Output Details	
            print(f'time = {t:0.4f}, min rho = {np.amin(rho_Solver.rho):0.4f}, ' )


            if False: # rho_Solver.jj %100 == 0:
                pl.figure(1); pl.clf()
                pl.plot(rho0,z,'k')
                pl.plot(rho_Solver.rho,z,label='T')
                pl.plot(-erf(z/np.sqrt(4*(rho_Solver.kappa0+mystery_kappa)*(t) + 0.3**2)),z,'k--')
                
                pl.legend()
                pl.xlabel(r'$\rho$')
                pl.ylabel('z')
                pl.title(f'Step = {rho_Solver.jj},Time = {t:0.4f}')
                pl.xlim([-1.1,1.1])
                pl.pause(0.1)

                # pl.waitforbuttonpress()
            if i%10 ==0:
                coeff1, coeff2 = 10*np.abs(np.random.randn()),np.random.randn()
            flux = 10*np.diff(z,append=0.5)*np.sin(2*t)*(-z**2 + 0.25)
            # flux = 10*np.diff(z,append=0.5)*np.sin(2*t)*(z**2)*10
            # flux = (z*0)

            # noise = np.random.normal(0.25, 0.05, rho_Solver.rho.shape)

            # flux[0] = 1*np.abs(np.sin(2*t))

            
            
            # fluxArr.append(flux)
            t_Array.append(rho_Solver.t)
            rho_Array.append(rho_Solver.rho)# + flux)

        

        

        with hdf5.File(outputFile,'a') as h5:
            
            ## Base Case 
            # save_data = {'t': np.array(t_Array)[::(len(t_Array)//512)],
            #         'z':z,
            #         'rho':np.array(rho_Array).T[:,::(len(t_Array)//512)],}
            save_data = {'t': np.array(t_Array),
                    'z':z,
                    'rho':np.array(rho_Array).T,
                    'trueF':np.array(fluxArr).T}

            h5.require_group(testCase)
            for key in save_data.keys():
                h5[testCase][key] = save_data[key]


            # ## Base Case + Noise 
            # save_data = {'t': np.array(t_Array)[::1024],
            #         'z':z,
            #         'rho':np.array(rho_Array).T[:,::1024] + np.random.normal(0,0.1,save_data['rho'].shape),
            #         }
            # h5.require_group(testCase+'_Noisy')
            # for key in save_data.keys():
            #     h5[testCase+'_Noisy'][key] = save_data[key]

            # ## Base Case + Periodic 
            # save_data = {'t': np.array(t_Array)[:512:],
            #         'z':z,
            #         'rho':np.array(rho_Array).T[:,:512:] + 0.5*np.sin(20*np.asarray(t_Array[:512:])),
            #         }
            # h5.require_group(testCase+'_Sin')
            # for key in save_data.keys():
            #     h5[testCase+'_Sin'][key] = save_data[key]

            # ## Base Case + Subsampled
            # save_data = {'t': np.array(t_Array),
            #         'z':np.append(z[::4],3),
            #         'rho':np.append(np.array(rho_Array).T[::4,:],np.array(rho_Array).T[-1:,:],axis=0) ,
            #         }
            # h5.require_group(testCase+'_Subsample')
            # for key in save_data.keys():
            #     h5[testCase+'_Subsample'][key] = save_data[key]


            # ## Base Case + Interpolated
            # interpolate = interp1d(save_data['z'],save_data['rho'],axis=0,kind='linear')
            # newRho = interpolate(z)
            # save_data = {'t': np.array(t_Array),
            #         'z':z,
            #         'rho':newRho ,
            #         'z_raw':save_data['z'],
            #         }
            # h5.require_group(testCase+'_Interpolated')
            # for key in save_data.keys():
            #     h5[testCase+'_Interpolated'][key] = save_data[key]



# np.savez('Fake_Data_Test2.npz',
#         t_Array=np.array(t_Array),
#         z = z,
#         rho_Array=np.array(rho_Array).T,
#         )