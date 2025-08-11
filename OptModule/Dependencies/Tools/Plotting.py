### file for indexing/error calcs functions/plotting etc


import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy import array,arange,exp
from matplotlib.widgets import Slider
import matplotlib.pyplot as pl
import h5py
# from ....Graphics.VideoGen import generate_video
# from ....Graphics.VideoGen import make_gif
import os
from ...Parameters import shared_parameters
from .Operators import CHEB

class tools():
    
    @staticmethod 
    def plot_rho(time,z):
        with h5py.File(shared_parameters.simDataName, 'r') as file:
            graph1 = file['iter_19']['rho']

            pl.figure()
            c = pl.contourf((time), z, graph1,np.linspace(0,6,21))
            pl.title('Temperature as a Function of Time and Water Depth - SolverData')
            pl.xlabel('time [s]')
            pl.ylabel('z [m]')
            pl.ylim(top = 0.2)
            pl.autoscale(False)
            pl.colorbar(c)
            pl.show()

    @staticmethod
    def plot_lab():
        with h5py.File('LabData_20201211.h5','r') as file:
            graph   = file['/Timeseries/Temperature/T_Chain']
            time    = file['/Timeseries/scales/time']
            z       = file['/Timeseries/scales/z']

            pl.figure()
            c = pl.contourf((time), z, graph,np.linspace(0,6,21))
            pl.title('Temperature as a Function of Time and Water Depth - LabData')
            pl.xlabel('time [s]')
            pl.ylabel('z [m]')
            pl.ylim(top = 0.2)
            pl.autoscale(False)
            pl.colorbar(c)
            pl.show()

    @staticmethod
    def plot_comp(time,z,labData,num_iter,method):
        '''
        overlayed contours of solver and labdata
        '''
        newlabData  = []
        newTime     = []
        time        = np.asarray(time)
        time        = time[shared_parameters.trunc_lowerB:(shared_parameters.trunc_upperB + shared_parameters.dt0)]

        if shared_parameters.dt0 != 1:
            for i in range(0,(shared_parameters.fTime + shared_parameters.dt0),shared_parameters.dt0):
                newTime.append(time[i])
                newlabData.append(labData[:,i])
            time    = np.asarray(newTime)
            labData = np.asarray(newlabData).T

        groupname = 'iter_' + str(num_iter-1)
        with h5py.File(shared_parameters.simDataName, 'r') as file:
            graph1 = file[groupname]['rho']

            pl.figure()
            c = pl.contourf((time/3600), z, graph1,np.linspace(0,6,25))
            cmap2 = pl.get_cmap('Reds')
            cL = pl.contour(time/3600,z,labData,np.linspace(0,6,25),alpha = 1.0,cmap = cmap2)
            #pl.title('Temperature as a Function of Time and Water Depth')
            pl.xlabel('time [hrs]')
            pl.ylabel('z [m]')
            #pl.legend()
            #pl.ylim(top = 0.2)
            pl.ylim(top = 0.2,bottom = 0.03)
            pl.autoscale(False)
            pl.colorbar(c,label=('$^\circ$ C'))
            pl.show()


    @staticmethod
    # def plot_kappaE(time,z,num_iter,path,show = True):

    #     #shared_parameters = shared_parameters()
    #     group   = ('iter_' + str(num_iter-1))   

    #     with h5py.File(path) as file:
    #         rho     = file[group]['rho'][...]
    #         rho     = np.asarray(rho)
    #         Fm1     = file[group]['F'][...]
    #         Fm1     = np.asarray(Fm1)
        

    #     time    = np.asarray(time)
    #     time    = time[shared_parameters.trunc_lowerB:(shared_parameters.trunc_upperB + shared_parameters.dt0)]

    #     newTime = []
    #     newFm1  = []

    #     if shared_parameters.dt0 != 1:
    #         for i in range(0,int((shared_parameters.fTime + shared_parameters.dt0)),shared_parameters.dt0):
    #             newTime.append(time[i])
    #             #newFm1.append(Fm1[:,i])
    #         time    = np.asarray(newTime)
    #         #Fm1     = np.asarray(newFm1).T

    #     rho_grad = []
    #     for i in range(rho.shape[1]):
    #         rho_grad.append(np.gradient(rho[:,i]))

    #     rho_grad    = np.asarray(rho_grad).T 
    #     zero_x,zero_y   = np.where(rho_grad == 0)
    #     kappas      = Fm1/rho_grad
        
    #     kappas[zero_x,zero_y] = None
    #     kappas  = np.nan_to_num(kappas)
    #     kappas  = kappas/shared_parameters.kappa0

    #     cont = pl.figure(figsize=[8,5])
    #     pl.rc('axes',labelsize = 10)
    #     c = pl.contourf((time/3600), z, kappas,np.linspace(-0.6,0.6,40))
    #     pl.xlabel('time [hrs]')
    #     pl.ylabel('z [m]')
    #     pl.ylim(top = 0.2,bottom = 0.03)
    #     pl.tick_shared_parameters(labelsize=8)
    #     pl.autoscale(False)
    #     cbar = pl.colorbar(c)
    #     cbar.set_label('Eddy Diffusivity')

    #     if show == False:
    #         pl.close()

    #     return kappas

    def plot_kappaE(time,z,rho,F,show = True):

        time    = np.asarray(time)
        #time    = time[shared_parameters.trunc_lowerB:(shared_parameters.trunc_upperB + shared_parameters.dt0)]

        Nz = len(z)
        Lz = z[Nz-1]
        D,z = CHEB.cheb(Nz -1,xmin=0,xmax=Lz)

        # if shared_parameters.dt0 != 1:
        #     for i in range(0,int((shared_parameters.fTime + shared_parameters.dt0)),shared_parameters.dt0):
        #         newTime.append(time[i])
        #         #newFm1.append(Fm1[:,i])
        #     time    = np.asarray(newTime)
        #     #Fm1     = np.asarray(newFm1).T
        

        rho_grad = []
        rho_grad = np.dot(D,rho)

        # for i in range(rho.shape[1]):
        #     rho_grad.append(np.gradient(rho[:,i]))

        # rho_grad    = np.asarray(rho_grad).T 
        zero_x,zero_y   = np.where(rho_grad == 0)
        kappas      = F/rho_grad
        
        kappas[zero_x,zero_y] = None
        kappas  = np.nan_to_num(kappas)
        kappas  = kappas/shared_parameters.kappa0

        cont = pl.figure(figsize=[8,5])
        pl.rc('axes',labelsize = 10)
        c = pl.contourf((time/3600), z, kappas,np.linspace(0,15,40))
        pl.xlabel('time [hrs]')
        pl.ylabel('z [m]')
        pl.ylim(top = 0.19,bottom = 0.12)
        #pl.tick_shared_parameters(labelsize=8)
        pl.autoscale(False)
        cbar = pl.colorbar(c)
        cbar.set_label('Eddy Diffusivity')

        if show == False:
            pl.close()

        return kappas

    @staticmethod
    def plot_diff(time,z,simData,labData):
        simData = np.asarray(simData)
        labData = np.asarray(labData)
        diff = np.abs(simData-labData)
        d = pl.contourf(time,z,diff,np.linspace(0,6,25))
        pl.title('Difference Between Simulation and Lab Data')
        pl.xlabel('time [s]')
        pl.ylabel('z [m]')
        pl.ylim(top = 0.2)
        pl.autoscale(False)
        pl.colorbar(d)

    @staticmethod
    def profile_plots(lab,time,z,num_iter,pngOutPath,vidOutPath,dt = 1,show = True):
        '''
        subplots depth profiles of rho, Fm1, phi showing progression in time
        '''

        group       = ('iter_' + str(num_iter - 1))
            

        with h5py.File(shared_parameters.simDataName) as file:
            rho     = file[group]['rho'][...]
            rho     = np.asarray(rho)
            phi     = file[group]['phi'][...]
            phi     = np.asarray(phi)
            Fm1     = file[group]['Fm1'][...]
            Fm1     = np.asarray(Fm1)

        labTrunc    = []
        ### take every x number of data points in time so we can see progression - but not overstuffed
        for i in range(0,(shared_parameters.trunc_upperB + shared_parameters.dt0)-shared_parameters.trunc_lowerB,shared_parameters.dt0):
            labTrunc.append(lab[:,i])    

        labTrunc    = np.asarray(labTrunc).T

        fig,ax = pl.subplots(1,3,sharey=True,gridspec_kw={'width_ratios':[3,1,1]},**{'figsize':[8,5]})
        pl.subplots_adjust(wspace=0.3)
        
        ind = 0

        for i in range(rho.shape[1]):
            
            currTime = (time[shared_parameters.trunc_lowerB] + (i*shared_parameters.dt0))/3600
            fig.suptitle(f'Time: {currTime:0.3f} hrs')

            if i % dt == 0 :
                ax[0].set_title('\u03C1')
                ax[0].set_xlim(0,5)
                ax[0].set_xlabel('$^\circ$ C')
                ax[0].set_ylim(0.03,0.2)
                ax[1].set_title('Fm1')
                ax[1].set_xlim(-1e-4,1e-4)
                ax[1].set_ylim(0.03,0.2)
                ax[1].tick_shared_parameters(labelsize = 8)
                ax[2].set_title('\u03A6')
                ax[2].set_xlim(-1000,1000)
                ax[2].set_ylim(0.03,0.2)
                ax[2].tick_shared_parameters(labelsize = 8)

                ax[0].plot(rho[:,i],z)
                ax[0].plot(labTrunc[:,i],z)
                ax[1].plot(Fm1[:,i],z)
                ax[2].plot(phi[:,i],z)

                pl.savefig(os.path.join(pngOutPath,(f'{ind:04}.png')))
                ind += 1

            if i % dt == 0:
                pl.pause(0.01)
            
                ax[0].cla()
                ax[1].cla()
                ax[2].cla()
        generate_video(pngOutPath,'ProfileAnim.mp4',vidOutPath,20)
        make_gif(vidOutPath,pngOutPath,'rho_anim.GIF',150)
        pl.close()

        if show == True:
            pl.show()
            pl.ion()  

    @staticmethod
    def error_plot(log = False,show = True):
        path = os.path.join(shared_parameters.error_path, (shared_parameters.outName + '_errDegrees.txt'))
        print (path)
        error = np.loadtxt(path)
        num_iter = len(error)

        iteration = list(range(1,(num_iter+1)))

        #pl.ion()


        fig = pl.figure()
        if log == True:
            pl.semilogy(iteration,error,'s-')
        else:
            pl.plot(iteration,error,'s-')
        pl.xlabel('Iterations')
        pl.ylabel('Error')

        if show == True:
            pl.show()  
        else:
            pl.close()

        return iteration,error  
        
            
    @staticmethod
    def plot_compF(time,z,labData,num_iter):
        '''
        overlayed contours of solver and labdata
        '''
        newlabData  = []
        newTime     = []
        time        = np.asarray(time)
        

        groupname = 'iter_' + str(num_iter-1)
        with h5py.File(shared_parameters.simDataName, 'r') as file:
            graph1 = file[groupname]['rho']
            SimData = graph1[...]
            
            SimData = SimData[:,:1001]
            

            pl.figure()
            c = pl.contourf((time), z, SimData,np.linspace(-1,1,25))
            cmap2 = pl.get_cmap('Reds')
            cL = pl.contour(time,z,labData,np.linspace(-1,1,25),alpha = 1.0,cmap = cmap2)
            #pl.title('Temperature as a Function of Time and Water Depth')
            pl.xlabel('time [hrs]')
            pl.ylabel('z [m]')
            #pl.legend()
            #pl.ylim(top = 0.2)
            pl.ylim(top = 0.2,bottom = 0.03)
            pl.autoscale(False)
            pl.colorbar(c,label=('$^\circ$ C'))
            pl.show()


    @staticmethod
    def plot_kappaE_F(time,z,num_iter,filename,show = True):
        #for fake data sets

        group   = ('iter_' + str(num_iter-1))   

        with h5py.File(filename) as file:
            rho     = file[group]['rho'][...]
            rho     = np.asarray(rho)
            rho     = rho[:,:1001]
            Fm1     = file[group]['F'][...]
            Fm1     = np.asarray(Fm1)
            Fm1     = Fm1[:,:1001]
        

        time    = np.asarray(time)

        rho_grad = []
        for i in range(rho.shape[1]):
            rho_grad.append(np.gradient(rho[:,i]))

        rho_grad    = np.asarray(rho_grad).T 
        zero_x,zero_y   = np.where(rho_grad == 0)
        kappas      = Fm1/rho_grad
        
        kappas[zero_x,zero_y] = None
        kappas  = np.nan_to_num(kappas)
        kappas  = kappas/shared_parameters.kappa0

        print(kappas)

        cont = pl.figure(figsize=[8,5])
        pl.rc('axes',labelsize = 10)
        c = pl.contourf((time), z, kappas,np.linspace(-1,1,40))
        pl.xlabel('time [hrs]')
        pl.ylabel('z [m]')
        pl.ylim(top = 0.2,bottom = 0.03)
        #pl.tick_shared_parameters(labelsize=8)
        pl.autoscale(False)
        cbar = pl.colorbar(c)
        cbar.set_label('Eddy Diffusivity')
        pl.show()

        if show == False:
            pl.close()

        return kappas
    
    @staticmethod
    def Kprofiles_timeAvg(z,time,rho,F,t1,t2,salinity):
        '''
        time avg profile plots for kappa
        *provide time bounds in hours*
        '''
        pl.ion()
        t1 = int((t1*60*60)/shared_parameters.dt0)
        t2 = int((t2*60*60)/shared_parameters.dt0)
    
        time    = np.asarray(time)

        Nz = len(z)
        Lz = z[Nz-1]

        D,z = CHEB.cheb(Nz -1,xmin=0,xmax=Lz)

        rho_grad = []
        rho_grad = np.dot(D,rho)

        zero_x,zero_y   = np.where(abs(rho_grad) < 1e-1)
        kappas      = F/rho_grad
        
        kappas[zero_x,zero_y] = 0
        kappas  = np.nan_to_num(kappas)
        kappas  = kappas/shared_parameters.kappa0
        kappas = np.asmatrix(kappas)

        kappas = kappas[:,t1:t2]
        # time = time[t1:t2]

        newKappas = []
        for row in range(kappas.shape[0]):
            kappas[row,:] = np.mean(kappas[row,:])
            newKappas.append(np.mean(kappas[row]))

        newKappas = np.asarray(newKappas)
        # newKappas[newKappas<0] = None

        # z[newKappas == None] = None
        normZ = (z - min(z))/(max(z) - min(z))
        print(normZ)

        pl.plot(newKappas,normZ)
        pl.xlim(0,max(newKappas))
        pl.ylim(0,1)
        pl.xlabel('\u03BA\u2091 / \u03BA\u2080')
        pl.ylabel('z (m)')
        

        return newKappas

    