import numpy as np 
# import h5py as h5 
import matplotlib.pyplot as pl 
import seaborn as sns
import os

import sys
from pathlib import Path

# Add parent directory to sys.path to avoid relative import error
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data
from OptModule.Dependencies.Tools.Error_Analysis import error_analysis
from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from OptModule.Dependencies.Tools.Utilities import tools
from OptModule.Dependencies.Tools.Utilities import shared_parameters
from logisticsManager import admin
# from Figures.LineStyles import RT_lineClrs
sns.set_context('paper', font_scale=1.5)

# fig,axs = pl.subplots(1,2,figsize=(10,6),constrained_layout=True)
# ax0 = axs[0]
# ax1 = axs[1]
fig1,ax = pl.subplots(1,1)
# pl.ion()
fig = pl.figure()
ax0 = pl.gca()
ax2 = ax0.twinx()
# ax.colorbar()
# simDataName = os.path.join('Analyses','IceOn_short_20201211.h5')
# simDataName = os.path.join('Analyses','Cosine_FluxTest2.h5')
# simDataName1 = os.path.join('Analyses','Cheb_k0-1_middleFlux.h5')
# simDataName2 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-3.h5')
# simDataName3 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-4.h5')
# simDataName4 = os.path.join('Analyses','Cheb_k0-1_BoundaryFlux.h5')


# simList = [simDataName1, simDataName2, simDataName3,simDataName4]

# simList = os.listdir(os.path.join('Analyses','time_perFlux_0402'))
# simDataName = os.path.join('Analyses', 'testFluxTest2.h5')
# simDataName2 = os.path.join('Analyses/LakeErieData','2024-02-28_12-00.h5')
# simDataName3 = os.path.join('Analyses/LakeErieData','2024-02-22_14-24.h5')

colourWheel = ['blue','orange','green','red','brown','purple','cyan']

inputDataName = 'FakeData//FakeData_SinForcing.h5'
inputData = admin.inData(inputDataName,var_path='SineForcing/Flux',t_path='SineForcing/t',z_path='SineForcing/z')
inputRho = admin.inData(inputDataName,var_path='SineForcing/rho',t_path='SineForcing/t',z_path='SineForcing/z')
trueF = (inputData.data)# - 1e-2*np.gradient(np.gradient(inputRho.data,axis=0),axis=0)


for i,simDataName in enumerate(['Saddles_kPhi1e-2_long_FakeData_SinForcing.h5']):#,'Analyses//perFlux_time_0423//FABLS_restart_kappaPhi_1e-1.h5']):

    # simDataName = os.path.join('Analyses','time_perFlux_0402',simDataName)

    labData = DataManager(simDataName,
                format="SimData",
                iteration=0,
                t_path= 'LabData/t',
                z_path=  'LabData/z',
                var_path= 'LabData/T',
                relative_path=False,
        )

    rho_lab = labData.data
    params = shared_parameters()
    params.zType = 'Cheb'
    params.construct_grid(len(labData.z),np.max(labData.z)-np.min(labData.z))
    # trueF = 2e-3*params.grid.get_grad(rho_lab)
    # trueF = 0.1*np.flipud(inputData.data)# - 2e-3*params.grid.get_grad(params.grid.get_grad(inputRho.data))
        ### Construct initial variables
    # clr = RT_lineClrs()
    I_list = []
    err_list = []
    err_list2 = []
    err_list3 = []
    smoothList = []
    Ferr_list = []

    Fz_arr = []

    # pl.ion()

    if '' in simDataName:

        for iteration in range(0,500):
                #### Simulated Data 
                try:
                    ErrorData = HDF5data(simDataName,  
                        format="SimData",
                        iteration=iteration  ,
                        var_path='rho'       
                        )

                    (t_sim,z_sim,rho_sim) = ErrorData.loadData()
                
                    Fdata = HDF5data(simDataName,  
                        format="SimData",
                        iteration=iteration  ,
                        var_path='F'       
                        )
                    
                    (t_sim,z_sim,Fsim) = Fdata.loadData()

                    # fig1,ax = pl.subplots(1,1,figsize=([7,6]))

                    # ctf = ax.contourf(t_sim,z_sim,trueF-Fsim,np.linspace(-0.001,0.001,20),extend='both')
                    # ax.set_xlabel(r'$t$ (s)')
                    # ax.set_ylabel(r'$z$ (m)')
                    # ax.set_title('Iteration ' + str(iteration))
                    # cb = pl.colorbar(ctf,ticks=[-0.001,-0.0005,0,0.0005,0.001])
                    # # cbt = cb.get_ticks()
                    # # cb.set_ticks([-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1],[-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1])
                    # cb.set_label(r'$F-F_{sim}$ ($^\circ$C/ms)')
                    # # pl.show()
                    # # fig.canvas.draw()
                    # fig1.savefig(f'FconvPNGs/figure_{iteration:03}.png')
                    # # pl.pause(0.1)
                    # pl.close('all')

                    ## Find error and plot 
                    errZ = np.trapz((rho_sim - rho_lab).T**2,z_sim)#/(np.amax(z_sim) - np.amin(z_sim))
                    err= (np.trapz(np.array(errZ),np.array(t_sim)))#/(np.amax(t_sim) - np.amin(t_sim))
                    Ferr = np.trapz(np.trapz((Fsim[:,:]-trueF[:,:])**2,z_sim,axis=0))#*np.mean(np.diff(t_sim))*np.mean(np.diff(z_sim))#/(np.amax(z_sim) - np.amin(z_sim)),t_sim)/(np.amax(t_sim) - np.amin(t_sim))
                    
                    # Fz_arr.append(np.trapz((Fsim[:,:-1]-trueF[:,:-1])**2,z_sim,axis=0))
                    
                    # Ferr = np.sum(np.abs(Fsim-trueF)
                except:
                    err = np.nan
                    Ferr = np.nan
                    
                if i ==0:
                    I_list.append(iteration); err_list.append(err); Ferr_list.append(Ferr);
                else:
                    I_list.append(iteration+88); err_list.append(err); Ferr_list.append(Ferr);


        I_list = np.asarray(I_list)
        



        err_list = np.asarray(err_list)
        # err_list = err_list#/np.nanmax(err_list)
        Ferr_list = np.asarray(Ferr_list)

        if 'Backtrack' in simDataName:
            ###diffusive
            ls = '-'
        elif ('FABLS' in simDataName) or ('GSS' in simDataName):
            ### artificial flux
            ls = '-'
        else:
            ls = '--'


        if '25' in simDataName:
            colour = 'orange'
        elif '50' in simDataName:
            colour = 'cyan'
        elif 'GSS' in simDataName:
            colour='m'
        else:
            colour = 'blue'

        # paths = simDataName.split('\\')
        # plotLabel = paths[2][:-3]
        # words = plotLabel.split('_')
        # plotLabel = words[:]
        # plotLabel = ' '.join(plotLabel)
        # plotLabel = plotLabel.replace('$-v',r'$\v')
        # plotLabel = plotLabel.replace('on ',r'on_')
        plotLabel = '1'

        if  ''  in simDataName:
            simDataName = simDataName.replace('Artificial','')
            ax0.semilogy(I_list,err_list,linestyle='-',label=plotLabel,color = colour)#//2])
            ax2.semilogy(I_list,Ferr_list,linestyle='--',color = 'r')#//2])
            
            # ax0.axhline((t_sim[1] - t_sim[0])**3,linestyle = ls,color=colour)
        # else:
            # ax1.semilogy(I_list,err_list,linestyle=ls,label=simDataName)
            
# ax.semilogy(I_list[0:1020],err_list2,'-')
# ax.semilogy(I_list[500:661],err_list3,'o-')

    ## Convergence Rate 
# p0 = np.polyfit(I_list[12:],np.log(err_list[12:]),1)

# ax.plot(x_eg,np.exp(p0[0]*x_eg + p0[1]),'k--')

## Labels
ax0.set_xlabel('Iteration number',fontsize=18)
ax0.set_ylabel(r'Error $J = \int\int(T-T_m)^2 dzdt$')


ax0.set_ylabel(r'Error $J = \int\int(T-T_{sim})^2 dzdt$',color='blue',fontsize=18)
ax0.tick_params(axis='y', labelcolor='blue')

ax2.set_ylabel(r'$\int\int(F-F_{sim})^2 dzdt$',color='red',fontsize=18)
ax2.tick_params(axis='y', labelcolor='red')
# ax.set_ylabel('Error Total')
# ax0.set_xlim([-1,130])
# ax.axvline(88,'k')
# ax0.set_ylim([1e-6,1.2])
# ax0.set_title(r'Artificial Data')

# Legend 
handles, labels = ax0.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
pl.show()
# ax0.legend(by_label.values(), by_label.keys(),fontsize=10,loc='upper right', bbox_to_anchor=(1, 1))
# plt.legend(loc='upper left', bbox_to_anchor=(1, 1))


# pl.show()