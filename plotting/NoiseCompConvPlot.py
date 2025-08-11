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

# from Figures.LineStyles import RT_lineClrs
sns.set_context('paper', font_scale=1.5)

fig,axs = pl.subplots(1,2,figsize=(14,8),constrained_layout=True)
ax0 = axs[0]
ax1 = axs[1]

fig = pl.figure(figsize=(12,8))
ax0 = pl.gca()
# ax.colorbar()
# simDataName = os.path.join('Analyses','IceOn_short_20201211.h5')
# simDataName = os.path.join('Analyses','Cosine_FluxTest2.h5')
# simDataName1 = os.path.join('Analyses','Cheb_k0-1_middleFlux.h5')
# simDataName2 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-3.h5')
# simDataName3 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-4.h5')
# simDataName4 = os.path.join('Analyses','Cheb_k0-1_BoundaryFlux.h5')


# simList = [simDataName1, simDataName2, simDataName3,simDataName4]

simList = os.listdir(os.path.join('Analyses','noiseTest'))
# simDataName = os.path.join('Analyses', 'testFluxTest2.h5')
# simDataName2 = os.path.join('Analyses/LakeErieData','2024-02-28_12-00.h5')
# simDataName3 = os.path.join('Analyses/LakeErieData','2024-02-22_14-24.h5')

colourWheel = ['blue','orange','green','red','brown','purple','cyan']

for i,simDataName in enumerate(simList[:]):

    if 'k-1' not in simDataName:

        evals = 0

        simDataName = os.path.join('Analyses','noiseTest',simDataName)

        labData = DataManager(simDataName,
                    format="SimData",
                    iteration=0,
                    t_path= 'LabData/t',
                    z_path=  'LabData/z',
                    var_path= 'LabData/T',
                    relative_path=False,
            )

        rho_lab = labData.data

            ### Construct initial variables
        # clr = RT_lineClrs()
        I_list = []
        err_list = []
        t_list = []
        eval_list = []
        # pl.ion()


        for iteration in range(0,100):
                #### Simulated Data 
                try:
                    ErrorData = HDF5data(simDataName,  
                        format="SimData",
                        iteration=iteration  ,
                        var_path='rho'       
                        )

                    (t_sim,z_sim,rho_sim) = ErrorData.loadData()

                    ## Find error and plot 
                    errZ = np.trapz((rho_sim - rho_lab).T**2,z_sim)/(np.amax(z_sim) - np.amin(z_sim))
                    
                    err= (np.trapz(np.array(errZ),np.array(t_sim)))/(np.amax(t_sim) - np.amin(t_sim))
                    # (np.trapz(np.array(err_int),np.array(self.sim_time)))
                    timeData = HDF5data(simDataName,  
                        format="SimData",
                        iteration=iteration  ,
                        var_path='evals'       
                        )

                    (t_sim,z_sim,evaluations) = timeData.loadData()
                    evals+=evaluations

                    timeData = HDF5data(simDataName,  
                        format="SimData",
                        iteration=iteration  ,
                        var_path='epsilon'       
                        )

                    (t_sim,z_sim,eps) = timeData.loadData()

                except:
                    err = np.nan
                    # tComp=np.nan
                    
                
                if eps != 0 or iteration ==0:
                    I_list.append(iteration); err_list.append(err); eval_list.append(evals)


        I_list = np.asarray(I_list)
        
        # simDataName = simDataName[25:-13]
        # simDataName = simDataName.replace('_',' ')
        # simDataName = simDataName.replace('k0',r'$\kappa_o$'+'$= 10^')
        # simDataName = simDataName.replace('Ke',r' $\kappa_e$ $= 10^')
        # simDataName = simDataName + 


        err_list = np.asarray(err_list)
        err_list = err_list #/np.nanmax(err_list)

        if 'Noisy' in simDataName:
            ###diffusive
            ls = '--'

        else:
            ls = '-'
        if 'k1' in simDataName:
            if 'tanh' in simDataName:
                colour = 'm'
            elif 'Parab' in simDataName:
                colour='c'

        else:
            if 'tanh' in simDataName:
                colour = 'orange'
            elif 'Parab' in simDataName:
                colour='g'
           

        paths = simDataName.split('\\')
        plotLabel = paths[2][:-3]
        plotLabel = plotLabel.replace('_',' ')
        plotLabel = plotLabel.replace('k1',r'$\kappa_0=1 m^2/s$')
        plotLabel = plotLabel.replace('k-3',r'$\kappa_0=0.001 m^2/s$')
        if  '' in simDataName:
            simDataName = simDataName.replace('Artificial','')
            ax0.semilogy(eval_list,err_list,linestyle=ls,label=plotLabel,color = colour)
            # ax0.axhline((t_sim[1] - t_sim[0])**3,linestyle = ls,color=colour)
        # else:
            # ax1.semilogy(I_list,err_list,linestyle=ls,label=simDataName)
            
# ax.semilogy(I_list[0:1020],err_list2,'-')
# ax.semilogy(I_list[500:661],err_list3,'o-')

    ## Convergence Rate 
# p0 = np.polyfit(I_list[12:],np.log(err_list[12:]),1)

# ax.plot(x_eg,np.exp(p0[0]*x_eg + p0[1]),'k--')

## Labels
ax0.set_xlabel('Forward Evaluations',fontsize=18)
ax0.set_ylabel(r'Error $J = \int\int(T-T_m)^2 dzdt$',fontsize=18)
# ax0.set_ylabel(r'$J/J_0$')
# ax0.set_xticks([0,1000,2000,3000,4000])
# ax.set_ylabel('Error Total')
ax0.set_xlim([-10,510])
# ax0.set_ylim([1e-6,1.2])
# ax0.set_title(r'Artificial Data')

# Legend 
handles, labels = ax0.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax0.legend(by_label.values(), by_label.keys(),fontsize=16,loc='upper right',ncol=2)#, bbox_to_anchor=(1.57, 0.8))
# plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

## Labels
ax1.set_xlabel('Iteration number')
# ax1.set_ylabel(r'Normalized Error $J/J_0$')
# ax.set_ylabel('Error Total')
# ax.set_xlim([0,5001])
ax1.set_ylim([1e-4,10])
ax1.set_title(r'Field and Laboratory Data')

# Legend 
# handles, labels = ax1.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# ax1.legend(by_label.values(), by_label.keys(),fontsize=5,loc='upper left', bbox_to_anchor=(1, 1))
# # ax.legend(['Variable top and bottom','Normalized'])

# for iteration in range(0,500):
#         if smoothList[iteration] == True:
#                 ax.axvline(iteration,ls = '--')
# pl.legend(loc='upper left', bbox_to_anchor=(1, 1))
pl.show()