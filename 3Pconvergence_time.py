import numpy as np 
# import h5py as h5 
import matplotlib.pyplot as pl 
import seaborn as sns
import os

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data
from OptModule.Dependencies.Tools.Error_Analysis import error_analysis
from OptModule.Dependencies.EquationSolvers.Loop_Data import DataManager
from OptModule.Dependencies.Tools.Utilities import tools

# from Figures.LineStyles import RT_lineClrs
sns.set_context('paper', font_scale=1.5)

# fig,axs = pl.subplots(1,2,figsize=(10,4),constrained_layout=True)
# ax0 = axs[0]
# ax1 = axs[1]

fig,axs = pl.subplots(1,3,figsize=[16,8])

# ax.colorbar()
# simDataName = os.path.join('Analyses','IceOn_short_20201211.h5')
# simDataName = os.path.join('Analyses','Cosine_FluxTest2.h5')
# simDataName1 = os.path.join('Analyses','Cheb_k0-1_middleFlux.h5')
# simDataName2 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-3.h5')
# simDataName3 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-4.h5')
# simDataName4 = os.path.join('Analyses','Cheb_k0-1_BoundaryFlux.h5')


# simList = [simDataName1, simDataName2, simDataName3,simDataName4]

simList = os.listdir(os.path.join('Analyses','time_extraDiff_ke-2+2e-3'))
# simDataName = os.path.join('Analyses', 'testFluxTest2.h5')
# simDataName2 = os.path.join('Analyses/LakeErieData','2024-02-28_12-00.h5')
# simDataName3 = os.path.join('Analyses/LakeErieData','2024-02-22_14-24.h5')

colourWheel = ['blue','orange','green','red','brown','purple','cyan']
for j,folder in enumerate(['time_extraDiff_ke-2+2e-3','perFlux_time_0423','time_randFlux_0402']):
    for i,simDataName in enumerate(os.listdir(os.path.join('Analyses',folder))):

        evals = 0

        simDataName = os.path.join('Analyses',folder,simDataName)

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


        for iteration in range(0,2500):
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
                    
                    err= (np.trapz(np.array(errZ),np.array(t_sim)))#/(np.amax(t_sim) - np.amin(t_sim))
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
                    I_list.append(iteration); err_list.append(err); eval_list.append(evals+iteration)


        I_list = np.asarray(I_list)
        
        # simDataName = simDataName[25:-13]
        # simDataName = simDataName.replace('_',' ')
        # simDataName = simDataName.replace('k0',r'$\kappa_o$'+'$= 10^')
        # simDataName = simDataName.replace('Ke',r' $\kappa_e$ $= 10^')
        # simDataName = simDataName + 


        err_list = np.asarray(err_list)
        err_list = err_list #/np.nanmax(err_list)

        if 'Backtrack' in simDataName:
            ###diffusive
            ls = '-'
        elif ('FABLS' in simDataName) or ('GSS' in simDataName):
            ### artificial flux
            ls = '-'
        else:
            ls = '--'

        if '25%' in simDataName:
            colour = 'firebrick'
        elif '50' in simDataName:
            colour='seagreen'
        elif '75' in simDataName:
            colour = 'orange'
        elif '80' in simDataName:
            colour = 'cyan'
        elif 'GSS' in simDataName:
            colour='fuchsia'
        else:
            colour = 'steelblue'

        paths = simDataName.split('\\')
        plotLabel = paths[2][:-3]
        words = plotLabel.split('_')
        plotLabel = words[:]
        plotLabel = ' '.join(plotLabel)
        plotLabel = plotLabel.replace('$-v',r'$\v')
        plotLabel = plotLabel.replace('on ',r'on_')



        if  '' in simDataName:
            simDataName = simDataName.replace('Artificial','')
            axs[j].semilogy(eval_list,err_list,linestyle=ls,label=plotLabel,color = colour)
            # ax0.axhline((t_sim[1] - t_sim[0])**3,linestyle = ls,color=colour)
        # else:
            # ax1.semilogy(I_list,err_list,linestyle=ls,label=simDataName)
        
# ax.semilogy(I_list[0:1020],err_list2,'-')
# ax.semilogy(I_list[500:661],err_list3,'o-')

    ## Convergence Rate 
# p0 = np.polyfit(I_list[12:],np.log(err_list[12:]),1)

# ax.plot(x_eg,np.exp(p0[0]*x_eg + p0[1]),'k--')

## Labels
axs[0].set_xlabel('Evaluations',fontsize=18)
axs[1].set_xlabel('Evaluations',fontsize=18)
axs[2].set_xlabel('Evaluations',fontsize=18)
# axs[0].set_ylabel(r'Error $J = \int\int(T-T_m)^2 dzdt$',fontsize=18)
axs[0].set_ylabel(r'$J\ \mathrm{( ^\circ C^2)}$',fontsize=18)

# ax0.set_xticks([0,1000,2000,3000,4000])
# ax.set_ylabel('Error Total')
# axs[0].set_xlim([-1,8000])
axs[0].set_ylim([5e-7,1])
# axs[1].set_xlim([-1,8000])
axs[1].set_ylim([5e-7,1])
# axs[2].set_xlim([-1,7500])
axs[2].set_ylim([5e-7,1])

axs[0].text(0,1e-6,'a)',fontsize=18)
axs[1].text(0,1e-6,'b)',fontsize=18)
axs[2].text(0,1e-6,'c)',fontsize=18)


axs[0].set_title('Case A: Diffusive',fontsize=20)
axs[1].set_title('Case B: Periodic Flux',fontsize=20)
axs[2].set_title('Case C: Random Flux',fontsize=20)


# Legend 
handles, labels = axs[0].get_legend_handles_labels()
by_label = dict(zip([r'CSS @ 25% $\varepsilon_{med}$',r'CSS @ 50% $\varepsilon_{med}$','FABLS', 'GSS'], handles))
axs[0].legend(by_label.values(), by_label.keys(),fontsize=16,loc='upper right')#, bbox_to_anchor=(1, 1))
# plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

# ## Labels
# ax1.set_xlabel('Iteration number')
# # ax1.set_ylabel(r'Normalized Error $J/J_0$')
# # ax.set_ylabel('Error Total')
# # ax.set_xlim([0,5001])
# ax1.set_ylim([1e-4,2])
# ax1.set_title(r'Field and Laboratory Data')

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