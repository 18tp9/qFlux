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

# fig,axs = pl.subplots(1,2,figsize=(10,6),constrained_layout=True)
# ax0 = axs[0]
# ax1 = axs[1]

fig,axs = pl.subplots(1,2,figsize=(12,7))
ax0 = axs[0]
ax1 = ax0.twinx()
# ax.colorbar()
# simDataName = os.path.join('Analyses','IceOn_short_20201211.h5')
# simDataName = os.path.join('Analyses','Cosine_FluxTest2.h5')
# simDataName1 = os.path.join('Analyses','Cheb_k0-1_middleFlux.h5')
# simDataName2 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-3.h5')
# simDataName3 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-4.h5')
# simDataName4 = os.path.join('Analyses','Cheb_k0-1_BoundaryFlux.h5')


# simList = [simDataName1, simDataName2, simDataName3,simDataName4]

simList = os.listdir(os.path.join('Analyses','HYBRID'))

# simDataName = os.path.join('Analyses', 'testFluxTest2.h5')
# simDataName2 = os.path.join('Analyses/LakeErieData','2024-02-28_12-00.h5')
# simDataName3 = os.path.join('Analyses/LakeErieData','2024-02-22_14-24.h5')

colourWheel = ['blue','firebrick','green','black','magenta','darkorange','lime','purple','deeppink','chartreuse','seagreen','darkorange','pink','turquoise']

# for i,simDataName in enumerate(simList[:1]):
#   simDataName = os.path.join('Analyses','HYBRID',simDataName)
# simDataName = os.listdir(os.path.join('Analyses','Saddles_subIter_middleFlux.h5'))
simDataName = os.path.join('Analyses','plotting_Parab','Saddles_ParabRho_ke-2.h5')
simDataName = 'Saddles_p_perflux2_ke-2.h5'



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
eps_list = []
err_list = []
gammList = []
smoothList = []

# pl.ion()
###parab
# iterList = [1,4,8,10,11,12,20]#,25,26,27]

###tanh
# iterList = [1,3,4,5,7,10,12]

# ###tanh+F
iterList = [2,4,6,12,14,19]
# iterList = [2,22,45,63,71,80]


if ''  in simDataName:

    for iteration in range(0,21):
            #### Simulated Data 
            try:
                ErrorData = HDF5data(simDataName,  
                    format="SimData",
                    iteration=iteration  ,
                    var_path='epsilon'       
                    )

                (t_sim,z_sim,eps) = ErrorData.loadData()


            except:
                eps = np.nan

            try:
                ErrorData = HDF5data(simDataName,  
                    format="SimData",
                    iteration=iteration  ,
                    var_path='gamma'       
                    )

                (t_sim,z_sim,gamma) = ErrorData.loadData()


            except:
                gamma = np.nan

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
    
            except:
                err = np.nan
                
        

                
            if eps != -2:
                I_list.append(iteration); eps_list.append(eps); gammList.append(gamma) 
                err_list.append(err)


    I_list = np.asarray(I_list)
    
    # simDataName = simDataName[18:-3]
    # simDataName = simDataName.replace('_',' ')
    # simDataName = simDataName.replace('k0',r'$\kappa_o$'+'$= 10^')
    # simDataName = simDataName.replace('Ke',r' $\kappa_e$ $= 10^')
    # simDataName = simDataName + '$'


    eps_list = np.asarray(eps_list)
    gammList = np.asarray(gammList)
    # err_list = err_list#/np.nanmax(err_list)

    if 'Backtrack' in simDataName:
        ###diffusive
        ls = '-'
    elif 'Saddles' in simDataName:
        ### artificial flux
        ls = '--'
    else:
        ls = '-'

    # if 'ke-2' in simDataName:
    #     colour = 'red'
    # elif 'ke-3' in simDataName:
    #     colour='green'
    # else:
    #     colour = 'blue'
    

    if  ''  in simDataName:
        simDataName = simDataName.replace('Artificial','')
        ax1.plot(I_list,eps_list*gammList,'-o',color='red',alpha=0.3,label=r'Step Size $\epsilon\gamma$')
        ax0.plot(I_list,err_list,label=r'Error J',color='tab:blue')
        
        for i,val in enumerate(iterList):
                
                ax1.plot(val-1,eps_list[val-1]*gammList[val-1],'o',markersize=6,color=colourWheel[i])
        
        # ax0.axhline((t_sim[1] - t_sim[0])**3,linestyle = ls,color=colour)


## Labels


# ax.set_ylabel('Error Total')
# ax0.set_xlim([-1,23])
# ax0.set_ylim([1e-2,0.15])
# ax0.set_title(r'contraction 0.6')


for i,iteration in enumerate(iterList):
                #### Simulated Data 

                ErrorData = HDF5data(simDataName,  
                    format="SimData",
                    iteration=iteration -1 ,
                    var_path='alphas'       
                    )

                (t_sim,z_sim,alpha) = ErrorData.loadData()



                ErrorData = HDF5data(simDataName,  
                    format="SimData",
                    iteration=iteration  -1,
                    var_path='Jvals'       
                    )

                (t_sim,z_sim,Jvals) = ErrorData.loadData()

                axs[1].plot(alpha,Jvals,color=colourWheel[i])
                
                minLoc = np.where(Jvals == np.min(Jvals))
                axs[1].plot(alpha[minLoc],Jvals[minLoc],'o',markersize=5,color=colourWheel[i])

                
                axs[1].set_xlabel('Step Size ' + r'$\varepsilon\gamma$',fontsize=18)
                axs[1].set_xlim([-0.1,1])
                axs[1].set_ylim([0.01,0.125])
                axs[0].set_ylim([0.01,0.125])


                # Move y-axis ticks to the right
                axs[1].yaxis.tick_right()
                

                # Optionally, you can also move the y-axis label to the right
                axs[1].yaxis.set_label_position("right")
                axs[1].set_ylabel(r'$J \ \mathrm{(^\circ C^2)}$',fontsize=18)
                # pl.ylim([1e3,2e5])
                axs[1].text(0.92, 0.07, 'b)', transform=axs[1].transAxes,
                        verticalalignment='top', horizontalalignment='left',
                        fontsize=18, color='k')

fig.subplots_adjust(wspace=0.3)

ax1.text(0.92, 0.07, 'a)', transform=axs[0].transAxes,
                        verticalalignment='top', horizontalalignment='left',
                        fontsize=18, color='k')

ax1.set_ylabel(r'Step Size ' + r'$\varepsilon\gamma$',color='tab:red',fontsize=18)
ax1.tick_params(axis='y', labelcolor='tab:red')
ax0.set_xlabel('Iteration Number',fontsize=18)
ax0.set_ylabel(r'$J \ \mathrm{(^\circ C^2)}$',color='tab:blue',fontsize=18)
ax0.tick_params(axis='y', labelcolor='tab:blue')
 



                    
            

                    
