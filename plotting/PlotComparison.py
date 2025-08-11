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

from OptModule.Dependencies.Tools.Error_Analysis import error_analysis
from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data
from OptModule.Dependencies.Tools.Utilities import tools
from matplotlib.colors import LogNorm

sns.set_context('paper', font_scale=1.2)


# surfFlux = Dusolt.heatFlux()





simDataName = os.path.join('Analyses','2024-08-09_11-29.h5')  #fake
# simDataName = os.path.join('Analyses','2024-08-09_17-36.h5')    #day204
# simDataName = os.path.join('Analyses','2024-08-10_11-36.h5')    #day203
# simDataName = os.path.join('Analyses','2024-08-12_14-07.h5')    #day204
# simDataName = os.path.join('Analyses','2024-08-12_16-08.h5')


# simDataName = os.path.join('Analyses','Erie2Best.h5')
simDataName = os.path.join('Analyses','Cheb_k0-6_ChebBase_ke-4.h5')
# simDataName = os.path.join('Analyses', 'testFluxTest2.h5')
simDataName = os.path.join('Analyses', 'Cheb_k0-2_IceOn_short_20201211.h5')
simDataName = os.path.join('Analyses','Saddles_Erie2.h5')
# simDataName = os.path.join('Analyses','time_middleFlux','Automated_Step_tanh_+_flux.h5')
# simDataName = os.path.join('Analyses','time_middleFlux','Automated_Step_middleFlux.h5')
simDataName = os.path.join('Analyses','time_tanhFlux','Optimized_Step_tanh_+_flux.h5')
simDataName = 'Analyses\\LakeData\\k1e-3Saddles_BIWA1TempInterp.h5'
simDataName = "Analyses\\time_randFlux2\\PABLO_tanh_+_Rand._Flux.h5"
simDataName = 'Analyses\\perFlux_time_0423\\PABLO_tanh_+_Per._Flux.h5'


# 2024-08-16_09-09




######################################
#####################################



### Test Case 
# caseName = "Constant_Noise_Case12"
CompData = HDF5data(simDataName,            
        t_path=('LabData' + '/t'),
        z_path= ('LabData' + '/z'),
        var_path=('LabData' + '/T'),
        # format='LabData',
        )

(t_lab,z_lab,rho_lab) = CompData.loadData()
# rho_lab = rho_lab[:,:-1]

#### iter 0 
SimData0 = HDF5data(simDataName,  
        format="SimData",
        iteration=0,
        var_path='rho'       
        )

(t,z,rho_sim0) = SimData0.loadData()

#### optimal 
SimData = HDF5data(simDataName,  
        format="SimData",
        iteration=83,
        var_path='rho'       
        )

(t,z,rho_sim) = SimData.loadData()

#### Simulated Data 
# SimDataRaw = HDF5data(simDataName,  
#         format="SimData",
#         iteration=99,
#         var_path='rho_labGrid'       
#         )

# (_,_,rho_sim_raw) = SimDataRaw.loadData()
# t = t[:-1]
# rho_sim = rho_sim[:,:-1]

# rho_lab_sim = tools.grid_change(z_lab,rho_lab,z)


# H = (surfFlux)*(1/(1000*4182)) #+ np.median(z)


fig = pl.figure(figsize=(15,10),constrained_layout=False)
axs = fig.add_gridspec(3,1)

ax0 = fig.add_subplot(axs[0])
ax1 = fig.add_subplot(axs[1])
ax2 = fig.add_subplot(axs[2])
fig.subplots_adjust(hspace=0.5)

z_lab = z_lab[::-1]
# z = z[::-1]

#get time in seconds
t_lab = (t_lab-np.min(t_lab))#*(24*60)*60

# #flip rhos
# rho_lab = np.flipud(rho_lab)
# rho_sim = np.flipud(rho_sim)


cax_l = ax0.contourf((t_lab),z-0.5,rho_lab,np.linspace(-1.05,1.05,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax0.set_title(r'$T$',fontsize=18)

# cax_s = ax1.contourf((t_lab),z,rho_sim,np.linspace(-1,1,50),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
# ax1.set_title('Sim Data')

error = (rho_lab - rho_sim)
error0 = rho_lab - rho_sim0

cax_s = ax1.contourf((t_lab),z-0.5,rho_sim,np.linspace(-1.05,1.05,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax1.set_title(r'$T_{sim}$',fontsize=18)


cax_e = ax2.contourf((t_lab),z-0.5,(error),np.linspace(-0.005,0.005,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax2.set_title(r'$T - T_{sim}$',fontsize=18)

# ax0.invert_yaxis()
# ax1.invert_yaxis()
# ax2.invert_yaxis()

ax2.set_xlabel('t [s]',fontsize=16)
ax0.set_ylabel('z [m]',fontsize=16)
ax1.set_ylabel('z [m]',fontsize=16)
ax2.set_ylabel('z [m]',fontsize=16)
fig.subplots_adjust(hspace=0.5)


from mpl_toolkits.axes_grid1 import make_axes_locatable

divider = make_axes_locatable(ax0)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
cb0 = pl.colorbar(cax_l,cax=ax_cb,format='%.1f')
cb0.set_label(label = r'$^\circ C$',fontsize=16)

divider = make_axes_locatable(ax1)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
cb1 = pl.colorbar(cax_s,cax=ax_cb,format='%.1f')
cb1.set_label(label = r'$^\circ C$',fontsize=16)

divider = make_axes_locatable(ax2)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
cb2 = pl.colorbar(cax_e, cax=ax_cb,format='%.0e')
cb2.set_label(label = r'$^\circ C$',fontsize=16)


# pl.clf()
# pl.ion()
# ax = pl.gca()
# for i in range(rho_lab.shape[1]):
#         pl.cla()
#         pl.plot(rho_lab[:,i],z)
#         pl.plot(rho_sim[:,i],z)
#         pl.xlabel('Temp (\u2103)')
#         pl.ylabel('z (m)')
#         pl.xlim([9,23])
#         pl.legend(['Field Data','Simulated Data'])
#         ax.invert_yaxis()
#         # pl.pause(0.01)
#         pl.savefig(f'Figures/profilePNGs/{i:04}.png')
        # pl.show()

pl.show()

# fig.savefig("../Figures//Comparison_Constant_Noise_Case12.png",dpi=300)
# fig.savefig("../Figures//Comparison_Constant_Subsample4.png",dpi=300)