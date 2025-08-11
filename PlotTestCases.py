import numpy as np 
# import h5py as h5 
import matplotlib.pyplot as pl 

import seaborn as sns
import os
from OptModule.Dependencies.Tools.Error_Analysis import error_analysis
from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data
from OptModule.Dependencies.Tools.Utilities import tools
from matplotlib.colors import LogNorm

sns.set_context('paper', font_scale=1.2)


simDataCheb = os.path.join('Analyses','time_extraDiff_ke-2+2e-3','FABLS_Diffusive.h5')
simDataParab = os.path.join('Analyses','perFlux_time_0423','FABLS_tanh_+_Per._Flux.h5')
simDataFlux = os.path.join('Analyses','time_randFlux_0402','FABLS_tanh_+_Rand._Flux.h5')

######################################
#####################################



### Test Cases
# tanh
tanhData = HDF5data(simDataCheb,            
        t_path=('LabData' + '/t'),
        z_path= ('LabData' + '/z'),
        var_path=('LabData' + '/T'),
        # format='LabData',
        )

(t_tanh,z_tanh,rho_tanh) = tanhData.loadData()
# rho_lab = rho_lab[:,:-1]

#### parab
parabData = HDF5data(simDataParab,            
        t_path=('LabData' + '/t'),
        z_path= ('LabData' + '/z'),
        var_path=('LabData' + '/T'),
        # format='LabData',
        )

(t_parab,z_parab,rho_parab) = parabData.loadData()

#### flux
fluxData = HDF5data(simDataFlux,            
        t_path=('LabData' + '/t'),
        z_path= ('LabData' + '/z'),
        var_path=('LabData' + '/T'),
        # format='LabData',
        )

(t_flux,z_flux,rho_flux) = fluxData.loadData()

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


fig = pl.figure(figsize=(9,7),constrained_layout=False)
axs = fig.add_gridspec(3,1)

ax0 = fig.add_subplot(axs[0])
ax1 = fig.add_subplot(axs[1])
ax2 = fig.add_subplot(axs[2])
fig.subplots_adjust(hspace=0.5)

# z_lab = z_lab[::-1]
# z = z[::-1]

#get time in seconds
# t_lab = (t_lab-np.min(t_lab))#*(24*60)*60

# #flip rhos
# rho_lab = np.flipud(rho_lab)
# rho_sim = np.flipud(rho_sim)

cax_l = ax0.contourf((t_tanh),z_tanh-0.5,rho_tanh,np.linspace(-1.2,1.2,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax0.set_title(r'a) Diffusive')

cax_s = ax1.contourf((t_parab),z_parab-0.5,rho_parab,np.linspace(-1.2,1.2,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax1.set_title(r'b) Periodic Flux')


cax_e = ax2.contourf((t_flux),z_flux-0.5,(rho_flux),np.linspace(-1.2,1.2,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax2.set_title(r'c) Random Flux')

# ax0.invert_yaxis()
# ax1.invert_yaxis()
# ax2.invert_yaxis()

ax2.set_xlabel('t (s)')
ax0.set_ylabel('z (m)')
ax1.set_ylabel('z (m)')
ax2.set_ylabel('z (m)')


from mpl_toolkits.axes_grid1 import make_axes_locatable

divider = make_axes_locatable(ax0)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
pl.colorbar(cax_l,cax=ax_cb,label = r'$^\circ C$',format='%.1f')

divider = make_axes_locatable(ax1)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
pl.colorbar(cax_s,cax=ax_cb,label = r'$^\circ C$',format='%.1f')

divider = make_axes_locatable(ax2)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
pl.colorbar(cax_e, cax=ax_cb,format='%.1f',label = r'$^\circ C$')


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