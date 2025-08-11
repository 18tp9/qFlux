import numpy as np 
# import h5py as h5 
import matplotlib.pyplot as pl 
from matplotlib import colors
import seaborn as sns
import os


import sys
from pathlib import Path

# Add parent directory to sys.path to avoid relative import error
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))

from OptModule.Dependencies.HDF5.HDF5_Generic import HDF5data
from OptModule.Dependencies.EquationSolvers.DiffusionSolver_base.Operators import Operator
from OptModule.Parameters import shared_parameters as params
from OptModule.Dependencies.Tools.Utilities import tools
# from OptModule.Dependencies.Tools.Dusolt_data import Dusolt
# from mpl_toolkits import mplot3d
from matplotlib import ticker
from logisticsManager import admin




sns.set_context('paper', font_scale=1.2)

# sW_I = Dusolt.getSW_I(0.2)*(1/4182)*(1/1000)*3600*24

# simDataName = os.path.join('Analyses','2024-08-09_11-29.h5')  #fake
# simDataName = os.path.join('Analyses','2024-08-09_17-36.h5')    #day204
# simDataName = os.path.join('Analyses','2024-08-12_09-10.h5')    #day204
simDataName = os.path.join('Analyses','2024-08-10_11-36.h5')    #day203
# simDataName = os.path.join('Analyses','2024-08-12_16-08.h5')
simDataName = os.path.join('Analyses','Erie2Best.h5')
simDataName = os.path.join('Analyses', 'IceOn+time_20201211.h5')
simDataName = os.path.join('Analyses','time_ChebBase','Automated_Step_tanh.h5')
simDataName = os.path.join('Analyses','Saddles_Erie2.h5')
simDataName = os.path.join('Analyses','time_middleFlux','Optimized_Step_tanh_+_flux.h5')
simDataName = 'Analyses\\time_perFlux_0402\\PABLO_tanh_+_Per._Flux.h5'
simDataName = 'Saddles_extraDiff2_ke-2+2e-3.h5'
######################################
#####################################
inputDataName = 'FakeData//perflux3_ke-2.h5'
inputData = admin.inData(inputDataName,var_path='Constant/trueF',t_path='Constant/t',z_path='Constant/z')

#### Simulated Data 
SimData = HDF5data(simDataName,  
        format="SimData",
        iteration=67 ,
        var_path='rho'       
        )

(t,z,rho) = SimData.loadData()

# for i in range(rho.shape[1]):
#     rho[:,i] = tools.sortProfile(rho[:,i])

# labGridInterp = HDF5data(simDataName,  
#         format="SimData",
#         iteration=10  ,
#         var_path='rho_labGrid'       
#         )

# (_,_,rho_sim) = labGridInterp.loadData()

CompData = HDF5data(simDataName,            
        t_path=('LabData' + '/t'),
        z_path= ('LabData' + '/z'),
        var_path=('LabData' + '/T'),
        # format='LabData',
        )

(t_lab,z_lab,rho_lab) = CompData.loadData()

fluxData = HDF5data(simDataName,  
        format="SimData",
        iteration=67 ,
        var_path='F'       
        )
(_,_,F) = fluxData.loadData()



# Hypso = Dusolt.spline(z)

# F = (F.T/Hypso).T
# F = F*(-1)

#subtract off SW radiation

# sW_I = tools.grid_change(z_lab,sW_I,z)



density = 999.842594 + 6.793952e-2*rho - (9.095290e-3)*rho**2 + (1.001685e-4)*rho**3 - (1.120083e-6)*rho**4 + (6.536332e-9)*rho**5
b = (9.81/999.842594)*(6.793952e-2*F - (9.095290e-3)*F**2 + (1.001685e-4)*F**3 - (1.120083e-6)*F**4 + (6.536332e-9)*F**5)

rhoGrad = np.asarray(np.gradient(rho,axis=0))
# maxRhoGrad = np.max(np.abs(rhoGrad),axis=0)

thermIndex = []
thermZ     = []

# for i in range(len(maxRhoGrad)):
#     thermLoc = np.where(np.abs(rhoGrad[:,i]) == maxRhoGrad[i])[0]
#     thermIndex.append(thermLoc[0])

#     thermoclineZ = z[np.where(np.abs(rhoGrad[:,i]) == maxRhoGrad[i])]
    
#     thermZ.append(thermoclineZ)

# thermZ = np.asarray(thermZ)


op = Operator(z[-1],len(z),'Cheb')
rho_z = op.get_grad_array(rho)

densGrad = op.get_grad_array((density))

# rho_z = tools.grid_change(z,rho_z,z_lab)
# rho_z = np.flipud(rho_z)

Nsquared = (((9.81/1000)*densGrad))#np.gradient(density,axis=0))

fig = pl.figure(figsize=(15,10),constrained_layout=False)
axs = fig.add_gridspec(3,1)

ax0 = fig.add_subplot(axs[0])
ax1 = fig.add_subplot(axs[1])
ax2 = fig.add_subplot(axs[2])
fig.subplots_adjust(hspace=0.5)

zP = z-0.5#[::-1]
# z_lab = z_lab[::-1]
# t = (t*6e-5)*24*3600

# pl.ion()
# pl.figure()
# for i in range(F.shape[0]):
#         Favg = np.mean(F)
#         Fz = F[i,:]
#         F_spectra = np.real(np.fft.fft(Fz-np.mean(Fz),axis=0,n=F.shape[1]))**2
#         fftfreqT = np.fft.fftfreq(F.shape[1],19)
#         pl.plot(fftfreqT[F_spectra.shape[0]//2:],F_spectra[F_spectra.shape[0]//2:])
#         pl.show()
#         pl.pause(0.1)
#         pl.cla()
        # fftfreqZ = np.fft.fftfreq(F.shape[0],zP[1]-zP[0])


# F_spectra[120:190,:] = 0



# osc1=F_spectra[:,:10]
# osc2=F_spectra[:,-10:]
# oscNew = np.zeros(F.shape)
# oscNew[:,:10] = F_spectra[:,:10]
# oscNew[:,-10:] = F_spectra[:,-10:]



# F_spectra[-4:-8] = 0
# F_filt = np.real(np.fft.ifft(oscNew))#,axis=0,n=F_spectra.shape[0]))
# F_filt = F_filt#+Favg

t_lab = (t_lab-np.min(t_lab))#*(24*60)*60

cax_l = ax0.contourf((t_lab),zP,(rho),np.linspace(-1.2,1.2,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax0.set_title(r'a) $T_{sim}$',fontsize=18)

###GLM
# cax_l = ax0.contourf((t_lab)/(3600*24) + 150,zP+2,(rho),np.linspace(np.min(rho),np.max(rho),50),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
# ax0.set_title('Sim Data')
# ax0.invert_yaxis()


# F = (F+(params.kappa0*rho_z))#/(t_lab[1]-t_lab[0])

cax_s = ax1.contourf((t_lab),zP,(F),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax1.set_title(r'b) $F$',fontsize=18)
# ax1.set_xlabel('t (mins)')

# for i in range(len(t_lab)):
#         ax1.axvline(t_lab[i]/60)

###GLM
# cax_s = ax1.contourf((t_lab)/(3600*24) + 150,zP+2,(-F),np.linspace(-1e-5,1e-4,30),extend='both',cmap=sns.color_palette("turbo", as_cmap=True))
# ax1.set_title('Flux F')

# cax_e = ax2.contourf((19*t_lab)/60,zP,np.flipud(F-F_filt),np.linspace(-0.03,3e-2,30),extend='both',cmap=sns.color_palette("turbo", as_cmap=True))
# ax2.set_title('F-oscillations')


# cax_s = ax1.contourf((t_lab),z+2,rho_z,np.linspace(-1,1,30),extend='both',cmap=sns.color_palette("rocket", as_cmap=True))
# ax1.set_title('rhoGrad -- Iteration 300  ')
# ax0.invert_yaxis()
# ax1.invert_yaxis()
# ax2.invert_yaxis()

# zero_x,zero_y = np.where(abs(rho_z) < 1e-1)
# rho_z[zero_x,zero_y] = None
# lowFx,lowFy = np.where(F<0)
# F[lowFx,lowFy] = 0

kappa = ((b/Nsquared))
kappa[(Nsquared)<1e-4] = 0

kappa = -F/rho_z
# kappa[F<3e-2] = 0
# kappa = kappa - sW_I


kappa  = np.nan_to_num(kappa)
# kappa = kappa+ params.kappa0


# ax2.axhline(1,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(5,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(9,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(13,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(19,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(25,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(31,linestyle='--',alpha=0.2,color='white')
# ax2.axhline(36,linestyle='--',alpha=0.2,color='white')

# cax_e = ax2.contourf((t),z,F_prime,np.linspace(-10,10,30),extend='both',cmap=sns.color_palette("turbo", as_cmap=True))
# ax2.set_title('Flux F - Shortwave')
# ax2.invert_yaxis()

cax_e = ax2.contourf((t_lab),zP,(kappa),np.linspace(1e-5,1,40),extend='both',cmap=sns.color_palette("turbo", as_cmap=True))
ax2.set_title(r'$\varkappa$',fontsize=18) #/\kappa_0$')
# ax2.invert_yaxis()

###GLM
# cax_e = ax2.contourf((t_lab)/(3600*24) + 150,zP+2,(kappa),np.linspace(-600,3000,50),extend='both',cmap=sns.color_palette("turbo", as_cmap=True))
# ax2.set_title(r'$\kappa_{e}/\kappa_0$')

error = (rho_lab - rho)


cax_e = ax2.contourf((t_lab),zP,(error),np.linspace(-1e-3,1e-3,40),extend='both',cmap=sns.color_palette("coolwarm", as_cmap=True))
ax2.set_title(r'c) $T-T_{sim}$',fontsize=18)
# ax2.invert_yaxis()


ax2.set_xlabel('t [s]',fontsize=16)
ax0.set_ylabel('z [m]',fontsize=16)
ax1.set_ylabel('z [m]',fontsize=16)
ax2.set_ylabel('z [m]',fontsize=16)
fig = pl.gcf()
fig.subplots_adjust(hspace=0.3)
# pl.tight_layout()

from mpl_toolkits.axes_grid1 import make_axes_locatable

divider = make_axes_locatable(ax0)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
cb0 = pl.colorbar(cax_l,cax=ax_cb,format='%.1f')
cb0.set_ticks([-1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2])
cb0.set_label(label=r'$C^\circ$',fontsize=16)

divider = make_axes_locatable(ax1)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
cb1 = pl.colorbar(cax_s,cax=ax_cb,format='%.1f')
cb1.set_ticks([-1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2])
cb1.set_label(label=r'$C^\circ m/s$',fontsize=16)

divider = make_axes_locatable(ax2)
ax_cb = divider.append_axes("right", size="2%", pad=0.05)
cb2 = pl.colorbar(cax_e,cax=ax_cb,format='%.0e')
cb2.set_ticks([-1e-3, -5e-4, 0.0, 5e-4, 1e-3])
cb2.set_label(label=r'$C^\circ$',fontsize=16)

pl.show()

# pl.clf()
# pl.ion()
# for i in range(rho.shape[1]):
#         pl.cla()
#         pl.plot(rho[:,i],z)
#         pl.pause(0.1)
#         pl.show()


# find depth averaged mixing rates over time

kappa_t = np.mean(kappa,axis=0)
F_t     = np.mean(F,axis=0)
rho_t   = np.mean(rho_lab,axis=0)

# pl.figure()
# # pl.plot(t,kappa_t)
# pl.plot(t,F_t)
# pl.legend('$\kappa$','F')

# pl.show()


# Calculate the 2D FFT
fft_2d = np.fft.fft2(F)

# Shift the zero-frequency component to the center of the spectrum
fft_2d_shifted = np.fft.fftshift(fft_2d)

# Calculate the magnitude spectrum
magnitude_spectrum = np.abs(fft_2d_shifted)

# Frequency vectors
freq_x = np.fft.fftfreq(F.shape[1])
freq_y = np.fft.fftfreq(F.shape[0],d=np.mean(np.diff(zP)))

# Shift the zero-frequency component to the center of the spectrum
freq_x_shifted = np.fft.fftshift(freq_x)
freq_y_shifted = np.fft.fftshift(freq_y)


# Create a filter (e.g., a low-pass filter)
rows, cols = F.shape
crow, ccol = rows // 2 , cols // 2  # center
mask = np.ones((rows, cols), np.uint8)
r = 3  # radius of the filter
center = [crow, ccol]
x, y = np.ogrid[:rows, :cols]
# mask_area = (x - center[0])**2 + (y - center[1])**2 <= r*r

mask[:,ccol-2:ccol+2] = 0
mask[crow-1:crow+1,:] = 0

# Apply the filter to the shifted FFT
filtered_fft_shifted = fft_2d_shifted * mask

# Shift the zero-frequency component back to the original position
filtered_fft = np.fft.ifftshift(filtered_fft_shifted)

# Compute the inverse 2D FFT
filtered_F = np.fft.ifft2(filtered_fft)

# Take the real part of the inverse FFT (since the result might be complex)
filtered_F = np.real(filtered_F)


# Plot the original array and its 2D FFT magnitude spectrum
fig, ([ax1, ax2], [ax3, ax4]) = pl.subplots(2, 2, figsize=(16, 12))
fig.subplots_adjust(wspace=0.25,hspace=0.22)

ctf1 = ax1.contourf(t,zP,F,np.linspace(-3,3,30),cmap='coolwarm',extend='both')
ax1.set_title('a)',fontsize=18)

ctf2 = ax2.contourf(np.log10(magnitude_spectrum), np.linspace(-1,5,50), extent=(-0.02, 0.02, -0.02, 0.02))
ax2.set_title('c)',fontsize=18)

ctf3 = ax3.contourf(t,zP,filtered_F,np.linspace(-3,3,30),cmap='coolwarm',extend='both')
ax3.set_title('b)',fontsize=18)

ctf4 = ax4.contourf(np.log10(np.abs(filtered_fft_shifted)), np.linspace(-1,5,50), extent=(-0.02, 0.02, -0.02, 0.02))
ax4.set_title('d)',fontsize=18)

ax1.set_xlabel('t (s)',fontsize=14)
ax1.set_ylabel('z (m)',fontsize=14)

ax3.set_xlabel('t (s)',fontsize=14)
ax3.set_ylabel('z (m)',fontsize=14)

ax2.set_xlabel('Temporal Frequency (Hz)',fontsize=14)
ax2.set_ylabel('Spatial Frequency (k)',fontsize=14)

ax4.set_xlabel('Temporal Frequency (Hz)',fontsize=14)
ax4.set_ylabel('Spatial Frequency (k)',fontsize=14)

cb1 = pl.colorbar(ctf1,format='%.1f')
cb1.set_label(label = r'$^\circ C$',fontsize=16)

cb3 = pl.colorbar(ctf3,format='%.1f')
cb3.set_label(label = r'$^\circ C$',fontsize=16)

cb2 = pl.colorbar(ctf2,format='%.1f')
cb2.set_label(label = r'$log_{10}$(PSD)',fontsize=16)

cb4 = pl.colorbar(ctf4,format='%.1f')
cb4.set_label(label = r'$log_{10}$(PSD)',fontsize=16)

pl.show()
