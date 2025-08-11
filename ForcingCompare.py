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


simDataName1 = os.path.join('Analyses','time_middleFlux','Automated_Step_middleFlux.h5')
simDataName2 = os.path.join('Analyses','time_middleFlux','40%_Optimal_Step_middleFlux.h5')



######################################
#####################################


#### Saddles 
SimData0 = HDF5data(simDataName1,  
        format="SimData",
        iteration=99,
        var_path='F'       
        )

(t,z,Fsadd) = SimData0.loadData()

#### Median Step
SimData = HDF5data(simDataName2,  
        format="SimData",
        iteration=493,
        var_path='F'       
        )

(t,z,Fmed) = SimData.loadData()

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

Fdiff = Fsadd - Fmed

pl.contourf(t,z,Fdiff,np.linspace(-0.01,0.01,20))
pl.ylabel('z (m)')
pl.xlabel('t(s)')
pl.title('F difference: Optimizing vs. 0.4 * optimal')
cb = pl.colorbar(label=r'$F_{diff}$')


pl.show()

# fig.savefig("../Figures//Comparison_Constant_Noise_Case12.png",dpi=300)
# fig.savefig("../Figures//Comparison_Constant_Subsample4.png",dpi=300)