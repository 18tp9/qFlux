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

fig,axs = pl.subplots(1,2,figsize=(10,6),constrained_layout=True)
ax0 = axs[0]
ax1 = axs[1]

fig = pl.figure()
ax0 = pl.gca()
ax1 = ax0.twinx()
# ax.colorbar()
# simDataName = os.path.join('Analyses','IceOn_short_20201211.h5')
# simDataName = os.path.join('Analyses','Cosine_FluxTest2.h5')
# simDataName1 = os.path.join('Analyses','Cheb_k0-1_middleFlux.h5')
# simDataName2 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-3.h5')
# simDataName3 = os.path.join('Analyses','Cheb_k0-1_ChebBase_ke-4.h5')
# simDataName4 = os.path.join('Analyses','Cheb_k0-1_BoundaryFlux.h5')


# simList = [simDataName1, simDataName2, simDataName3,simDataName4]

simList = os.listdir(os.path.join('Analyses','time_ChebBase'))
# simDataName = os.path.join('Analyses', 'testFluxTest2.h5')
# simDataName2 = os.path.join('Analyses/LakeErieData','2024-02-28_12-00.h5')
# simDataName3 = os.path.join('Analyses/LakeErieData','2024-02-22_14-24.h5')

colourWheel = ['blue','orange','green','red','brown','purple','cyan']

numIter = {}

for i,simDataName in enumerate(simList[:]):

    simDataName = os.path.join('Analyses','time_ChebBase',simDataName)

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

    if ''  in simDataName:

        for iteration in range(0,2000):
                #### Simulated Data 
                try:
                    ErrorData = HDF5data(simDataName,  
                        format="SimData",
                        iteration=iteration  ,
                        var_path='epsilon'       
                        )

                    (t_sim,z_sim,eps) = ErrorData.loadData()


                except:
                    eps = None




                    
            

                if eps != None:
                    I_list.append(iteration); eps_list.append(eps); 



        I_list = np.asarray(I_list)
        
        numIter[simDataName] = I_list[-1]

            # ax0.axhline((t_sim[1] - t_sim[0])**3,linestyle = ls,color=colour)


