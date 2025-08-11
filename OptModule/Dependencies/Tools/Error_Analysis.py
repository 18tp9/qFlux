import numpy as np
import matplotlib.pyplot as pl
import os 


class error_analysis():

    

    def __init__(self,error_path):
        '''
        Class to hangle to calculation of the errors
        '''
        self.error_path = error_path


    @staticmethod
    def _create_error_files(error_path):
        ###clear error.txt file
        for path in ['error',]:
            with open(os.path.join(error_path,path+'.txt'),'a') as f:
                f.write(path + '\n')


    @staticmethod
    def write_error(fname,err):
        '''
        Function writes error value to txt file
        '''
        with open(fname,'a') as f:
            f.write(str(err)+'\n')


    @staticmethod
    def error_at_fixed_points(rho_sim,z_sim,rho_lab,rawZ):
        '''
        Define the error at certain fixed points
        --- rawZ is an array 
        '''


        ## If equal, then rho_lab may be nan
        if len(rho_lab)==len(rawZ):
            rawZ = rawZ[~np.isnan(rho_lab)]
        else:
            assert np.amax(np.isnan(rho_lab))==0, "lab data is nan"
            assert len(rho_sim)==len(rho_lab), "Dimension mismatch"

        ## Find sim data grid point closest to lab data
        # raw_location = [np.argmin(np.abs(z_sim - value)) for value in rawZ]

       
        rho_sim_lab = np.zeros(rawZ.shape)

        ### Weighted interpolation
        for i,zL0 in enumerate(rawZ):

            arg_upper = np.argmax(z_sim>zL0) if np.argmax(z_sim>zL0)>0 else -1
            z1        = z_sim[arg_upper] 
            z0        = z_sim[arg_upper-1] if arg_upper >0 else z_sim[0]
            
            w1 = (z1 - zL0)/(z1-z0) 
            w0 = (zL0 - z0)/(z1-z0) 
            rho_sim_lab[i] = rho_sim[arg_upper-1]*w1 + rho_sim[arg_upper]*w0
            

        #     interpolated = (rho_sim[raw_location[i]]*w2 + rho_sim[raw_location[i]+1]*w1)

        #     z = z_sim[raw_location[i]]  #sim grid point nearest to ith lab grid point

        #     if (z < rawZ[i]) and ((i+1) < (len(raw_location))): #closest point is less than labZ

        #         zp1 = z_sim[raw_location[i+1]]  #opposite side of closest point

        #         #weights for interpolation
        #         w1 = abs(rawZ[i] - z)/abs(z-zp1)
        #         w2 = abs(rawZ[i] - zp1)/abs(z-zp1)

        #         interpolated = (rho_sim[raw_location[i]]*w2 + rho_sim[raw_location[i]+1]*w1)
        #         rho_sim_raw.append(interpolated)


        #     elif (z > rawZ[i]):     # closest point is less than labZ

        #         zm1 = z_sim[raw_location[i-1]]

        #         w1 = abs(rawZ[i] - z)/abs(z-zm1)
        #         w2 = abs(rawZ[i] - zm1)/abs(z-zm1)

        #         interpolated = (rho_sim[raw_location[i]]*w1 + rho_sim[raw_location[i]-1]*w2)
        #         rho_sim_raw.append(interpolated)

            
        #     else:   #just append closest point
        #         rho_sim_raw.append(rho_sim[raw_location[i]])


        # print(rho_sim_raw)
        # print(rho_sim)
        # pl.ion(); pl.cla()
        # pl.plot(rho_sim_lab,rawZ,'o')
        # # pl.plot(z_sim[raw_location],rho_sim[raw_location],'x')
        # pl.plot(rho_lab,rawZ,'^')
        # pl.xlabel('Temp C')
        # pl.ylabel('z(m)')
        # # pl.legend('Interpolated','Simulated')
        # # pl.savefig(f'')
        # # pl.show()
        # # pl.pause(0.1)


        # # print(rho_sim_raw)
        # # print(rho_sim[raw_location])

        # rho_sim_raw = np.asarray(rho_sim_raw)



        # # rho_lab_raw = rho_lab[raw_location] if len(z_sim)==len(rawZ) else rho_lab[~np.isnan(rho_lab)]
        

        # n = len(rawZ)
        # error = np.trapz((rho_sim_lab - rho_lab)**2,rawZ)
        error = rho_sim_lab - rho_lab
        return error
        
    @staticmethod
    def error_integrated(rho_sim,z_sim,rho_lab,rawZ) -> float:
        '''
        Define the error, integrated over z
        '''

        rho_sim_lab = np.zeros(rawZ.shape)

        ### Weighted interpolation
        for i,zL0 in enumerate(rawZ):

            arg_upper = np.argmax(z_sim>zL0) if np.argmax(z_sim>zL0)>0 else -1
            z1        = z_sim[arg_upper] 
            z0        = z_sim[arg_upper-1] if arg_upper >0 else z_sim[0]
            
            w1 = (z1 - zL0)/(z1-z0) 
            w0 = (zL0 - z0)/(z1-z0) 
            rho_sim_lab[i] = rho_sim[arg_upper-1]*w1 + rho_sim[arg_upper]*w0
        # ## Subsample if not the same size 
        # if (len(z_sim) != len(z_lab)) or (np.amax(z_sim - z_lab)>1e-6):
        #     #print(z_lab)
        #     #print(z_sim)
            
        #     z_closest = np.abs(z_sim - value)).argmin() for value in (z_lab[~np.isnan(rho_lab)
        #     print(z_closest.shape)
        #     print(rho_sim.shape)
        #     rho_sim = rho_sim[z_closest]
        
        error = np.trapz((rho_sim_lab - rho_lab)**2,rawZ)/(np.amax(rawZ) - np.amin(rawZ)) 
        return error