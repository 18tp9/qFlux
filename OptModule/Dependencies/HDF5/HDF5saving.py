### Save Data to HDF5

from cmath import isnan
import os
import h5py as hdf
import numpy as np

from ...Parameters import shared_parameters


class HDF5data():
    def __init__(self,filename = shared_parameters.simDataName):
        self.parameters = shared_parameters()
        self.filename   = filename
        self.datafile   = shared_parameters.labDataName
        self.dataPath   = (r'/Timeseries/Temperature/T_Chain')
        self.tArrPath   = (r'/Timeseries/scales/time')
        self.RawZ_path  = (r'/Timeseries/scales/z_raw')
        self.RawT_path  = (r'/Timeseries/Temperature/T_raw')
        self.outpath    = self.parameters.outputPath
        self.inpath     = self.parameters.inputPath
        self.groupnames = []
        self.precision  = 1e-2
        self.sigdigs    = 2
        self.SimSize    = [self.parameters.Nz + 1,self.parameters.datasize]     #plus one because cheb mod
        
        self.z          = self.parameters.grid.z
        self.date       = self.parameters.labDate

        self.initialize_main()
        self.loadMeasured(self.parameters.labDate)

        
    
    def initialize_main(self):
        self.main           = hdf.File(self.filename,'w')
        
        initial = self.main.create_group('initial')

        forcing = np.zeros(self.SimSize)
        phis    = np.zeros(self.SimSize)

        initial['Fm1'] = forcing
        initial['phi'] = phis

        self.MeasuredData = {}

    
    def new_SimData(self,iter_num,type,data):
        '''
        eg. iter_num = 1 
        type = 'phi'/'rho'/'Fm1'
        '''
        
        groupname = 'iter_' + str(iter_num)
        try:
            new_group = self.main.create_group(groupname)
        except:
            new_group = self.main[groupname]
        new_group[type] = data

    def load_SimData(self,iter_num,type):
        if iter_num != None:
            groupname = 'iter_' + str(iter_num)
        else:
            groupname = 'initial'
        return self.main[groupname][type][...]


    def loadMeasured(self,dataset):
        with hdf.File(self.datafile) as file:
            self.MeasuredData[dataset] = file[self.dataPath][...]
        return self.MeasuredData[dataset][...]

    def load_TimeArray(self):
        with hdf.File(self.datafile) as file:
            self.data_TimeArray = file[self.tArrPath][...]
        return self.data_TimeArray

    def loadRawData(self):
        '''
        loads raw Data and raw z grid
        '''
        with hdf.File(self.datafile) as file:
            self.RawData    = file[self.RawT_path][...]
            self.RawZ       = file[self.RawZ_path][...]
            self.RawData    = np.asarray(self.RawData)
        '''for i in range(self.RawData.shape[0]):
            self.RawData[i,:] = np.delete(self.RawData[i,:], np.where(np.isnan(self.RawData)))'''
        return self.RawData


    def pullMatrix(self,group,datatype):
        grp = self.main[group]
        return grp[datatype][:,:]

    @staticmethod
    def data_index(data,index):
        '''
        index into measured/lab data
        '''
        #data = group[dataset][...]
        #self.curr_rho = data[:,index]
        return data[:,index]


    @staticmethod 
    def get_hData(ind_func,dataset,t,precision = 1, sigdigs = 2):
        '''
        t is index -- to be fixed
        '''
        #index = np.argmax(self.time==t)
        ### get data at time t
        #data = self.data_index(int(t/self.precision))


        '''
        Jason --- There is a bug here, It t= 0.001 , the try statement could be fine but the data will be from some other time  
        '''
        # try:
        #     ### get data at time t
        if (t/precision) % 1 == 0:
            data = ind_func(dataset,int(t/precision))

        # except:
        else:
             ###if no data point at time t, round
             roundedT    = round(t,sigdigs)
             index1      = int(roundedT/precision)

             ###make sure we rounded down
             if roundedT > t:
                 index1 = index1 - 1
             ###add 1 to get second index and average the two values
             index2 = index1 + 1
             if index1 < 10**(sigdigs):          #make sure index 2 is not out of bounds
                 data   = (ind_func(dataset,index1) + ind_func(dataset,index2))/2.0
             else:
                 data   = (ind_func(dataset,(index1 - 1)))
        return data
    
    def error_calc(self,rho):
        '''
        *** TO FILL 
        '''
        frho = (rho-self.curr_rho)**2
        err1 = np.trapz(frho,self.z)
        return err1
    
    def save_data(self,testName,datatype,data):
        try:
            self.SimData[testName][datatype] = data
        except:
            self.SimData[testName][datatype][...] = data

    def labData_npz(self):
        data = self.MeasuredData[self.date][...]
        npz = {'rho_Array':data,'z':self.z}
        filename = self.npzName
        np.savez(os.path.join(self.inpath,filename),**npz)
        


    
    


'''
testData = np.ones([5,5])

test.new_SimData('Test1','RhoData',testData)
test.get_groups()
print (test.pullMatrix('Test1','RhoData'))
#test.loadMeasured(test.dataPath,dataset='20201211')
index = (test.measured_index(dataset='20201211',index = 500))
print (index)'''