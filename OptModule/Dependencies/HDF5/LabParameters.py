###shared parameters
'''
class to share parameters such as grid, total time, data size, etc.
will use operator to create grid in this class and then pass it to phi, simloop etc.
that way if we change grid size, time - don't have to do it in each file
'''

import os
import numpy as np

from ..EquationSolvers.DiffusionSolver_base.Operators import Operator


class lab_parameters():


    #time parameters   
    rawZ           = np.array([0.003, 0.023, 0.043, 0.063, 0.083, 0.123, 0.143, 0.163, 0.183])

    def __init__(self,filename):


        if os.path.basename(filename) == 'LabData_20201211.h5':
            self.trunc_upperB    = 18000
            self.trunc_lowerB    = 16000

