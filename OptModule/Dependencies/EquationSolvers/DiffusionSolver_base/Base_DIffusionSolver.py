import numpy as np 
import matplotlib.pyplot as pl
import time
from numba import int32,float32
from numba.experimental import jitclass



from .Operators import Operator 
from .Operators import CHEB
    

# spec = [
#     ('value', int32),               # a simple scalar field
#     ('array', float32[:]),          # an array field
# ]
# @jitclass(spec)
class timestepper(object):

    kappa0 = 1.0e-3           ## Base Diffusivity
    do_filter = True     ## Filter the solution

    def __init__(self,rho,operator,kappa0=1e-3):

        ### 
        # Save the background diffusivity 
        self.kappa0=kappa0

        ##### 
        # Save rho, F, t, timestepping coefficients 
        self.rho   = rho
        self.rhom1 = -1
        self.rhom2 = -1 

        self.F   = 0
        self.Fm1 = 0
        self.Fm2 = 0 

        self.t    = 0
        self.dt   = -1
        self.dtm1 = 0  
        self.dtm2 = 0

        self.alpha0, self.alpha1 = np.nan,np.nan
        self.alpha2, self.alpha3 = np.nan,np.nan
        self.beta1,  self.beta2  = np.nan,np.nan
        self.beta3               = np.nan

        ##### 
        # Save operator, which does the differencing
        self.op = operator

        ##### 
        # Set the number of iterations 
        self.jj = 0


    def update_relabel(self):
        '''
        Re-label dt, U, and rho
        '''
        self.dtm2  = self.dtm1;  self.dtm1  = self.dt
        self.rhom2 = self.rhom1; self.rhom1 = self.rho
        # self.Um2   = self.Um1;   self.Um1   = self.U
        self.Fm2   = self.Fm1;   self.Fm1   = self.F


    def update_timeStep(self,dt0,order=3):
        '''
        Get the time Step and update all of the integration coefficients
        - Input - self.dtm1,self.dtm2 - previous two timesteps
            - order = {1,2,3} - order of the time stepping
        '''	
        import random 	
        self.dt = dt0; #self.dt /= random.randint(1,4)
        # self.dt = max(min(self.dt,1.0/75),1e-6)	

        if order==1:
            # self.dt = max((self.dt**3)/4,1e-8);
            self.alpha0 = 1.0; self.alpha1 = -1.0;
            self.alpha2 = 0; self.alpha3 = 0; 
            
            self.beta0 = 1.0;
            self.beta1 = 0; self.beta2 = 0; 

        elif order==2:
            # self.dt = max((self.dt**(3))/4,1e-8);
            # self.dt = max((self.dt**(3))/4,1e-8);
            self.alpha0 = 1.0 + self.dt/(self.dt+self.dtm1); self.alpha1 = -1.0*(self.dt+self.dtm1)/(self.dtm1);
            self.alpha2 = (self.dt*self.dt)/(self.dtm1)/(self.dt + self.dtm1); self.alpha3 = 0;
            
            self.beta0 = (self.dt + self.dtm1)/(self.dtm1); 
            self.beta1 = -1.0*(self.dt)/(self.dtm1);
            self.beta2 = 0;

        elif order==3:
            self.alpha0 = (1.0 + self.dt/(self.dt + self.dtm1) + self.dt/(self.dt + self.dtm1 + self.dtm2))
            self.alpha1 = -1.0*(self.dt + self.dtm1)/self.dtm1*(self.dt+self.dtm1+self.dtm2)/(self.dtm1+self.dtm2)
            self.alpha2 = (self.dt*self.dt)/(self.dt+self.dtm1)/(self.dtm1)*(self.dt+self.dtm1+self.dtm2)/(self.dtm2)
            self.alpha3 = -1.0*(self.dt*self.dt)/(self.dt+self.dtm1+self.dtm2)/(self.dtm1+self.dtm2)*(self.dt+self.dtm1)/(self.dtm2)
            
            self.beta0 = (self.dt+self.dtm1)/(self.dtm1)*(self.dt+self.dtm1+self.dtm2)/(self.dtm1+self.dtm2)
            self.beta1 = (-self.dt)/(self.dtm1)*(self.dt+self.dtm1+self.dtm2)/(self.dtm2)
            self.beta2 = (self.dt)/(self.dtm1+self.dtm2)*(self.dt+self.dtm1)/(self.dtm2)

        else:
            print('Time Step ordering not implemented')
        
        return self.dt 

    def get_explicit_forcing(self,F):
        raise NotImplementedError


    def SOlVE_IMPLICIT(self,rhs):
        '''
        Invert the time operator -- constant diffusivity
        '''


        if self.op.zType == 'Cosine': 
            
            rhs = self.op.flipud_even(rhs)
       
            #Invert Step 
            Frho = np.fft.fftn(rhs)
            #Invert Step 

            if isinstance(self.kappa0,np.float64) or isinstance(self.kappa0,float):
                Frho /= (1 - self.dt*self.kappa0/self.alpha0*(-self.op.K**2 ))

            else:
                raise NotImplementedError("Kappa isnt a float")

            # Filter once initial passes have been performed
            if self.jj >= 3 and self.do_filter: 
                Frho = self.op.apply_spec_filter(Frho)

            rho = np.real(np.fft.ifftn(Frho))	
            rho = self.op.flipud_crop(rho)

        elif self.op.zType == 'Cheb': 
            Dz = self.op.Dz
            if isinstance(self.kappa0,np.float64) or isinstance(self.kappa0,float):
                
                D2_Implicit_Op = np.eye(Dz.shape[0]) - self.dt*self.kappa0/self.alpha0*np.dot(Dz,Dz) 
                
                # # Neuman Bottom BC
                # D2_Implicit_Op[0,:] = Dz[0,:]

                # Dirichlet Bottom BC
                D2_Implicit_Op[0,:] = 0
                D2_Implicit_Op[0,0] = 1

                rhs[0] = self.op.BC1

                # # Neuman Top BC
                # D2_Implicit_Op[-1,:] = Dz[-1,:]

                # Dirichlet Top Bc
                D2_Implicit_Op[-1,:] = 0
                D2_Implicit_Op[-1,-1] = 1

                rhs[-1] = self.op.BC2

                rho = np.linalg.solve(D2_Implicit_Op,rhs)

            else:
                raise NotImplementedError("Kappa isnt a float")

        return rho

    def updateField(self,stepOrder):
        '''
        Time step the field rho
        * Solve Drho/self.dt = kappa Lap rho 
        Use third order Time stepping as in Subich (2013)
        Requires the use of two previous time steps to get convergence
        - Input - rho,U - density and velocity tuple at current time step
            - rhom1,Um1 - density at previous time step
            - rhom2,Um2 - density at previous previous time step 
            - stepOrder - order of timestep, stepOrder = {1,2,3}  
            - params - fluid parameters 
        *** Step ORDER ***
        1/(self.dt)*(sum alpha u^j) = L[u^(n+1)] + sum beta N[U^j]
            - L - Linear (diffusion) term
            - N - Nonlinear (advection) term

        '''
        #global self.alpha0 # Required for GMRES_FUNCTION 


        rhs = 1.0*self.rho[:]
        # First order method is simply Forward Euler
        if stepOrder == 1 :
            #print('First Order time-Stepping')
            #self.alpha0 = 1; self.alpha1 = -1
            #self.beta0 = 1;	
            
            F_uGrho = self.F 
            
            rhs *= -1*self.alpha1
            rhs -= self.beta0*self.dt*(F_uGrho)
        # Second Order	
        elif stepOrder == 2:
            #print('Second Order time-Stepping')
            #self.alpha0 = 3.0/2; self.alpha1 = -2; self.alpha2 = 0.5
            #self.beta0 = 2; self.beta1 = -1;		

            F_uGrho    = self.F
            F_uGrho_m1 = self.Fm1
            
            rhs *= -1*self.alpha1
            rhs -= self.beta0*self.dt*(F_uGrho)
            
            rhs -= self.alpha2*self.rhom1
            rhs -= self.beta1*self.dt*(F_uGrho_m1)
        # Third Order
        elif stepOrder == 3:
            #print('Third Order time-Stepping')
            #self.alpha0 = 11.0/6; self.alpha1 = -3; self.alpha2 = 1.5; self.alpha3 = -1.0/3
            #self.beta0 = 3; self.beta1 = -3; self.beta2 = 1;		


            F_uGrho    = self.F
            F_uGrho_m1 = self.Fm1
            F_uGrho_m2 = self.Fm2


            rhs *= -1*self.alpha1
            rhs -= self.beta0*self.dt*(F_uGrho)
            
            rhs -= self.alpha2*self.rhom1
            rhs -= self.beta1*self.dt*(F_uGrho_m1)
            
            rhs -= self.alpha3*self.rhom2
            rhs -= self.beta2*self.dt*(F_uGrho_m2)
        # Error
        else:
            print('Not impletmented order to time stepping'); 
        



        #Back solve for the density   
        #Invert Step 
        rho = self.SOlVE_IMPLICIT(rhs/self.alpha0)

        #Return 
        return rho




    def step(self,F,dt0):
        '''
        Take a time step
        Input --    F   - np.array  - Explicit forcing 
            # --  kappa   - float - CURRENTLY DOES NOT DO ANYTHING
            --  dt0     - float - previous timestep
        '''
        self.jj += 1
        # self.U  = U
        # self.kappa = kappa
        self.F = F


        ### Set the timestepping order 
        if self.jj == 1:
            stepOrder = 1
        elif self.jj == 2:
            stepOrder = 2
        else:
            stepOrder = 3
            
                
        _       = self.update_timeStep(dt0,order = stepOrder)
        rho1    = self.updateField(stepOrder = stepOrder)

        self.t = self.t + self.dt
        self.update_relabel()
        self.rho = rho1

        if np.any(np.isnan(rho1)):
                raise ValueError

        return self.t



