import numpy as np 
import matplotlib.pyplot as pl
import scipy.sparse.linalg as linal
from scipy.interpolate import interp1d
import time

from scipy.special import erf


class CHEB(object):
    @staticmethod
    def clencurt_weights(N):
        theta = np.pi*np.arange(0,N+1)/N
        x = np.cos(theta)
        w = np.zeros(N+1)
        ii = np.arange(1,N)
        v = np.ones(N-1)
        if np.mod(N,2)==0:
            w[0] = 1./(N**2-1)
            w[N] = w[0]
            for k in np.arange(1,int(N/2.)):
                v = v-2*np.cos(2*k*theta[ii])/(4*k**2-1)
            v = v - np.cos(N*theta[ii])/(N**2-1)
        else:
            w[0] = 1./N**2
            w[N] = w[0]
            for k in np.arange(1,int((N-1)/2.)+1):
                v = v-2*np.cos(2*k*theta[ii])/(4*k**2-1)
        w[ii] = 2.0*v/N
        return x, w


    @classmethod
    def int_cheb(cls,f,Nz = 1,Lz=1,axis=0):
        ''' 
        Integrate using a set of computed weights for a cheb grid
        Only an approximation
        - Specify the axis along which to integate. Default is axis=0
        !!!!! only the axis = 0 is implemented
        '''
        z0,dz_weights = cls.clencurt_weights(Nz-1);
        dz_weights *= Lz/2
        if axis>0:
            print('Not implemented')
        #
        integrand = np.dot(np.diag(dz_weights),f)
        return np.sum(integrand,axis=0);


    @staticmethod
    def cheb(N,xmin=-1,xmax=1):
        if (N == 0):
            D = 0
            x = 1
        x = np.cos(np.pi*np.arange(0,N+1)/N)
        c = np.hstack([2, np.ones(N-1), 2])*(-1)**np.arange(0,N+1)
        X = np.tile(x,(N+1,1))
        dX = X.T - X
        D = (c[:,np.newaxis]*(1.0/c)[np.newaxis,:])/(dX+(np.identity(N+1)))       # off-diagonal entries
        D = D - np.diag(D.sum(axis=1))              # diagonal entries

        # Shift and scale coordinates as appropriate
        xr = -(x -1 -xmin)/2 * (xmax-xmin)
        D /= -(xmax-xmin)/2 
        return D, xr

    @classmethod
    def cheb_diff(cls,var,N,Lz=2):
        '''
        Take the deriavtive of var
        Derivative must be along axis=0
        '''
        D,z = cls.cheb(N-1)
        D /= Lz/2
        return np.dot(D,var)



class Operator(object):


    def __init__(self,Lz,Nz,zType = 'Periodic',offset = None,BC1=None,BC2=None,gridMin = 0):
        '''
        Lz,Nz -- float -- domain parameters
        zType -- str   -- Type of domain boundaries - 'Periodic','Cosine','Cheb'
        offset -- float -- shift of z from 0
        BC1,BC2 -- float -- Only implemented on 'Cheb', n={0,-1} boundary conditions respectively. 
        '''

        self.Lz = Lz
        self.Nz = Nz
        self.zType =  zType

        if (zType == 'Periodic') :
            self.z = np.linspace(gridMin,Lz,Nz,endpoint=False)
            dz          = Lz/(Nz)
            self.k      = self.get_wavenumbers(Nz,dz,zType)
            self.K      = np.abs(self.k)	
            self.K_filt = np.abs(self.k/np.amax(self.k))	
    	
        elif (zType == 'Cosine'):
            self.z = np.linspace(gridMin,Lz,Nz,endpoint=True)
            dz          = Lz/(Nz-1)
            self.k      = self.get_wavenumbers(Nz,dz,zType)
            self.K      = np.abs(self.k)	
            self.K_filt = np.abs(self.k/np.amax(self.k))	

        elif zType == 'Cheb':
            self.Dz,self.z = CHEB.cheb(Nz-1,gridMin,Lz)
            self.k = self.Dz
            self.BC1 = BC1
            self.BC2 = BC2

        if not offset is None:
            self.z = self.__offset(self.z,offset)

    def __offset(self,z,offset):
        '''
        Shift the z by a given amount
        '''
        return z - offset

    @classmethod
    def get_grad1D(cls,u,k,dim = 'z',zType='Periodic',order =1):
        ''' 
        Find the gradient using a FFT
        - Input -- u - field to be differentiated
            - k - vector of wavenumber (in fft order)
            - index - String - axis on which to differentiate
                - options ('x','y','z') -> Dim (2,1,0)
            - zType - String - continuation type in vertical 
                - Option = Periodic or Cosine, Default = 'Periodic
            - order - order of derivative
        - Output - u_x - Array of size(u) which has been differentiated
        '''
        if dim == 'x':
            axis = 2
            Fu = np.fft.fft(u,axis = axis)
            Fu_x = (np.sqrt(-1 +0j))**order*np.multiply(Fu,k**order)
            u_x = np.real(np.fft.ifft(Fu_x,axis = axis,n=u.shape[axis]))
        elif dim == 'y':
            axis = 1
            Fu = np.fft.fft(u,axis = axis)
            Fu_x = (np.sqrt(-1 +0j))**order*np.multiply(Fu,k**order)
            u_x = np.real(np.fft.ifft(Fu_x,axis = axis,n=u.shape[axis]))
        elif dim == 'z':
            axis = 0
            if zType == 'Periodic':
                Fu = np.fft.fft(u,axis = axis)
                Fu_x = (np.sqrt(-1 +0j))**order*np.multiply(Fu,k**order)
                u_x = np.real(np.fft.ifft(Fu_x,axis = axis,n=u.shape[axis]))
            elif zType == 'Cosine':
                '''
                I can't figure out how to use one of the implemented DCTs to back out a derivative. For not I'm goin to flip perform the derivative that way. not pretty but it'll do the job for now.
                '''
                '''
                Fu = sp.dct(u,axis = axis,type=2)
                Fu_x = np.multiply(-k*Fu)
                #Fu_x = Fu
                #Fu_x = np.multiply(-Fu,k/np.cos(k+k[1,0,0]/2))
                #Fu_x[:-1,:,:] = Fu_x[1:,:,:]
                #Fu_x[u.shape[axis]-1,:,:] = 0
                
                u_x = np.real(sp.idst(Fu_x,axis = axis,n=u.shape[axis],type=2,norm='ortho'))
                '''
                temp = cls.flipud_even(u)
                Fu = np.fft.fft(temp,axis = axis)
                # Fu_x = (np.sqrt(-1 +0j))**order*np.dot(Fu.T,(k**order)).T
                Fu_x = (np.sqrt(-1 +0j))**order*np.multiply(Fu,(k**order))
                # Fu_x = (np.sqrt(-1 +0j))**order*np.apply_along_axis(lambda x:x*(k**order),arr=Fu,axis=0)
                u_x = np.real(np.fft.ifft(Fu_x,axis = axis,n=temp.shape[axis]))
                u_x = cls.flipud_crop(u_x)			
                
            elif zType == 'Sine':
                '''
                Sine differentiation - note cosine notice.
                '''
                temp = cls.flipud_odd(u)
                Fu = np.fft.fft(temp,axis = axis)
                Fu_x = (np.sqrt(-1 +0j))**order*np.multiply(Fu,k**order)
                u_x = np.real(np.fft.ifft(Fu_x,axis = axis,n=u.shape[axis]*2))
                u_x = cls.flipud_crop(u_x)
                
            elif zType == 'Cheb':
                Dz = k                  ### Hack to pass Dz matrix 
                u_x = np.dot(Dz,u)
                for _ in range(1,order):
                    u_x = np.dot(Dz,u_x)
                    
            else:
                print('Dim in FFT is not the correct type')

        return u_x


    @staticmethod
    def flipud_even(A,dim=0):
        ''' 
        flipud assuming even symmetry 
        -Input A - assumed to be a three dimensional array 
            - input - dimension on which to flipud
        -Ouput A extended evenly in dim dir
            size(A) = 2*n_dim -1 
        '''
        if dim == 0:
            #temp = np.append(A,np.flipud(A)[:,:,:],axis=0)   # double domain
            temp = np.append(A,np.flipud(A)[1:-1,],axis=0) # Trick flip 
        else:
            print('Not Implemented')

        return temp


    @staticmethod
    def flipud_odd(A,dim=0):
        ''' 
        flipud assuming odd symmetry 
        -Input A - assumed to be a three dimensional array 
            - input - dimension on which to flipud
        -Ouput A extended evenly in dim dir
            size(A) = 2*n_dim -1 
        '''
        if dim == 0:
            #temp = np.append(A,0*A[1:2,:,:],axis=0)
            #temp = np.append(temp,-1*np.flipud(A)[1:-1,:,:],axis=0)
            temp = np.append(A,-1*np.flipud(A)[1:-1,:,:],axis=0)
        else:
            print('Not Implemented')

        return temp

    @staticmethod
    def flipud_crop(A,dim=0):
        ''' 
        flipud assuming even symmetry 
        -Input A - assumed to be a three dimensional array 
            - input - dimension on which to flipud
        -Ouput A extended evenly in dim dir
            size(A) = 2*n_dim -1 
        '''
        if dim == 0:
            N = int(A.shape[0]/2)
            if A.ndim == 1:
                return A[:N+1]
            elif A.ndim==3:
                return A[:N+1,:,:]
        else:
            print('Not Implemented')

        # return temp

    #def interp_rho(self,u)

    def get_grad(self,u):
        ''' 
        Find the gradient using a FFT
        - Input -- u - field to be differentiated
            - k,l,m - vector of wavenumber
        - Output - (u_x,u_y,u_z) 
        '''

        u_z = self.get_grad1D(u,self.k,dim = 'z',zType=self.zType,order =1)

        return u_z

    def get_grad_array(self,A):
        '''
        Compute the gradient of A along is axis=0
        '''
        return np.apply_along_axis(self.get_grad,0,A)

    def get_lap(self,u):
        ''' 
        Find the laplacian using a FFT
        - Input -- u - field to be differentiated
            - k,l,m - vector of wavenumber
        - Output - u_xx + u_yy + u_zz - array of size(u) 
        '''
        #lap = get_grad1D(u,k,'x',order=2)
        #lap += get_grad1D(u,l,'y',order=2)
        #lap += get_grad1D(u,m,'z',zType,order=2)

        u_z = self.get_grad(u)
        lap  = self.get_grad(u_z)
        return lap

    @staticmethod
    def get_wavenumbers(Nz,dz,zType='Periodic'):
        '''
        Construct the wavenumber vectors
        Here, it is assumed that the only real values variables are used for 
            the FFTs 
        - Input - Nx,Ny,Nz - domain size
            - dx,dy,dz - physical grid spacing
        - Output - tuple k,l,m
        '''
        '''
        # Assuming the the field hold real values
        k = np.fft.rfftfreq(Nx,dx/np.pi/2)
        l = np.fft.rfftfreq(Ny,dy/np.pi/2)
        if (Nx%2 == 0) : k[-1] = 0
        if (Ny%2 == 0) : l[-1] = 0

        if (zType == 'Periodic'):
            m = np.fft.rfftfreq(Nz,dz/np.pi/2)
            if (Nz%2 == 0) : m[-1] = 0
        elif (zType == 'Cosine'):
            m = np.fft.rfftfreq(Nz*2-2,dz/np.pi/2)
            m[-1] = 0
            #m = np.fft.rfftfreq(Nz*2,0.5)+0.5/Nz
            #m = m[0:Nz]; 
            #if (Nz%2==0): m[-1] = 0
            #m /= dz/np.pi;

        k = np.array(k,ndmin=3).reshape(1,1,len(k))
        l = np.array(l,ndmin=3).reshape(1,len(l),1)
        m = np.array(m,ndmin=3).reshape(len(m),1,1)
        '''
        # Assuming the the field hold real values
        
        if (zType == 'Periodic'):
            k = np.fft.fftfreq(Nz,dz/np.pi/2)
            #if (Nz%2 == 0) : m[int(Nz/2)] = 0
        elif (zType == 'Cosine'):
            k = np.fft.fftfreq(Nz*2-2,dz/np.pi/2)
            #m[Nz-1] = 0
            #m = np.fft.rfftfreq(Nz*2,0.5)+0.5/Nz
            #m = m[0:Nz]; 
            #if (Nz%2==0): m[-1] = 0
            #m /= dz/np.pi;
        elif zType == 'Cheb':
            k = -1
        
        return k

    


    def apply_spec_filter(self,u):
        '''
        Apply the filter in spectral space 
        - Input - u - field to be filtered in spectral space
            - K_norm - normalized wavenumber array
        - Output - u_f filtered field in spectral space
        '''
        
        # Subich filter -- Subich 2013	
        f_cutoff = 0.6#0.0   # wavenumber cutoff for the filter
        f_decay  = 10.0#10.0  # decay rate of the Nyquist frequency
        f_order  = 6.0 #6   # order of the spectral filter

        
        filt = 1.0*(self.K_filt<=f_cutoff)
        filt += np.exp(-f_decay * ((self.K_filt - f_cutoff)/(1 - f_cutoff))**f_order)*(self.K_filt>f_cutoff)
        
        u_f = u*filt
        
        return u_f


    def filter_spec(self,u,zType='Periodic'):
        ''' 
        Filter a given input using a sepectral filter 
        - Input -- u - field to be filtered
            - K_norm - normalized wavenumber array (sqrt(k^2 + l^2 + m^2)(in fft order)
            - zType - String - continuation type in vertical 
                - Option = Periodic or Cosine, Default = 'Periodic
        - Output - u_f - Array of size(u) which has been filtered
        '''

        if self.zType == 'Cosine':
            u = self.flipud_even(u)

        u_f = np.real(np.fft.ifftn(self.apply_spec_filter(np.fft.fftn(u))))

        if self.zType == 'Cosine':	
            u_f = self.flipud_crop(u_f)
                
        return u_f




        




