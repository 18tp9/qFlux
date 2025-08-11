import numpy as np 
import math
import matplotlib.pyplot as pl
import csv
import random

class minimization():

    def __init__(self):
        self.previous_error = 0.0
        self.current_min = 1.0
        self.MethodNames = []
        self.Errors = []
        self.Times = []

        

    @staticmethod
    def get_exponent(num):
        '''
        gets the magnitude of the exponent from scientific notation 
        '''
        numSci = np.format_float_scientific(num,exp_digits=1)
        numStr = str(numSci)
        exp = numStr.split("e")
        exp = int(exp[1])
        return exp
    
    def performance_metrics(self):
        with open('performanceMetrics.csv') as file:
            reader = csv.reader(file)
            for row in reader:
                self.MethodNames.append(row[0])
                self.Errors.append(float(row[1]))
                self.Times.append(float(row[2]))
                print (row)


    def efficiency_plots(self):
        self.performance_metrics()

        pl.show()
        pl.ion()
        x=np.arange(len(self.MethodNames))

        fig,ax = pl.subplots()
        ax2 = ax.twinx()
        ax.set_title("Method Efficiency: Error and Speed (20 Iterations)")
        ax.set_ylabel("Error")
        ax2.set_ylabel("Speed (s)")
        ax.set_xticks(x,self.MethodNames)
        errorMag = self.get_exponent(max(self.Errors))
        errorLim = 10**(-(errorMag-1))
        ax.set_ylim(0,errorLim)
        ax2.set_ylim(0,100)
        ax.set_xlim(-0.5,len(x)-0.5)
        ax2.set_xlim(-0.5,len(x)-0.5)
        ax.autoscale(False)
        ax2.autoscale(False)
        
        ax.bar((x-0.2),self.Errors,width=0.4)
        ax2.bar((x+0.2),self.Times,width=0.4,color='orange')


    @staticmethod
    def Brent(ax,cx,bx,func):
            '''
            From Press 2007
            Given a function or functor f, and given a bracketing triplet of abscissas ax, bx, cx (such
            that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine
            isolates the minimum to a fractional precision of about tol using Brent’s method. The
            abscissa of the minimum is returned as xmin, and the function value at the minimum is
            returned as min, the returned function value.
            '''
            ITMAX = 100
            CGOLD = 0.3819660
            epsilon = 2.220446049250313E-016
            ZEPS = np.sqrt(epsilon)*10         # const Doub ZEPS=numeric_limits<Doub>::epsilon()*1.0e-3;
            tol   = 3e-6


            # Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio;
            # and ZEPS is a small number that protects against trying to achieve fractional accuracy
            # for a minimum that happens to be exactly zero.

            # Doub a,b,d=0.0,etemp,fu,fv,fw,fx;
            # Doub p,q,r,tol1,tol2,u,v,w,x,xm;
            # Doub e=0.0; This will be the distance moved on
            # the step before last.

            d,e = 0,0
            # a=(ax < cx ? ax : cx)  # a and b must be in ascending order,
            # b=(ax > cx ? ax : cx)  #  but input abscissas need not be.
            a = ax if ax<cx else cx
            b = ax if ax>cx else cx

            x = w = v = bx

            fw=fv=fx=func(x);

            for iter in range(0,ITMAX):
                ## Define the midpoint 
                xm=0.5*(a+b);
                #tol2 = 2.0*(tol1=tol*abs(x)+ZEPS);
                tol1 = (tol*abs(x) + ZEPS)
                tol2 = 2.0*tol1  ### Avoid a lot of meaningless oscillations

                fmin=fx;
                xmin=x;
                ### If the difference between x and the mean is less than the tolerance
                if (abs(x-xm) <= (tol2-0.5*(b-a))):
                    # Test for done here.
                    
                    # print(f'Converged after {iter} iterations')
                    return (xmin,fmin)

                if (abs(e) > tol1) : #Construct a trial parabolic fit.
                    r=(x-w)*(fx-fv);
                    q=(x-v)*(fx-fw);
                    p=(x-v)*q-(x-w)*r;
                    q=2.0*(q-r);

                    if (q > 0.0): p = -p;
                    
                    q=abs(q);
                    etemp=e;
                    e=d;
                    
                    if (abs(p) >= abs(0.5*q*etemp)) or (p <= q*(a-x)) or (p >= q*(b-x)):
                        # The above conditions determine the acceptability of the parabolic fit. Here
                        # we take the golden section step into the larger of the two segments.
                        # d=CGOLD*(e=(x >= xm ? a-x : b-x));
                        e = a-x if x>xm else b-x
                        d = CGOLD * e


                    else:
                        d = p/q  #Take the parabolic step.
                        u = x+d;
                        ## Push away from the boundaries
                        if (u-a < tol1) or (b-u < tol1):
                            if (x<xm):
                                d = tol1
                                #u = x+push_d
                            else:
                                d = -tol1
                                #u = x-push_d

                    

                else:       ###GSS with GSS constant CGOLD
                    # d=CGOLD*(e=(x >= xm ? a-x : b-x));
                    e = a-x if x>xm else b-x
                    d = CGOLD*e
                
                ## Location of past guess? 
                #u=(abs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
                if abs(d) >= tol1:
                    u = x+d  
                elif d<0:
                    u = x - tol1 
                else:
                    u = x + tol1
                
                fu=func(u);
                # This is the one function evaluation per iteration.

                #### Housekeeping 
                if (fu <= fx) : 
                    # Now decide what to do with our function evaluation.
                    if (u >= x): a=x
                    else: b=x
                    
                    #Housekeeping follows:
                    v,w,x,= w,x,u     
                    fv,fw,fx = fw,fx,fu
                
                else :
                    if (u < x):  a=u
                    else:  b=u;
                    if (fu <= fw) or (w == x):
                        v=w;
                        w=u;
                        fv=fw;
                        fw=fu;
                    elif (fu <= fv) or (v == x) or (v == w):
                        v=u;
                        fv=fu;
                
                # Done with housekeeping. Back for
                # another iteration.


            raise ValueError("Too many iterations in brent")

    @staticmethod
    def gss(f, a, b, tol=0.1):
        """Golden-section search.

        Given a function f with a single local minimum in
        the interval [a,b], gss returns a subset interval
        [c,d] that contains the minimum with d-c <= tol.
        """
        invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
        invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2

        (a, b) = (min(a, b), max(a, b))
        h = b - a
        if h <= tol:
            return (a, b)

        # Required steps to achieve tolerance
        n = int(math.ceil(math.log(tol / h) / math.log(invphi)))
        # n=30

        c = a + invphi2 * h
        d = a + invphi * h
        yc = f(c)
        yd = f(d)

        for k in range(n-1):
            if yc < yd:  # yc > yd to find the maximum
                b = d
                d = c
                yd = yc
                h = invphi * h
                c = a + invphi2 * h
                yc = f(c)
            else:
                a = c
                c = d
                yc = yd
                h = invphi * h
                d = a + invphi * h
                yd = f(d)

        if yc < yd:
            minLoc = (a+d)/2.0
        else:
            minLoc = (c+b)/2.0
        return (minLoc, f(minLoc))

    @staticmethod
    def BackTrack(func,alpha0,rho,c1,norm2gradPhi):
        
        # print('Getting local gradient...')
        f0 = func(0)
        
        # fc1 = func(c1)

        # norm2gradPhi = (fc1-f0)/c1

        # if norm2gradPhi<0:
        #     return 1,1
        

        # if f1<(f0):
        #     alpha = 5
        #     f1 = func(alpha)
        #     return alpha,f1

        max_iter = 100; # Max number of iterations
        alpha = alpha0
        alphaMin = 0
        fMin     = f0

        # if fc1<f0:
        #     alphaMin = 1
        #     fMin = fc1

        # convCriterion = c1*f0
        convCriterion = f0 - rho*c1*norm2gradPhi
        # if convCriterion<0:

        #     return 0,1000
        #     convCriterion = f0 - rho*(c1/10)*norm2gradPhi
        #     print('Armijo-Wolfe Negative, c modified = ' + str(convCriterion))
        # else:
        print('Functional Gradient:  ' + str(norm2gradPhi))
        print('Armijo Wolfe J: ' + str(convCriterion))

        alphaArr = []
        Jarr = []

        for i in range(max_iter):

            

            # Get estimate
            fN = func(alpha)

            alphaArr.append(alpha)
            Jarr.append(fN)
            # with open('OptModule/Dependencies/Minimization/minimization.txt','a') as f:
            #         f.write(str(alpha)+','+str(f1)+'\n')
            
            # Test Amijjo / Wolfe condition
            
            if (fN <= convCriterion) and fN<f0:
                
            # if (f1 <= f0*c1):
                pl.plot(alphaArr,Jarr)
                pl.show()
                return (alpha,fN)  # Armijo/Wolfe Condition
            

            
            else:
                alpha *= rho


            # if fN<fMin:
            #     # print("yes" + str(fN))
            #     alphaMin = alpha/rho
            #     fMin     = fN

        # if fMin>=fN:
        #     # print("TRUE" + str(fN))
        #     alphaMin = 0
        #     fMin = f0

        raise ValueError("Did not converge.")
    
    @staticmethod
    def SaddleFinder(func,rho,c1):
        '''
        PABLS: Passive and Adaptive Backtracking Line Search
        '''

        f0 = func(0)
        f1 = func(1)

        max_iter = 10; # Max number of iterations
        alpha = 1
        alphaMin = 0
        fMin     = f0

        # if f1<f0:
        #     alphaMin = 1
        #     fMin = f1
        #     return alphaMin,fMin
        
        fM1 = f1
        alpha*=rho

        for i in range(max_iter):

            

            # Get error
            fN = func(alpha)

            
            # check for inflection point after 1st sub-iter
            
            if i > 0 and (fN > fM1):
                
                print("Inflection point.")
                return (alpha/rho,fM1) 
            
            #return 1 to increase gamma if second step error is greater than 1st
            elif fN>fM1:
                alphaMin = 1
                fMin = f1
                return alphaMin,fMin
            
            #apply contraction ratio
            else:
                alpha *= rho
                fM1 = fN


            if fN<fMin:
                
                alphaMin = alpha/rho
                fMin     = fN

        if fMin>=fN:
            
            alphaMin = 0
            fMin = f0

        return (alphaMin,fMin)
    
    @staticmethod
    def BigSteps(func,rho,c1 = 0.7):
        
        # print('Getting local gradient...')
        f0 = func(0)
        # fP = 1000

        f1 = func(1,printing=False)
        

        max_iter = 50; # Max number of iterations
        alpha = 1
        alphaMin = 0
        fMin     = f0

        if f1<f0:
            alphaMin = 1
            fMin = f1
        
        # c1 = 0.7
        convCriterion = c1*f0


        print('Target J: ' + str(convCriterion))

        for i in range(max_iter):

            

            # Get estimate
            fN = func(alpha)
            # with open('OptModule/Dependencies/Minimization/minimization.txt','a') as f:
            #         f.write(str(alpha)+','+str(f1)+'\n')
            
            # Test Amijjo / Wolfe condition
            
            if (fN <= convCriterion):
                
            # if (f1 <= f0*c1):

                return (alpha,fN)  # Armijo/Wolfe Condition
            

            
            else:
                alpha *= rho


        return (alphaMin,fMin)
    
    @staticmethod
    def SimpleStep(func,step):

        f0 = func(0)
        fP = 1000

        f1 = func(step)
        
        if f1>f0:
            return 0,0

        return step,f1
    
    @staticmethod
    def Plotting(func,gamma,iter,init=1,rho=0.8):
        

        

        max_iter = 50; # Max number of iterations

        print('Plotting Run')
        
        f0 = func(0)

        alpha = init
        alphaArr = []
        Jarr = []

        for i in range(max_iter):
            
            if alpha <0:
                break
            
            print("PR")
            # Get estimate
            fN = func(alpha)

            alphaArr.append(alpha*gamma)
            Jarr.append(fN)
            # with open('OptModule/Dependencies/Minimization/minimization.txt','a') as f:
            #         f.write(str(alpha)+','+str(f1)+'\n')
            
            # Test Amijjo / Wolfe condition
            
            alpha  = alpha - 0.025*init


            # if fN<fMin:
            #     # print("yes" + str(fN))
            #     alphaMin = alpha/rho
            #     fMin     = fN

        # if fMin>=fN:
        #     # print("TRUE" + str(fN))
        #     alphaMin = 0
        #     fMin = f0
        alphaArr.append(0)
        Jarr.append(f0)

        # pl.plot(alphaArr,Jarr,'-o',markersize=5,markerfacecolor=None)
        # pl.ylabel(r'$J = \int\int(T-T_m)^2 dzdt$')
        # pl.xlabel(r'$\epsilon$')
        # pl.title('Iteration ' + str(iter))
        # pl.show()

        return alphaArr,Jarr
        
        

