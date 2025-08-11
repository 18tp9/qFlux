### backtracking mock up

'''
Find the gradient at an initial point xc, take a large step in the direction of gradient
take iteratively smaller steps until function value is a specified amount below f(xc)
'''

'''
NOTE: two-way backtracking - changing alpha0 for subsequent iterations dependent on whether or not condition is reached
'''

from re import T
from termios import FF0
from Dependencies.Operators import Operator
from Dependencies.Parameters import shared_parameters

def backtrackMinimize(grad,func,eps0 = 0.05,alpha0  = 0.3,max_iter = 10):

    alpha   = alpha0 #initial/max step size
    c       = 0.05   #minimum function decrease factor
    tao       = 0.7   #contraction factor

    for n in range(max_iter):
    
        f0      = func(eps0)
        new_eps = eps0 + grad*alpha #take a step - need to confirm this math will work
        new_f   = func(new_eps)
        t       = c*grad    #not super sure if should use grad here

        ### Check if function evaluation satisfies decrease condition
        if (f0 - new_f) >= alpha*t:
            return (new_eps,new_f)
        else:
            alpha = alpha*tao   #decreases step size by contraction factor

    return ("Condition not satisfied in maximum iterations")

    






