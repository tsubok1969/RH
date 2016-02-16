import math
def Hau(beta, theta, ax2):
    pi = math.pi
    tanth2 = math.tan(theta*pi/180.)**2
    costh2 = math.cos(theta*pi/180.)**2
    gamma = 5./3.

    x1 = ax2*( (gamma-1.)/gamma * ( (gamma+1.)/(gamma-1.) - tanth2 ) * \
                   (ax2-1.)**2 + tanth2*( (gamma-1.)/gamma*ax2 - 1. ) * \
                   ( ax2-2.) ) - beta/costh2*(ax2-1.)**2
    x2 = (gamma-1.)/gamma * (ax2-1.)**2/costh2 - ax2*tanth2*( (gamma-1.)/gamma*ax2 - 1. )

    return x1/x2
