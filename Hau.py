import math
def ratio(beta, theta, ax2):
    pi = math.pi
    tanth2 = math.tan(theta*pi/180.)**2
    costh2 = math.cos(theta*pi/180.)**2
    gamma = 5./3.

    x1 = ax2*( (gamma-1.)/gamma * ( (gamma+1.)/(gamma-1.) - tanth2 ) * \
                   (ax2-1.)**2 + tanth2*( (gamma-1.)/gamma*ax2 - 1. ) * \
                   ( ax2-2.) ) - beta/costh2*(ax2-1.)**2
    x2 = (gamma-1.)/gamma * (ax2-1.)**2/costh2 - ax2*tanth2*( (gamma-1.)/gamma*ax2 - 1. )

    return x1/x2

def bysection(beta, theta, ax2, mode, xacc=1.e-6, iteration=100):
    # mode => 1: fast mode, 2: intermediate mode 3: slow mode
    if mode != 2:
        if mode == 3:
            x1 = 0.0
            x2 = 1.0
        elif mode == 1:
            x1 = 1.0 + xacc
            x2 = 1.d+10

        for j in range(iteration):
            f2 = ratio(beta, theta, x2) - x2
            f1 = ratio(beta, theta, x1) - x1
            if f1*f2 > 0:
                print('root must be bracketed in bisection')
                break
            xmid = 0.5*(x2 + x1)
            fmid = ratio(beta, theta, xmid) - xmid
            if f1*fmid > 0.0:
                x1 = xmid
            else:
                x2 = xmid
            if math.abs(fmid) < xacc:
                return xmid
    else:
        x1 = 0.0
        x2 = 1.0
        for j in range(iteration):
            f2 = ratio(beta, theta, x2) - 1.0
            f1 = ratio(beta, theta, x1) - 1.0
            if f1*f2 > 0.0:
                print('root must be bracketed in bisection')
                break
            xmid = 0.5*(x2 + x1)
            fmid = ratio(beta, theta, xmid) - 1.0
            if f1*fmid > 0.0:
                x1 = xmid
            else:
                x2 = xmid
            if math.abs(fmid) < xacc:
                return xmid
    print('lack of iteration!!')
