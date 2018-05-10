import numpy as np

class MHDphase:
    def __init__(self, beta, theta, xacc=1.e-6):
        # mf: fast Mach number, ms: slow Mach number, mi: another branch of Alfven Mach number
        # normalized by normal Alfven velocity 
        self.beta = beta
        self.theta = theta
        self.mf = np.sqrt(self.rhsol(1.0+xacc, 1.e+10))
        self.ms = np.sqrt(self.rhsol(0.0, 1.0-xacc))
        self.mi = np.sqrt(self.rhsol(0.0, 1.0, alf_mode=True))

    def ax1(self, ax2):
        # determine the square of the upstream Alfven Mach number from that in the downstream
        beta = self.beta
        theta = self.theta
        tanth2 = np.tan(theta*np.pi/180.)**2
        costh2 = np.cos(theta*np.pi/180.)**2
        gamma = 5./3.

        x1 = ax2*( (gamma-1.)/gamma * ( (gamma+1.)/(gamma-1.) - tanth2 ) * \
                       (ax2-1.)**2 + tanth2*( (gamma-1.)/gamma*ax2 - 1. ) * \
                       ( ax2-2.) ) - beta/costh2*(ax2-1.)**2
        x2 = (gamma-1.)/gamma * (ax2-1.)**2/costh2 - ax2*tanth2*( (gamma-1.)/gamma*ax2 - 1. )

        return x1/x2

    def rhsol(self, x1, x2, alf_mode=False, xacc=1.e-6, iteration=100):
        beta = self.beta
        theta = self.theta
        for j in range(iteration):
            if(alf_mode):
                x1ref = 1.0
                x2ref = 1.0
            else:
                x1ref = x1
                x2ref = x2
            f2 = self.ax1(x2) - x2ref
            f1 = self.ax1(x1) - x1ref
            if f1*f2 > 0:
                print('root must be bracketed in bisection')
                break
            xmid = 0.5 * (x2 + x1)
            if(alf_mode):
                xref = 1.0
            else:
                xref = xmid
            fmid = self.ax1(xmid) - xref
            if f1*fmid > 0.0:
                x1 = xmid
            else:
                x2 = xmid
            if np.absolute(fmid) < xacc:
                return xmid
        print('lack of iteration')

    def critical_Alfven(self, step = 1.e-5):
        ymax = -1.e+10
        x = self.ms
        while x < 1.0:
            y = self.ax1(x)
            if y > ymax:
                xmax = x
                ymax = y
            x = x + step
        return xmax
