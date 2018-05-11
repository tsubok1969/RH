import sys
import numpy as np

class MHDphase:
    def __init__(self, beta, theta, xacc=1.e-6):
        # beta and theta : upstream quantities
        # squared mf: fast MS, ms: slow MS, mi: another branch of Alfven velocity
        # normalized by the upstream normal Alfven velocity 
        self.beta = beta
        self.theta = theta
        self.mf = self.rhsol(1.0+xacc, 1.e+10)
        self.ms = self.rhsol(0.0, 1.0-xacc)
        self.mi = self.rhsol(0.0, 1.0, alf_mode=True)
        self.mcd, self.mcu = self.critical_Alfven()

    def ax1(self, ax2, gamma=5./3.):
        # determine the square of the upstream Alfven Mach number from that in the downstream
        beta = self.beta
        theta = self.theta
        deg2rad = np.pi/180.
        tanth2 = np.tan(theta*deg2rad)**2
        costh2 = np.cos(theta*deg2rad)**2

        x1 = ax2*( (gamma-1.)/gamma * ( (gamma+1.)/(gamma-1.) - tanth2 ) * \
                       (ax2-1.)**2 + tanth2*( (gamma-1.)/gamma*ax2 - 1. ) * \
                       ( ax2-2.) ) - beta/costh2*(ax2-1.)**2
        x2 = (gamma-1.)/gamma * (ax2-1.)**2/costh2 - ax2*tanth2*( (gamma-1.)/gamma*ax2 - 1. )

        return x1/x2

    def rhsol(self, x1, x2, alf_mode=False, xacc=1.e-6, iteration=1000):
        beta = self.beta
        theta = self.theta
        for j in range(iteration):
            if(alf_mode):
                x1ref = x2ref = 1.0
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
        # determine the critical Alfven Mach number
        # Ma_down > xmax: weak intermediate shock
        # Ma_down < xmax: strong intermediate shock
        # ymax: critical Alfven Mach number of upstream flow
        ymax = -1.e+10
        x = (self.ms)**2
        while x < 1.0:
            y = self.ax1(x)
            if y > ymax:
                xmax = x
                ymax = y
            x = x + step
        return xmax, ymax

    def solsearch(self, xdown, mode):
        # mode 1: fast shock 2: slow shock 3: IS
        ax2 = xdown**2
        deg2rad = np.pi/180.
        while(True):
            if mode == 1:
                if ax2 > self.mf:
                    break
            elif mode == 2:
                if self.ms < ax2 < self.mi:
                    break
            else:
                if self.mi < ax2 < 1.0:
                    break
            print('value is not suitable for RH')
            sys.exit()

        ax1 = self.ax1(ax2)
        r = ax1/ax2
        thetad = np.rad2deg( np.arctan( (ax1-1.)/(ax2-1.)*np.tan(self.theta*deg2rad) ) )
        bd = np.cos(self.theta*deg2rad)/np.cos(thetad*deg2rad)
        bt = bd * np.sin(thetad*deg2rad)
        uin = (1.-1./r)*np.sqrt(ax1)*np.cos(self.theta*deg2rad)
        vshock = np.sqrt(ax2/r)*np.cos(self.theta*deg2rad)
        print('compressional ratio    : ' + str(r) )
        print('downstream field angle : ' + str(thetad) )
        print('downstream field mag.  : ' + str(bd) )
        print('injection velocity     : ' + str(uin) )
        print('shock velocity         : ' + str(vshock) )
