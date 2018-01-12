from tdd.twindyn import Dislocation, Simulation
import numpy as np
import matplotlib.pyplot as p

class analyticalThickness(object):
    def __init__(self, isomatprops, b, h):
        self.isomatprops = isomatprops
        self.b = b
        self.h = h

    def thinPts(self, x, tau, twoSided=True):
        '''Return twin points as a function of x. If twoSided, will return pos and negative values for plotting 
           x = [-1,0,1] may return x[0,1,0,0,-1,0]
        '''
        width = float(x[-1])-float(x[0])
        const = self.h*tau*(1.-self.isomatprops.v)/(self.isomatprops.mu*self.b)        
        t = np.where(abs(x)<width/2., const*np.sqrt((width/2.)**2 - x**2),0.0)

        if twoSided:
            x = np.array([x,x[::-1]]).flatten()
            t = np.array([t,-t[::-1]]).flatten()
            
        return x,t
             
    def thickConstEZ(self, width, tau):
        return tau*width*self.h/(self.isomatprops.mu*self.b)
        
        
    def thickConstHard(self, width, tau, nmax=10000):
        tb = [Dislocation(width/2., 0.0, 0.0, self.b, False), Dislocation(-width/2., 0.0, np.pi, self.b, False)]
        tbsim = Simulation(tb, self.isomatprops, tau=tau, vertsym=True)

        jprev=0
        for j in xrange(jprev, nmax):
            ycur = self.h*(j+1)
            tbsim.adddislo(Dislocation(width/2, ycur, 0.0, self.b))
            tbsim.adddislo(Dislocation(-width/2, ycur, np.pi, self.b))
            sigxx, sigzz, sigxz = tbsim.calc_stresses_point(0.0, ycur)
            if sigxz < 0.0:
                break
        return 2.0*ycur

