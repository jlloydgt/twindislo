#dislocation angles go ccw
import numpy as np
from math import pi, sin, cos
from staticsolns import sigstat
import itertools
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as p
try:
    import cPickle as pickle
except:
    import pickle
    
    
class IsoMaterialProps(object):
    def __init__(self, lam, mu, rho):
        self.lam = lam
        self.mu = mu
        self.v = lam/(2.*(lam+mu))
        self.e = 2.*mu*(1.+self.v)
        self.k = lam+2.*mu/3.
        self.rho = rho
        
    @staticmethod
    def mg():
        return IsoMaterialProps(24.1333e3, 16.9e3, 1.74e-3)

class Dislocation:
    newid = itertools.count().next
    def __init__(self, x, z, ang, burg, state='m'):
        self.id = Dislocation.newid()
        self.x = x
        self.z = z
        self.ang = ang
        self.burg = burg
        self.state = state #'m', 'gb', 's' 

        #utility properties to store during time step
        self.dx = None
        self.dz = None

    def __str__(self):
        if self.state=='m':
            return ('m: '  + str(self.x) + ', ' + str(self.z))
        elif self.state=='s':
            return ('s: '  + str(self.x) + ', ' + str(self.z))
        elif self.state=='gb':
            return ('gb: ' + str(self.x) + ', ' + str(self.z))
        else:
            return ('??: '  + str(self.x) + ', ' + str(self.z))

    def calc_force_disp(self, dislolist, matprops, dt, tau):
        dislolistc = dislolist[:]
        try:
            dislolistc.remove(self)
        except:
            pass
        
        #print 'dislo used to calc force on ' + str(self.x) + ', ' + str(self.z)
        #for dislo in dislolistc:
        #    print dislo.x, dislo.z
               
        #calculation force on a dislocation at this point in the local dislocation frame
        sigxx, sigzz, sigxz = 0.0, 0.0, 0.0
        for dislo in dislolistc:
            sigxxadd, sigzzadd, sigxzadd = sigstat(self.x - dislo.x, self.z - dislo.z, dislo.ang, dislo.burg, matprops)
            sigxx, sigzz, sigxz = sigxx + sigxxadd, sigzz + sigzzadd, sigxz + sigxzadd 
        sx, sz, mx, mz = cos(self.ang), sin(self.ang), -sin(self.ang), cos(self.ang)
        force = self.burg*(sigxx*(sx*mx)+sigzz*(sz*mz)+sigxz*(sx*mz+sz*mx)+tau*(sx*mz)) #, -self.burg*sigxx
        
        if self.state=='m':
            return force, force*dt  
        else:
            return force, 0.0 

    def combineAngle(self, angle):
        '''Return burgers vector and angle assuming angle is the same dislocation type 
        '''

        b1x, b1y = self.burg*np.cos(self.ang), self.burg*np.sin(self.ang) 
        b2x, b2y = self.burg*np.cos(angle), self.burg*np.sin(angle) 

        bnx, bny = b1x+b2x, b1y+b2y
        if abs(bnx)/self.burg < 1e-5:
            if bny > 0:
                bnang = pi/2.0
            else:
                bnang = -pi/2.0
        else:
            bnang = np.arctan(bny/bnx)
        bnmag = np.sqrt(bnx**2+bny**2)

        return bnmag, bnang    

    def combineAngleReturnDislo(self, angle):
        '''Return burgers vector and angle assuming angle is the same dislocation type 
        '''

        b1x, b1y = self.burg*np.cos(self.ang), self.burg*np.sin(self.ang) 
        b2x, b2y = self.burg*np.cos(angle), self.burg*np.sin(angle) 
        
        bnx, bny = b1x+b2x, b1y+b2y
        bnang = np.arctan(bny/bnx)
        bnmag = np.sqrt(bnx**2+bny**2)

        return Dislocation(self.x, self.z, bnang, bnmag, state=self.state)

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
        tb = [Dislocation(width/2., 0.0, 0.0, self.b, False)]
        tbsim = Simulation(tb, self.isomatprops, tau=tau, vertsym=True, horsym=True)

        jprev=0
        for j in xrange(jprev, nmax):
            ycur = self.h*(j+1)
            tbsim.adddislo(Dislocation(width/2, ycur, 0.0, self.b))
            sigxx, sigzz, sigxz = tbsim.calc_stresses_point(0.0, ycur)
            if sigxz < 0.0:
                break
        return ycur, tbsim

    def thickConstHardAngle(self, width, tau, rangle, langle, nmax=10000, savefile=None):
        '''rangle is ccw in pos x dir, langle is ccw in neg x dir
        '''
        rinit = Dislocation(width/2., 0.0, 0.0, self.b, 'gb')
        br, angr = rinit.combineAngle(rangle+pi)
        dxr, dyr = width*cos(rangle), width*sin(rangle)

        linit = Dislocation(-width/2., 0.0, pi, self.b, 'gb')
        bl, angl = linit.combineAngle(langle)        
        dxl, dyl = width*cos(pi+langle), width*sin(pi+langle)        
        
        tb = [Dislocation(width/2., 0.0, angr, br, 'gb'), Dislocation(-width/2., 0.0, pi+angl, bl, 'gb'), 
              Dislocation(width/2.+dxr, dyr, rangle, self.b, 'gb'), Dislocation(-width/2.+dxl, dyl, pi+langle, self.b, 'gb')]
        tbsim = Simulation(tb, self.isomatprops, tau=tau)

        for j in xrange(nmax):
            ycur = self.h*(j+1)
            tbsim.adddislo(Dislocation(width/2, ycur, angr, br))
            tbsim.adddislo(Dislocation(width/2, -ycur, angr, br))
            tbsim.adddislo(Dislocation(-width/2, ycur, pi+angl, bl))
            tbsim.adddislo(Dislocation(-width/2, -ycur, pi+angl, bl))

            tbsim.adddislo(Dislocation(width/2+dxr, ycur+dyr, rangle, self.b))
            tbsim.adddislo(Dislocation(width/2+dxr, -ycur+dyr, rangle, self.b))
            tbsim.adddislo(Dislocation(-width/2+dxl, ycur+dyl, pi+langle, self.b))
            tbsim.adddislo(Dislocation(-width/2+dxl, -ycur+dyl, pi+langle, self.b))
            
            sigxx, sigzz, sigxz = tbsim.calc_stresses_point(0.0, ycur)
            if sigxz < 0.0:
                if savefile:
                    #f = p.figure(figsize=(6,6))
                    #pp = tbsim.plotalldislocations([-2.0*width, 2.0*width, -2.0*ycur, 2.0*ycur])
                    #f.savefig(savefile, dpi=300) 
                    tbsim.plotsig([-5.0*width, 5.0*width, -width, width], pmax=50.0e3, disloscale=1.0, fname=savefile)
                break
        return ycur

    def thickConstHardAngle2(self, width, tau, angle, nmax=10000, savefile=None, ngrains=1):
        '''angle is ccw in pos x dir
        '''
        
        #right side
        rinit = Dislocation(width/2., 0.0, 0.0, self.b, 'gb')
        rinterior = Dislocation(width/2., 0.0, angle, self.b, 'gb')        
        br, angr = rinit.combineAngle(angle+pi)
        brint, angrint = rinterior.combineAngle(pi-angle)
        dxr, dyr = width*cos(angle), width*sin(angle)
        
        #left side
        linit = Dislocation(-width/2., 0.0, pi, self.b, 'gb')
        linterior = Dislocation(width/2., 0.0, pi+angle, self.b, 'gb')                       
        bl, angl = linit.combineAngle(angle)
        blint, anglint = linterior.combineAngle(-angle)        
        dxl, dyl = width*cos(pi+angle), width*sin(pi+angle)        

        tb = []
        tbsim = Simulation(tb, self.isomatprops, tau=tau)
        
        for j in xrange(nmax):
            ycur = self.h*j
            
            #grain: 0
            tbsim.adddislo(Dislocation(width/2, ycur, angr, br))
            tbsim.adddislo(Dislocation(-width/2, ycur, pi+angl, bl))
            if j!=0:
                tbsim.adddislo(Dislocation(width/2, -ycur, angr, br))
                tbsim.adddislo(Dislocation(-width/2, -ycur, pi+angl, bl))            
            
            #grains: 1,ngrains-1
            for k in xrange(1,ngrains):
                angsign = -1.0+2.0*(k%2) #1 odd -1 even
                vertdisp = float(k%2) #0 even 1 odd
    
                #add inner dislocations
                tbsim.adddislo(Dislocation(width/2+k*dxr, ycur+dyr*vertdisp, angsign*angrint, brint))
                tbsim.adddislo(Dislocation(-width/2+k*dxl, ycur+dyl*vertdisp, pi+angsign*anglint, blint))
                if j!=0:
                    tbsim.adddislo(Dislocation(width/2+k*dxr, -ycur+dyr*vertdisp, angsign*angrint, brint))
                    tbsim.adddislo(Dislocation(-width/2+k*dxl, -ycur+dyl*vertdisp, pi+angsign*anglint, blint))

            #grain: ngrains
            k = ngrains
            angsign = -1.0+2.0*(k%2)  #1 odd -1 even
            vertdisp = float(k%2) #0 even 1 odd              
            
            tbsim.adddislo(Dislocation(width/2+k*dxr, ycur+dyr*vertdisp, angle*angsign, self.b))
            tbsim.adddislo(Dislocation(-width/2+k*dxl, ycur+dyl*vertdisp, pi+angle*angsign, self.b))
            if j!=0:
                tbsim.adddislo(Dislocation(width/2+k*dxr, -ycur+dyr*vertdisp, angle*angsign, self.b))
                tbsim.adddislo(Dislocation(-width/2+k*dxl, -ycur+dyl*vertdisp, pi+angle*angsign, self.b))

            sigxx, sigzz, sigxz = tbsim.calc_stresses_point(0.0, ycur)
            if sigxz < 0.0:
                if savefile:
                    tbsim.plotsig([-5.0*width, 5.0*width, -width, width], pmax=50.0e3, disloscale=1.0, fname=savefile)
                break
        return ycur, tbsim

class Simulation:
    #dislolist are the dislocations that are used in computations, more may exist in getalldislo
    def __init__(self, dislolist, matprops, t=0.0, nsteps=0, tau=0.0, vertsym=None, horsym=None, rwall=None, lwall=None, rang=None, lang=None,tauxfer=None, dxxfer=None):
        self.dislolist = dislolist
        self.matprops = matprops
        self.t = t
        self.nsteps = nsteps 
        self.tau=tau
        self.lwall = lwall
        self.rwall = rwall
        self.lang = lang
        self.rang = rang
        self.tauxfer=tauxfer
        self.dxxfer=dxxfer
        if vertsym is None:
            self.vertsym = None
        else:            
            self.vertsym = vertsym
        if horsym is None:
            self.horsym = None
        else:            
            self.horsym = horsym
        
        
    def __str__(self):
        ans = ['--------------------\n']
        ans.append('Step: ' + str(self.nsteps) + '\n')
        ans.append('Time: ' + str(self.t) + '\n')
        ans.append('Shear stress: ' + str(self.tau) + ' MPa\n')
        for dislo in self.getalldislo():
            ans.append(dislo.__str__()+'\n') #appends str representation of dislocations
        ans.append('--------------------\n')                    
        return ''.join(ans)        

    def save(self, filename):
        d = {}
        #sim attrs
        d['t'], d['nsteps'], d['tau'], d['vertsym'], d['horsym'], d['rwall'], d['lwall'], d['lang'], d['rang'], d['tauxfer'], d['dxxfer'] = self.t, self.nsteps, self.tau, self.vertsym, self.horsym, self.rwall, self.lwall, self.lang, self.rang, self.tauxfer, self.dxxfer
        d['lam'], d['mu'], d['rho'] = self.matprops.lam, self.matprops.mu, self.matprops.rho     
        d['x'], d['z'], d['ang'], d['burg'], d['state']  = [], [], [], [], []
        for dislo in self.dislolist:
            d['x'].append(dislo.x)
            d['z'].append(dislo.z)
            d['ang'].append(dislo.ang)
            d['burg'].append(dislo.burg)
            d['state'].append(dislo.state)
        f = open(filename, 'w')
        pickle.dump(d, f)
        f.close()               
    
    @staticmethod
    def load(filename):
        f = open(filename, 'r')
        d = pickle.load(f)
        f.close()
        
        matprops = IsoMaterialProps(d['lam'], d['mu'], d['rho'])
        dislolist = []
        for i in xrange(len(d['x'])):
            dislolist.append(Dislocation(d['x'][i], d['z'][i], d['ang'][i], d['burg'][i], state=d['state'][i]))
        return Simulation(dislolist, matprops, tau=d['tau'], vertsym=d['vertsym'], horsym=d['horsym'], t=d['t'], nsteps=d['nsteps'], rwall=d['rwall'], lwall=d['lwall'], rang=d['rang'], lang=d['lang'], tauxfer=d['tauxfer'], dxxfer=d['dxxfer'])
        
    #helper function for getting both regular and mirrored dislocations 
    def getalldislo(self, excldislo=None, poszonly=False, negzonly=False):
        dislolistc = self.dislolist[:]

        #add vertical symmetry 
        if self.vertsym == True:
            tempdislolist = dislolistc[:]
            for dislo in tempdislolist:
                if dislo.z > 0:
                    newdislo = Dislocation(dislo.x, -dislo.z, dislo.ang, dislo.burg, state=dislo.state)
                    dislolistc.append(newdislo)
        #add horizontal symmetry to vertical symmetry
        if self.horsym == True: 
            tempdislolist = dislolistc[:]
            for dislo in tempdislolist:
                if dislo.x > 0:
                    newdislo = Dislocation(-dislo.x, dislo.z, dislo.ang+pi, dislo.burg, state=dislo.state)
                    dislolistc.append(newdislo)

        #remove dislocation (and all it's mirrored dislocations) from the list 
        if excldislo:
            dislolistc.remove(excldislo)

        #remove dislocations based on pos or neg z
        if poszonly:
            tempdislolist = dislolistc[:]
            for dislo in tempdislolist:
                if dislo.z < 0:
                    dislolistc.remove(dislo)
                    
        if negzonly:
            tempdislolist = dislolistc[:]
            for dislo in tempdislolist:
                if dislo.z > 0:
                    dislolistc.remove(dislo)                    
                    
        return dislolistc
        
    def getxmaxnegxminpos(self, big=1.e8):
        xmaxneg, xminpos = -big, big
        for dislo in self.dislolist:
            if ((dislo.x < 0.0) and (dislo.x > xmaxneg)): xmaxneg = dislo.x
            if ((dislo.x > 0.0) and (dislo.x < xminpos)): xminpos = dislo.x         
        return xmaxneg, xminpos

    def getxminpos(self, big=1.e8):
        xminpos = big
        for dislo in self.dislolist:
            if ((dislo.x > 0.0) and (dislo.x < xminpos)): xminpos = dislo.x         
        return xminpos

    def getzmax(self):
        zmax = -1e-8
        for dislo in self.getalldislo():
            if dislo.z > zmax:
                zmax = dislo.z
        return zmax

    def getzmin(self):
        zmin = 1e8
        for dislo in self.getalldislo():
            if dislo.z < zmin:
                zmin = dislo.z
        return zmin

    def finddislo(self, nid):
        for dislo in self.getalldislo():
            if dislo.id == nid:
                return dislo

    def xinsertpos(self, dx0=0.01, dxfactor=1.2, big=1.e8, small=1e-8, dxmin=False):
        ''' Calculate where to insert next positive and negative dislocation        
            For example in pos case - if x = 2,3,5,7, dx=1.1, xmin = 2, xmin2=3, so xinsert = 2 - 1.1*(3-2) = 0.9
            
            If xpos becomes negative, or xneg becomes positive, converged will return True
        '''
                
        converged = False                
        xminpos = big
        xmin2pos = big*1.001        
        npos = 0

        #calculate min, and second min pos values, max and second max neg values
        for dislo in self.dislolist:                                
            if dislo.x > 0.0:
                npos = npos + 1
                if dislo.x < xminpos: 
                    xmin2pos = xminpos
                    xminpos = dislo.x    
                elif dislo.x < xmin2pos:
                    xmin2pos = dislo.x
        
        #calculate next xminpos, xmaxneg
        if npos == 0:
            raise AttributeError('xinsertnegpos requires at least 1 dislocation in pos and neg x, exiting')

        if npos == 1:
            xinsertpos = xminpos - dx0
        else:
            dxpos = xmin2pos - xminpos #a positive number
            
            #handle small values
            if dxmin:
                if abs(dxpos) < dxmin: dxpos = dxmin
            else:
                if abs(dxpos) < small: dxpos = dx0
                
            xinsertpos = xminpos - dxfactor*dxpos

        if xinsertpos < 0.0:
            converged = True
                  
        return xinsertpos, converged

    def xinsertnegpos(self, dx0=0.01, dxfactor=1.2, big=1.e8, small=1e-8, dxmin=False, dxmax=False, poszonly=False, negzonly=False, xmid=0.0):
        ''' Calculate where to insert next positive and negative dislocation        
            For example in pos case - if x = 2,3,5,7, dx=1.1, xmin = 2, xmin2=3, so xinsert = 2 - 1.1*(3-2) = 0.9
            
            If xpos becomes negative, or xneg becomes positive, converged will return True
            
            xmid is used to determine what to consider a "positive" and "negative" dislocationS
        '''
                
        converged = False                
        xmaxneg, xminpos = -big, big
        xmax2neg, xmin2pos = -big*1.001, big*1.001        
        npos, nneg = 0,0

        #calculate min, and second min pos values, max and second max neg values
        for dislo in self.getalldislo(poszonly=poszonly, negzonly=negzonly):
            if dislo.x < xmid:
                nneg = nneg + 1
                if dislo.x > xmaxneg: 
                    xmax2neg = xmaxneg
                    xmaxneg = dislo.x
                elif dislo.x > xmax2neg:
                    xmax2neg = dislo.x                                        
            if dislo.x > xmid:
                npos = npos + 1
                if dislo.x < xminpos: 
                    xmin2pos = xminpos
                    xminpos = dislo.x    
                elif dislo.x < xmin2pos:
                    xmin2pos = dislo.x

        #calculate next xminpos, xmaxneg
        if nneg == 0 or npos == 0:
            raise AttributeError('xinsertnegpos requires at least 1 dislocation in pos and neg x, exiting')
        
        if nneg == 1:
            xinsertneg = xmaxneg + dx0
        else:
            dxneg = xmaxneg - xmax2neg #a positive number

            #handle small values
            if dxmin:
                if abs(dxneg) < dxmin: dxneg = dxmin
            else:       
                if abs(dxneg) < small: dxneg = dx0
            #handle large values
            if dxmax:
                if abs(dxneg) > dxmax: dxneg = dxmax
            
            xinsertneg = xmaxneg + dxfactor*dxneg
            
        if npos == 1:
            xinsertpos = xminpos - dx0
        else:
            dxpos = xmin2pos - xminpos #a positive number
            
            #handle small values
            if dxmin:
                if abs(dxpos) < dxmin: dxpos = dxmin
            else:
                if abs(dxpos) < small: dxpos = dx0
            #handle large values            
            if dxmax:
                if abs(dxpos) > dxmin: dxpos=dxmax
            
            xinsertpos = xminpos - dxfactor*dxpos

        if xinsertneg > xinsertpos: converged = True
    
        return xinsertneg, xinsertpos, converged

    '''
    def xinsertnegpos2(self, dx0=0.01, dxfactor=1.2, big=1.e8, small=1e-8, dxmin=False, dxmax=False, poszonly=False, negzonly=False):
                
        converged = False                
        zmin, zmax = 1e8, -1e8

        #calculate min, and second min pos values, max and second max neg values
        for dislo in self.getalldislo(poszonly=poszonly, negzonly=negzonly):
            if dislo.x <= self.rwall and dislo.x >=self.lwall: #it is in the grain
                if dislo.z >= zmax: zmax = dislo.z
                if dislo.z < zmin: zmin = dislo.z
        
        xlist = []      
        
        if poszonly:
            for dislo in self.getalldislo(poszonly=poszonly, negzonly=negzonly):
                if dislo.x <= self.rwall and dislo.x >=self.lwall: #it is in the grain
                    if dislo.z == zmax: 
                        xlist.append(dislo.x)   
        if negzonly:
            for dislo in self.getalldislo(poszonly=poszonly, negzonly=negzonly):
                if dislo.x <= self.rwall and dislo.x >=self.lwall: #it is in the grain
                    if dislo.z == zmin: 
                        xlist.append(dislo.x)                         

        if len(xlist != 2): raise AttributeError('for some reason xlist is not 2, exiting')
        xlist.sort()
        mid = (x[0]+x[1])/2.0
    '''        
        

    def step(self, dxmax, big=1e10):

        #calculate max f, then associated dt (only for mobile dislocations). calc forces due to sym dislocations too
        fmax = 0.0

        for dislo in self.dislolist:
            if dislo.state=='m':
                force, disp = dislo.calc_force_disp(self.getalldislo(), self.matprops, 1.0, self.tau)
                if abs(force)>fmax: 
                    fmax = abs(force)
                    maxdislo = dislo    

        if abs(fmax) > 1e-10:
            dt = dxmax / fmax
        else:
            print 'simulation probably converged, fmax is ', fmax
            print 'setting dt to 1.0'
            dt = 1.0        


        #keep xferring dislocations in grain boundary
        if (self.rang is not None) or (self.lang is not None):        
            ltrans, rtrans = True, True
        else:
            ltrans, rtrans = False, False
        
        while ltrans or rtrans: #transfer until 
            ltrans, rtrans = False, False
            lf, rf = -big, -big

            for dislo in self.dislolist:
                if dislo.state=='gb':
                    
                    #xfer dislocation with highest tau right and left of wall
                    if abs(dislo.x-self.rwall)<1e-6:              
                        tmp = Dislocation(dislo.x+self.dxxfer*cos(self.rang), dislo.z+self.dxxfer*sin(self.rang), self.rang, dislo.burg)
                        force, disp = tmp.calc_force_disp(self.getalldislo(), self.matprops, dt, self.tau)
                        tau = force/dislo.burg
                        if tau > self.tauxfer:
                            rtrans = True
                            if tau > rf:
                                rf, rtid = tau, dislo.id
                                
                    elif abs(dislo.x-self.lwall)<1e-6: 
                        tmp = Dislocation(dislo.x+self.dxxfer*cos(self.lang+pi), dislo.z+self.dxxfer*sin(self.lang+pi), self.lang+pi, dislo.burg)
                        force, disp = tmp.calc_force_disp(self.getalldislo(), self.matprops, dt, self.tau)
                        tau = force/dislo.burg
                        if tau > self.tauxfer:
                            ltrans = True
                            if tau > lf:
                                lf, ltid = tau, dislo.id
                    else:
                        print 'Warning: gb dislocation at ', dislo.x, ' not in either wall'

            #transmute dislocation with highest stress
            if rtrans:
                for dislo in self.dislolist:
                    if dislo.id == rtid:
                        #make new sessile dislocation
                        dxr, dzr = 2.0*self.rwall*cos(self.rang), 2.0*self.rwall*sin(self.rang)
                        self.adddislo(Dislocation(dislo.x+dxr, dislo.z+dzr, self.rang, dislo.burg, 's'))
                        
                        #trasnmute existing dislocation
                        b, ang = dislo.combineAngle(pi+self.rang) #remember you are combining with the negative dislocation
                        dislo.burg, dislo.ang, dislo.state = b, ang, 's'
                        #print 'transmuted on right: ', dislo

            if ltrans:
                for dislo in self.dislolist:
                    if dislo.id == ltid:
                        #make a new sessile dislocation
                        dxl, dzl = 2.0*abs(self.lwall)*cos(pi+self.lang), 2.0*abs(self.lwall)*sin(pi+self.lang)
                        self.adddislo(Dislocation(dislo.x+dxl, dislo.z+dzl, pi+self.lang, dislo.burg, 's'))                        
                        
                        #transmute existing dislocation
                        b, ang = dislo.combineAngle(self.lang) #remember you are combining with the negative dislocation
                        dislo.burg, dislo.ang, dislo.state = b, ang+pi, 's'
                        #print 'transmuted on left: ', dislo

        
        #keep putting dislocations in grain boundary based on furthest past boundary
        if (self.rwall is not None) or (self.lwall is not None):        
            rimob, limob = True, True
        else:
            rimob, limob = False, False
        while rimob or limob: #immobilize dislocations until 
            rimob, limob =  False, False
            lx, rx = big, -big
            for dislo in self.dislolist:         
                if dislo.state=='m':
                    force, disp = dislo.calc_force_disp(self.getalldislo(dislo), self.matprops, dt, self.tau)

                    #find the id of the dislocation that goes furthest past the wall                    
                    newx = dislo.x + disp*cos(dislo.ang)
                    if self.rwall is not None:
                        if newx > self.rwall:
                            rimob = True
                            if newx > rx:
                                rx, rid = newx, dislo.id

                    if self.lwall is not None:
                        if newx < self.lwall:
                            limob = True
                            if newx < lx:
                                lx, lid = newx, dislo.id

            #immobilize dislocations
            if rimob:
                for dislo in self.dislolist:
                    if dislo.id == rid:
                        #print 'immobilized on right: ', dislo
                        dislo.x = self.rwall
                        dislo.state = 'gb'
            if limob:
                for dislo in self.dislolist:
                    if dislo.id == lid:
                        #print 'immobilized on left: ', dislo
                        dislo.x = self.lwall
                        dislo.state = 'gb'                        

        #calculate new positions (use original positions and not updated ones)
        fmax=0.0
        for dislo in self.dislolist:                                    
            if dislo.state == 'm':
                force, disp = dislo.calc_force_disp(self.getalldislo(dislo), self.matprops, dt, self.tau)
                dislo.dx = disp*cos(dislo.ang)
                dislo.dz = disp*sin(dislo.ang)
                if abs(force) > fmax: 
                    fmax = abs(force)
                    maxforce = force          
                    maxdislo = dislo

        #update all new positions
        for dislo in self.dislolist:
            if dislo.state == 'm':
                dislo.x = dislo.x + dislo.dx
                dislo.z = dislo.z + dislo.dz
                        
        #record simulation characteristics                          
        self.maxfdislo = maxdislo
        self.t = self.t + dt
        self.nsteps = self.nsteps + 1        
        self.fmax = fmax
        #print fmax, maxdislo.id

    def adddislo(self, dislo): 
        self.dislolist.append(dislo)

    def steptot(self, ttot, dxmax):
        while (self.t < ttot):
            self.step(dxmax)
    
    def stepton(self, nsteps, dxmax):
        while (self.nsteps < nsteps):
            self.step(dxmax)
            #print self.fmax, self.t    

    def stepton_fcriterion(self, nsteps, dxmax, fcutoff):
        while (self.nsteps < nsteps):
            self.step(dxmax)
            print 'step ', str(self.nsteps), self.maxfdislo, self.fmax
            #print '-------------'             
            if self.fmax < fcutoff:    
                break

    def stepton_fcriterion2(self, nsteps, dxmax, fcutoff):
        fmaxprev = 1e8
        while (self.nsteps < nsteps):
            self.step(dxmax)
            if self.fmax < fcutoff:    
                break
            divby = 1.0
            while self.fmax > fmaxprev:
                
                print 'subcycling'
                #reset simulation characteristics
                for dislo in self.dislolist:
                    if dislo.state == 'm':
                        dislo.x = dislo.x - dislo.dx
                        dislo.z = dislo.z - dislo.dz
                self.nsteps = self.nsteps - 1       
                                        
                #record simulation characteristics                          
                divby = divby*2.0
                self.step(dxmax/divby)
            fmaxprev = self.fmax


    def stepton_fcriterion3(self, nsteps, dxmax, fcutoff, ncutstart=4, ncutint=100, dxcut=2.0, fcut=2.0e-5):
        maxflist = []
        while (self.nsteps < nsteps):
            self.step(dxmax)
            maxflist.append(self.fmax)
            if self.fmax < fcutoff:    
                break
            else:
                if self.nsteps >= ncutstart: 
                    f1, f2 = abs(maxflist[-1]/maxflist[-3]-1.0), abs(maxflist[-2]/maxflist[-4]-1.0)
                    #print 'step ', str(self.nsteps), self.maxfdislo, self.fmax, f1, f2                            
                    if self.nsteps % ncutint == 0:
                        if f1 < fcut and f2 < fcut: 
                            dxmax = dxmax/dxcut
                            print 'cutting step'
                else:
                    pass #print 'step ', str(self.nsteps), self.maxfdislo, self.fmax



    def stepton_fcriterion_relf(self, nsteps, dxmax, fcutoff, neval=4, dxcut=2.0, fcut=1e-7, finc=1e-7):
        ''' same as stepton_fcriterion except that it includes logic so that if a single dislocation flutters 
            back and forth neval times in a row, then it divides dxmax by dxcut
        '''
        maxflist = []
        maxflistid = []
        maxflistodd = []
        dxorig = dxmax
        
        while (self.nsteps < nsteps):
            self.step(dxmax)            
            if self.fmax < fcutoff:    
                break

            maxflist.append(self.fmax)
            maxflistid.append(self.maxfdislo.id)
            if self.nsteps % 2 == 1: maxflistodd.append(self.fmax)

            if len(maxflist) == neval:
                multstep = True
                if all(maxflistid[0] == item for item in maxflistid) and dxmax<dxorig: #double step routine - only up to old dxmax
                    for i in xrange(1,neval-1):
                        if multstep == True:
                            multstep = abs(maxflist[i]/maxflist[-1] - 1.0) < finc  
                    if multstep == True:
                        #print 'inc step'
                        dxmax = dxmax*dxcut                  
                else: #cut step routine
                    for i in xrange(1,neval-1):
                        if multstep == True:
                            multstep = abs(maxflist[i]/maxflist[-1] - 1.0) < fcut               
                    if multstep == True:
                        #print 'cutting step'
                        dxmax = dxmax/dxcut 

                maxflist = []
                maxflistid = []
            
            if len(maxflistodd) == neval: #only perform evaluation every neval
                cutstep = True
                for i in xrange(1,neval-1):
                    if cutstep == True:
                        cutstep = abs(maxflistodd[i]/maxflistodd[-1] - 1.0) < fcut
                if cutstep == True:                                        
                    #print 'cutting step odd'
                    dxmax = dxmax/dxcut
                maxflistodd = []
                
    def stepton_fcriterion_autoconverge(self, nsteps, dxmax, fcutoff, neval=4, dxcut=2.0):
        ''' same as stepton_fcriterion except that it includes logic so that if a single dislocation flutters 
            back and forth neval times in a row, then it divides dxmax by dxcut
        '''
        maxdislolistx = []
        maxdislolistid = []
        while (self.nsteps < nsteps):
            self.step(dxmax)            
            if self.fmax < fcutoff:    
                break

            #if the same dislocation is always the maxf dislocation, then if it changes direction each iteration, cut dxmax by dxcut
            maxdislolistx.append(self.maxfdislo.x)
            maxdislolistid.append(self.maxfdislo.id)
            if len(maxdislolistx) == neval: #only perform evaluation every neval
                if all(maxdislolistid[0] == item for item in maxdislolistid): #if all ids are the same, count number of sign flips
                    if abs(maxdislolistx[1]) > abs(maxdislolistx[0]): 
                        inc = True
                    else:
                        inc = False
                    iflip = 0                        
                    for i in xrange(1,neval-1):
                        if abs(maxdislolistx[i+1]) > abs(maxdislolistx[i]): #increasing
                            if inc == True:
                                pass
                            else:
                                iflip = iflip + 1
                            inc = True
                        else: #decreasing
                            if inc == True:
                                iflip = iflip + 1
                            else:
                                pass
                            inc = False
                    if iflip >= neval-2: #sign is flipping each evaluation, cut dxmax
                        print 'cutting step'
                        dxmax = dxmax/dxcut
                maxdislolistx = []
                maxdislolistid = []        

    def stepton_vcriterion(self, nsteps, dxmax, vcutoff):
        while (self.nsteps < nsteps):
            tstart = self.t
            self.step(dxmax)
            tend = self.t
            dt = tend-tsart
            v = dxmax/dt
            if v < vcutoff:    
                break

    def stepton_dtcriterion(self, nsteps, dxmax, dtcutoff):
        while (self.nsteps < nsteps):
            tstart = self.t
            self.step(dxmax)
            tend = self.t
            dt = tend-tsart
            if dt < dtcutoff:    
                break

    def calc_width(self, big=1e8):

        xmax, xmin = -big, big
        for dislo in self.dislolist:
            if dislo.x < xmin:
                xmin = dislo.x
            if dislo.x > xmax:
                xmax = dislo.x
        width = xmax - xmin
        if self.horsym: width = xmax*2.0
        
        return width

    def settau(self, tau):
        self.tau = tau       
    
    def reset(self):
        self.t = 0.0
        self.nsteps = 0
        self.maxfdislo = None     

    def calc_stresses_point(self, X, Z):
        sigxx, sigzz, sigxz = 0.0, 0.0, 0.0
        for dislo in self.getalldislo():
            sigxxadd, sigzzadd, sigxzadd = sigstat(X-dislo.x, Z-dislo.z, dislo.ang, dislo.burg, self.matprops)
            sigxx, sigzz, sigxz = sigxx + sigxxadd, sigzz + sigzzadd, sigxz + sigxzadd
        return sigxx, sigzz, sigxz+self.tau

    def calc_stresses_line(self, X, Z):
        sigxx, sigzz, sigxz = 0.0, 0.0, 0.0
        for dislo in self.getalldislo():
            sigxxadd, sigzzadd, sigxzadd = sigstat(X-dislo.x, Z-dislo.z, dislo.ang, dislo.burg, self.matprops)
            sigxx, sigzz, sigxz = sigxx + sigxxadd, sigzz + sigzzadd, sigxz + sigxzadd
        return sigxx, sigzz, sigxz+self.tau

    def calc_stresses_square(self, dims, res):
        xmin, xmax, zmin, zmax = dims[0], dims[1], dims[2], dims[3]
        x = np.linspace(xmin, xmax, res)
        z = np.linspace(zmin, zmax, res)
        X,Z = np.meshgrid(x,z)

        sigxx, sigzz, sigxz = np.zeros((res, res)), np.zeros((res, res)), np.zeros((res, res))
        for dislo in self.getalldislo():
            sigxxadd, sigzzadd, sigxzadd = sigstat(X-dislo.x, Z-dislo.z, dislo.ang, dislo.burg, self.matprops)             
            sigxx, sigzz, sigxz = sigxx + sigxxadd, sigzz + sigzzadd, sigxz + sigxzadd
        return X, Z, sigxx, sigzz, sigxz+self.tau

    def plotsig(self, dims, res=400, pmax=0.1, plevels=61,  \
                cbtitle=None, disloscale=None, plotres=300, fname='test.jpg', justxz=True):

        xmin, xmax, zmin, zmax = dims[0], dims[1], dims[2], dims[3]    
        X, Z, sigxx, sigzz, sigxz = self.calc_stresses_square(dims, res)
        
        #make grid for contour plot
        plotLevels = np.linspace(-pmax, pmax, plevels)
        thiscm = cm.bwr
        thiscm.set_over('red')
        thiscm.set_under('blue')
        
        #correct file format
        fstart = fname[:-4]
        fend = fname[-4:]

        #sigxz
        f = p.figure(figsize=(6,6))
        p.contourf(X, Z, sigxz, levels=plotLevels, cmap=thiscm, extend='both')
        if disloscale: self.plotalldislocations(dims, disloscale=disloscale)
        p.axis([xmin, xmax, zmin, zmax])
        f.savefig(fstart+'-sigxz'+fend, dpi=plotres)
        if cbtitle:
            drawColorbar(pmax, fstart+'-cb'+fend, plevels=plevels, plotres=500, title=cbtitle)
        p.close(f)            
        
        if not justxz:
            #sigxx
            f = p.figure(figsize=(6,6))
            p.contourf(X, Z, sigxx, levels=plotLevels, cmap=thiscm, extend='both')
            if disloscale: self.plotalldislocations(dims, disloscale=disloscale)
            p.axis([xmin, xmax, zmin, zmax])
            f.savefig(fstart+'-sigxx'+fend, dpi=plotres)
            p.close(f)            
            
            #sigzz
            f = p.figure(figsize=(6,6))
            p.contourf(X, Z, sigzz, levels=plotLevels, cmap=thiscm, extend='both')
            if disloscale: self.plotalldislocations(dims, disloscale=disloscale)
            p.axis([xmin, xmax, zmin, zmax])
            f.savefig(fstart+'-sigzz'+fend, dpi=plotres)
            p.close(f) 


    def plotsigline(self, xline, z, pmax=0.1, plotres=300, fname='test.jpg', justxz=True):

        sigxx, sigzz, sigxz = self.calc_stresses_line(xline, z)

        #correct file format
        fstart = fname[:-4]
        fend = fname[-4:]

        #sigxz
        f = p.figure(figsize=(6,6))
        p.plot(xline, sigxz)
        p.axis([xline[0], xline[-1], -pmax, pmax])
        f.savefig(fstart+'-sigxz'+fend, dpi=plotres)
        p.close(f)            

        if not justxz:
            #sigxx
            f = p.figure(figsize=(6,6))
            p.plot(xline, sigxx)
            p.axis([xline[0], xline[-1], -pmax, pmax])
            f.savefig(fstart+'-sigxx'+fend, dpi=plotres)
            p.close(f)     
            
            #sigzz
            f = p.figure(figsize=(6,6))
            p.plot(xline, sigzz)
            p.axis([xline[0], xline[-1], -pmax, pmax])
            f.savefig(fstart+'-sigzz'+fend, dpi=plotres)
            p.close(f)  

    def plotalldislocations(self, dims, disloscale=2.0, dislocolor='0.0'):
        xmin, xmax, zmin, zmax = dims[0], dims[1], dims[2], dims[3]            
        for dislo in self.getalldislo():
            xpos = dislo.x
            zpos = dislo.z
            ang = dislo.ang
            xlen = disloscale*(xmax-xmin)/100.0
            zlen = disloscale*(zmax-zmin)/100.0
            
            #top left, then bot right
            x1, z1 = xpos-xlen/2.0*np.cos(ang), zpos-zlen/2.0*np.sin(ang)
            x2, z2 = xpos+xlen/2.0*np.cos(ang), zpos+zlen/2.0*np.sin(ang)
            p.plot([x1,x2],[z1,z2], color=dislocolor, linewidth=1*disloscale)
            
            x3, z3 = xpos-xlen*np.sin(ang), zpos+zlen*np.cos(ang)
            p.plot([xpos,x3],[zpos,z3], color=dislocolor, linewidth=1*disloscale)
        return p

    def compareAnalytical(self, savefile, b=0.245e-4, d = 1.9e-4, disloscale=1.0, npts=1000, plotres=500, width=None, plotlegend=False):
        if not width: width = self.calc_width()
        hw = width/2.0
        x = np.linspace(-hw, hw, npts)
        
        
        
        if self.tau == None:
            raise AttributeError('compareAnalytical only works if tau is defined, exiting')

        a = analyticalThickness(self.matprops, b, d)
        xthin,tthin = a.thinPts(x, self.tau)
        tthickez = a.thickConstEZ(width, self.tau)/2.0
        tthickhard = a.thickConstHard(width, self.tau)
        maxz = max(tthickhard, tthickez)

        f = p.figure(figsize=(6,6))

        pp = self.plotalldislocations([-hw, hw, -maxz, maxz], disloscale=disloscale)
        line1, = pp.plot(xthin, tthin, 'r-')
        line2, = pp.plot([-hw, hw], [tthickhard, tthickhard], 'k--')
        pp.plot([-hw, hw], [-tthickhard, -tthickhard], 'k--')
        line3, = pp.plot([-hw, hw], [tthickez, tthickez], 'k-.')
        pp.plot([-hw, hw], [-tthickez, -tthickez], 'k-.')
        line4, = pp.plot([hw, hw], [-maxz*2, maxz*2], 'k-')
        pp.plot([-hw, -hw], [-maxz*2, maxz*2], 'k-')
        
        ax = pp.gca()
        ax.set_ylim([-1.2*maxz,1.2*maxz])

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        
        f.savefig(savefile, dpi=plotres)        

        if plotlegend:
            figlegend = p.figure(figsize=(6.0,0.5))
            figlegend.legend([line1, line2, line3, line4], ['Eqn 2', 'Eqn 4', 'Eqn 5', 'Grain boundary'], 'center', numpoints=1, ncol=4, frameon=False, prop={'size':12})
            figlegend.savefig('paper/legend.jpg', dpi=600)

    def drawColorbar(pmax, figname, plevels=31, plotres=500, cmap=cm.bwr, title=None):
        f = p.figure(figsize=(6.5, 1.0))
        ax = f.add_axes([0.05, 0.5, 0.9, 0.15])

        # Set the colormap and norm to correspond to the data for which
        # the colorbar will be used.

        norm = mpl.colors.Normalize(vmin=-pmax, vmax=pmax)
        bounds = np.linspace(-pmax, pmax, plevels)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                        norm=norm,
                                        # to use 'extend', you must
                                        # specify two extra boundaries:
                                        boundaries=bounds,
                                        ticks=[bounds[0], 0.0, bounds[-1]],  # optional
                                        spacing='proportional',
                                        orientation='horizontal')
        if title: cb.set_label(title)

        f.savefig(figname, dpi=1000)            
            
if __name__ == "__main__":

    
    aMg = 0.321e-3
    b1012Mg, d1012Mg = 0.088*aMg, 0.633*aMg #d is vertical spacing bw dislocations
    mgProps = IsoMaterialProps(24.1333e3, 16.9e3, 1.74e-3)

        
