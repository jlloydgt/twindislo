from tdd.twindyn import IsoMaterialProps, Dislocation, Simulation
from math import pi, sin
import os, shutil
import numpy as np


#material properties
b, d = 0.245e-3, 1.9e-3
mgProps = IsoMaterialProps.mg()

for ang in np.linspace(10, 20, 3):

    print ang

    xmidpos = 0.0
    xmidneg = 0.0

    #twin and sim properties
    gbwidth = 1.0 ##
    tauexternal = 50.0 ##
    rang, lang = ang*pi/180., ang*pi/180.
    tauxfer = 200.0
    dxxfer = 10.0e-3 #5nm

    #initial iteration properties
    dxinit = 0.05 
    dxmin = dxinit
    nstepsusingdxmin = 5

    #iteration properties
    nmax = 2000
    dxmaxperstep = 1.e-4
    fmaxperstep = 1.e-4
    dxfactor = 2.0
    dxfactor2 = 1.3
    ndisloiter = 500

    #plotting and saving
    destfolder = 'runs/1um/tau50-taucrit200-xfer-cg-'+str(int(ang))+'/' ##
    maxstressplot = 150.0 ##
    yplot = 100.0e-3 
    plotdim = [-1*gbwidth/2.-0.1, 1*gbwidth/2.+0.1, -yplot, yplot]
    plotres = 200
    disloscale = 2.0

    #make destination, copy source file to destination
    if not os.path.exists(destfolder):
        os.makedirs(destfolder)
    srcscript = os.path.basename(__file__)
    shutil.copy(srcscript, destfolder+srcscript)
    logfile = open(destfolder+'logfile.txt', 'w')

    #make initial twin, set up 
    tb = [Dislocation(gbwidth/2., 0.0, 0.0, b, 'gb'), Dislocation(-gbwidth/2., 0.0, pi, b, 'gb')]
    sim = Simulation(tb, mgProps, tau=tauexternal, rwall=gbwidth/2., lwall=-gbwidth/2., rang=rang, lang=lang, tauxfer = tauxfer, dxxfer = dxxfer)

    #run simulation
    for i in xrange(ndisloiter):
        if i > nstepsusingdxmin:
            dxmin = False
        
        #positive z dislocation dipole   
        xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor, dxmin=dxmin, poszonly=True, dxmax=0.05, xmid=xmidpos)
        if converged: #try smaller increment
            xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor2, dxmin=dxmin, poszonly=True, xmid=xmidpos)
            if converged:
                break        
        xmidpos = (xneg+xpos)/2.0
        sim.adddislo(Dislocation(xpos, float(i+1)*d, 0.0, b))
        sim.adddislo(Dislocation(xneg, float(i+1)*d, pi, b))
        sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'a.png')        
        sim.stepton_fcriterion3(nmax, dxmaxperstep, fmaxperstep,ncutstart=1e8)

        sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'b.png')
        print 'Step ' + str(i+1) +  ' part 1 used ' +  str(sim.nsteps) + ' iterations'
        logfile.write('Step ' + str(i+1) +  ' part 1 used ' +  str(sim.nsteps) + ' iterations\n')
        
        sim.reset()

        #positive z dislocation dipole   
        xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor, dxmin=dxmin, negzonly=True, dxmax=0.05, xmid=xmidneg)
        if converged: #try smaller increment
            xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor2, dxmin=dxmin, negzonly=True, xmid=xmidneg)
            if converged:
                break            
        xmidneg = (xneg+xpos)/2.0
        sim.adddislo(Dislocation(xpos, -float(i+1)*d, 0.0, b))
        sim.adddislo(Dislocation(xneg, -float(i+1)*d, pi, b))
        sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'c.png')
        sim.stepton_fcriterion3(nmax, dxmaxperstep, fmaxperstep,ncutstart=1e8)

        sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'d.png')    
        print 'Step ' + str(i+1) +  ' part 2 used ' +  str(sim.nsteps) + ' iterations'
        logfile.write('Step ' + str(i+1) +  ' part 2 used ' +  str(sim.nsteps) + ' iterations\n')
        
        logfile.flush()
        sim.save(destfolder+'zlog'+str(i+1)+'.p')
        sim.reset()
    logfile.close()
