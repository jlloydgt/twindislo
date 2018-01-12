from tdd.twindyn import IsoMaterialProps, Dislocation, Simulation
from math import pi, sin
import os, shutil

#material properties
b, d = 0.245e-4, 1.9e-4
mgProps = IsoMaterialProps.mg()

#twin and sim properties
gbwidth = 1.0 ##
tauexternal = 10.0 ##

#initial iteration properties
dxinit = 0.05 
dxmin = dxinit
nstepsusingdxmin = 3

#iteration properties
nmax = 20000
dxmaxperstep = 1.e-4
fmaxperstep = 1.e-4
dxfactor = 2.0
dxfactor2 = 1.3
ndisloiter = 10

#plotting and saving
destfolder = 'runs/1um/tau10-nosym-dip/' ##
maxstressplot = 50.0 ##
yplot = 20.0e-3
plotdim = [-gbwidth/2.-0.1, gbwidth/2.+0.1, -yplot, yplot]
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
sim = Simulation(tb, mgProps, tau=tauexternal, rwall=gbwidth/2., lwall=-gbwidth/2.)

#run simulation
for i in xrange(ndisloiter):
    if i > nstepsusingdxmin:
        dxmin = False
    
    #positive z dislocation dipole   
    xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor, dxmin=dxmin, poszonly=True)
    if converged: #try smaller increment
        xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor2, dxmin=dxmin, poszonly=True)
        if converged:
            break    
    sim.adddislo(Dislocation(xpos, float(i+1)*d, 0.0, b))
    sim.adddislo(Dislocation(xneg, float(i+1)*d, pi, b))
    sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'a.png')        
    sim.stepton_fcriterion3(nmax, dxmaxperstep, fmaxperstep)

    sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'b.png')
    print 'Step ' + str(i+1) +  ' part 1 used ' +  str(sim.nsteps) + ' iterations'
    logfile.write('Step ' + str(i+1) +  ' part 1 used ' +  str(sim.nsteps) + ' iterations\n')
    
    sim.reset()

    #positive z dislocation dipole   
    xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor, dxmin=dxmin, negzonly=True)
    if converged: #try smaller increment
        xneg, xpos, converged = sim.xinsertnegpos(dx0=dxinit, dxfactor=dxfactor2, dxmin=dxmin, negzonly=True)
        if converged:
            break            
    sim.adddislo(Dislocation(xpos, -float(i+1)*d, 0.0, b))
    sim.adddislo(Dislocation(xneg, -float(i+1)*d, pi, b))
    sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'c.png')
    try:
        sim.stepton_fcriterion3(nmax, dxmaxperstep, fmaxperstep)
    except:
        print 'Simulation threw an error, exiting'
        logfile.write('Simulation threw an error, exiting\n')
        break
    sim.plotsig(plotdim, disloscale=disloscale, res=plotres, pmax=maxstressplot, fname=destfolder+'step'+str(i+1)+'d.png')    
    print 'Step ' + str(i+1) +  ' part 2 used ' +  str(sim.nsteps) + ' iterations'
    logfile.write('Step ' + str(i+1) +  ' part 2 used ' +  str(sim.nsteps) + ' iterations\n')
    
    logfile.flush()
    sim.save(destfolder+'zlog'+str(i+1)+'.p')
    sim.reset()
logfile.close()
