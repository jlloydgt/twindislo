import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as p

def plotalldislocations(dislolist, xmin, xmax, zmin, zmax, disloscale=2.0):
    
    for dislo in dislolist:
        xpos = dislo.x
        zpos = dislo.z
        ang = dislo.ang
        
        xlen = disloscale*(xmax-xmin)/100.0
        zlen = disloscale*(zmax-zmin)/100.0
        
        #top left, then bot right
        x1, z1 = xpos-xlen/2.0*np.cos(ang), zpos-zlen/2.0*np.sin(ang)
        x2, z2 = xpos+xlen/2.0*np.cos(ang), zpos+zlen/2.0*np.sin(ang)
        p.plot([x1,x2],[z1,z2], color='0.0', linewidth=1*disloscale)
        
        x3, z3 = xpos-xlen*np.sin(ang), zpos+zlen*np.cos(ang)
        p.plot([xpos,x3],[zpos,z3], color='0.0', linewidth=1*disloscale)

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

def plotsigstatic(sim, dims, res=400, pmax=0.1, plevels=61,  \
            cbtitle=None, disloscale=None, plotres=300, fname='test.jpg', justxz=True):

    xmin, xmax, zmin, zmax = dims[0], dims[1], dims[2], dims[3]    
    X, Z, sigxx, sigzz, sigxz = sim.calc_stresses_square(dims, res)

    
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
    if disloscale: plotalldislocations(sim.dislolist, xmin, xmax, zmin, zmax, disloscale=disloscale)
    p.axis([xmin, xmax, zmin, zmax])
    f.savefig(fstart+'-sigxz'+fend, dpi=plotres)
    if cbtitle:
        drawColorbar(pmax, fstart+'-cb'+fend, plevels=plevels, plotres=500, title=cbtitle)
    p.close(f)            
    
    if not justxz:
        #sigxx
        f = p.figure(figsize=(6,6))
        p.contourf(X, Z, sigxx, levels=plotLevels, cmap=thiscm, extend='both')
        if disloscale: plotalldislocations(sim.dislolist, xmin, xmax, zmin, zmax, disloscale=disloscale)
        p.axis([xmin, xmax, zmin, zmax])
        f.savefig(fstart+'-sigxx'+fend, dpi=plotres)
        p.close(f)            
        
        #sigzz
        f = p.figure(figsize=(6,6))
        p.contourf(X, Z, sigzz, levels=plotLevels, cmap=thiscm, extend='both')
        if disloscale: plotalldislocations(sim.dislolist, xmin, xmax, zmin, zmax, disloscale=disloscale)
        p.axis([xmin, xmax, zmin, zmax])
        f.savefig(fstart+'-sigzz'+fend, dpi=plotres)
        p.close(f)            

def plotsigstaticline(dislolist, dims, matprops, pmax,  \
            divby=None, burg=1.0, npts=1000, fname='test.jpg', plotres=500, justxz=True):
    
    xmin, xmax, zmin, zmax = dims[0], dims[1], dims[2], dims[3]
    d = np.sqrt((xmax-xmin)**2+(zmax-zmin)**2)
    plotpts = np.linspace(0,d,npts)
    
    #make grid for contour plot
    X = np.linspace(xmin, xmax, npts)
    Z = np.linspace(zmin, zmax, npts)
    
    #correct file format
    fstart = fname[:-4]
    fend = fname[-4:]

    sigxx, sigzz, sigxz = np.zeros((npts)), np.zeros((npts)), np.zeros((npts))
    for dislo in dislolist:
        sigxxadd, sigzzadd, sigxzadd = sigstat(X-dislo.x, Z-dislo.z, dislo.ang, dislo.burg, matprops)             
        if divby:  #used to be if not divby: divby = dislo.mu/burg
            sigxxadd, sigzzadd, sigxzadd = sigxxadd / divby, sigzzadd / divby, sigxzadd / divby     
        sigxx, sigzz, sigxz = sigxx + sigxxadd, sigzz + sigzzadd, sigxz + sigxzadd

    #sigxz
    f = p.figure(figsize=(6,6))
    if True:
        p.plot(X, sigxz)
        p.axis([X[0], X[-1], -pmax, pmax])
    else:
        p.plot(plotpts, sigxz)
        p.axis([plotpts[0], plotpts[-1], -pmax, pmax])
    f.savefig(fstart+'-sigxz'+fend, dpi=plotres)
    
    if not justxz:
        #sigxx
        f = p.figure(figsize=(6,6))
        if True:
            p.plot(X, sigxx)
            p.axis([X[0], X[-1], -pmax, pmax])
        else:
            p.plot(plotpts, sigxx)
            p.axis([plotpts[0], plotpts[-1], -pmax, pmax])
        f.savefig(fstart+'-sigxx'+fend, dpi=plotres)

        #sigzz
        f = p.figure(figsize=(6,6))
        if True:
            p.plot(X, sigzz)
            p.axis([X[0], X[-1], -pmax, pmax])
        else:
            p.plot(plotpts, sigzz)
            p.axis([plotpts[0], plotpts[-1], -pmax, pmax])
        f.savefig(fstart+'-sigzz'+fend, dpi=plotres)
