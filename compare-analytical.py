from tdd.twindyn import IsoMaterialProps
from tdd.twindyn import analyticalThickness
from math import pi, sin
import numpy as np

#material properties
b, d = 0.245e-4, 1.9e-4
mgProps = IsoMaterialProps.mg()
a = analyticalThickness(mgProps, b, d)

#loading properties
gbwidth = 1.0 ##
tauexternal = 50.0 ##

t = []
anglist = np.linspace(2.5, 87.5, 35)

#anglist = np.linspace(5.0, 60, 12)
print 'noxfer', a.thickConstHard(gbwidth, tauexternal)
print 0.0, a.thickConstHard(3*gbwidth, tauexternal)
for ang in anglist:
    print ang, a.thickConstHardAngle(gbwidth, tauexternal, ang*pi/180.0, -ang*pi/180.0)


#a.thickConstHardAngle(gbwidth, tauexternal, -anglist[-1]*pi/180.0, -anglist[-1]*pi/180.0, savefile='test.jpg')

