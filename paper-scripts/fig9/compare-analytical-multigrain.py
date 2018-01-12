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

anglist = np.linspace(2.5, 87.5, 35)
print 'noxfer', a.thickConstHard(gbwidth, tauexternal)
for grains in xrange(2,7):
    print 0.0, a.thickConstHard((1.0+2.0*grains)*gbwidth, tauexternal)
    for ang in anglist:
        print ang, a.thickConstHardAngle2(1.0, 50.0, ang*pi/180.0, ngrains=grains)

#a.thickConstHardAngle(1.0, 2.0, 5.0*pi/180.0, 5.0*pi/180.0, savefile='test-orig.jpg')
#a.thickConstHardAngle2(1.0, 2.0, 5.0*pi/180.0, ngrains=3, savefile='test.jpg')
