import numpy

from ixpeobssim.srcmodel.roi import xPointSource
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant

ra, dec = 45., 45.
spec = power_law(10., 2.)
pol_deg = constant(0.5)
pol_ang = constant(numpy.radians(65.))
src = xPointSource('Point source', ra, dec, spec, pol_deg, pol_ang)

print(src)
