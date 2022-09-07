import os

from ixpeobssim.srcmodel.roi import xChandraObservation, xChandraROIModel
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.astro import read_ds9
from ixpeobssim import IXPEOBSSIM_CONFIG

# This is the input event file taken from the Chandra database (obs id 8489)
# http://cda.harvard.edu/pop/mainEntry.do
# ftp://cdaftp.cfa.harvard.edu//pub/science/ao08/cat7/8489/primary/acisf08489N003_evt2.fits.gz
EVENT_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG, 'fits', 'cena_evt.fits')

# cena_jet+core contains the definition of jet and core in the Chandra map
REG_SOURCE_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG, 'fits',
                                    'cena_jet+core.reg')

regs = read_ds9(REG_SOURCE_FILE_PATH)

ROI_MODEL = xChandraROIModel(EVENT_FILE_PATH, acis='I')

polarization_degree = constant(0.)
polarization_angle = constant(0.)

# jet region
jet = xChandraObservation('Jet', polarization_degree, polarization_angle, regs[0])
ROI_MODEL.add_source(jet)
# core region
core = xChandraObservation('Core', polarization_degree, polarization_angle, regs[1])
ROI_MODEL.add_source(core)
# the remaining part
core = xChandraObservation('Others', polarization_degree, polarization_angle)
ROI_MODEL.add_source(core)
