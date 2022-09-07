
# core region, exclude the core region by passing the exclude=True
core = xChandraObservation('Core', polarization_degree, polarization_angle,
                           regions[1], exclude=True)
ROI_MODEL.add_source(core)
# the remaining part
core = xChandraObservation('Others', polarization_degree, polarization_angle)
ROI_MODEL.add_source(core)

#Add a point source at a given location in the roi.
RA = 201.43
DEC = -43.00
PL_NORM = 0.0004,
PL_INDEX = 2.
spec = power_law(PL_NORM, PL_INDEX)
src = xPointSource('Point source', RA, DEC, spec, polarization_degree, 
                   polarization_angle)
ROI_MODEL.add_source(src)
