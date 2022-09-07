PHOTONS extension
-----------------

Header keywords
~~~~~~~~~~~~~~~

* ``TELESCOP``: `Telescope name`
* ``INSTRUME``: `Instrument name`
* ``DETNAM``: `Name of the logical detector unit`
* ``DET_ID``: `Name of the physical detector unit`
* ``TSTART``: `Observation start time in MET`
* ``TSTOP``: `Observation end time in MET`
* ``DATE-OBS``: `Observation start datetime`
* ``DATE-END``: `Observation end datetime`
* ``TELAPSE``: `TSTOP-TSTART`
* ``TIMESYS``: `Time system`
* ``TIMEUNIT``: `Time units`
* ``TIMEREF``: `Time reference`
* ``MJDREFI``: `MJD ref day 01 Jan 2017 00:00:00 UTC`
* ``MJDREFF``: `Frac part of MJD ref (32.184secs+37leapsecs)`
* ``TIMEZERO``: `Zero time`
* ``ONTIME``: `On source time`
* ``LIVETIME``: `On source time corrected for dead time`
* ``DEADC``: `Dead time correction`
* ``DEADAPP``: `Has DEADC been applied to data`
* ``RA_OBJ``: `[deg] R.A. Object`
* ``DEC_OBJ``: `[deg] Dec Object`
* ``RA_PNT``: `[deg] RA pointing`
* ``DEC_PNT``: `[deg] Dec pointing`
* ``OBJECT``: `Name of observed object`

Columns
~~~~~~~

* ``SEC`` (J): `Integral part of event time (MET)` [s]
* ``MICROSEC`` (J): `Fractional part of event time (MET)` [us]
* ``TIME`` (D): `Event time in seconds (MET)` [s]
* ``ENERGY`` (E): `Event energy in keV` [keV]
* ``RA`` (E): `Right Ascension of the photon` [deg]
* ``DEC`` (E): `Declination of the photon` [deg]
* ``DETX`` (E): `X position at the top of the Be window (GPD frame)` [mm]
* ``DETY`` (E): `Y position at the top of the Be window (GPD frame)` [mm]
* ``POL_DEG`` (E): `Polarization degree` []
* ``POL_ANG`` (E): `Polarization angle` [rad]
* ``SRC_ID`` (I): `Monte Carlo source identifier` [None]

