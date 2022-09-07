Primary extension
-----------------

Header keywords
~~~~~~~~~~~~~~~

* ``OBS_ID``: `Observation ID`
* ``CONTNUMB``: `Contact number`
* ``OBS_MODE``: `Observation mode`
* ``SRC_CONF``: `Source configuration`
* ``ORIGIN``: `Organization responsible for the data`
* ``CALDBVER``: `CALDB version`
* ``CLOCKCOR``: `Whether the time has been corrected`
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
* ``FILE_LVL``: `File level`
* ``LV2_VER``: `Version of the LV2 data format`


EVENTS extension
----------------

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
* ``EQUINOX``: `equinox of celestial coord system`
* ``RADECSYS``: `celestial coord system`
* ``FILE_LVL``: `File level`
* ``LV2_VER``: `Version of the LV2 data format`

Columns
~~~~~~~

* ``TRG_ID`` (J): `Trigger identifier` [None]
* ``SEC`` (J): `Integral part of event time (MET)` [s]
* ``MICROSEC`` (J): `Fractional part of event time (MET)` [us]
* ``TIME`` (D): `Event time in seconds (MET)` [s]
* ``LIVETIME`` (J): `Live time since the previous event in microseconds` [us]
* ``PHA`` (J): `Event pulse height` [None]
* ``PI`` (E): `Event pulse invariant` [None]
* ``ENERGY`` (E): `Event energy in keV` [keV]
* ``NUM_CLU`` (I): `Number of clusters in the event` []
* ``DETX`` (E): `Reconstructed absorption point X (GPD frame)` [mm]
* ``DETY`` (E): `Reconstructed absorption point Y (GPD frame)` [mm]
* ``RA`` (E): `Event right ascension` [deg]
* ``DEC`` (E): `Event declination` [deg]
* ``X`` (E): `Event X position (SKY frame)` [pixel]
* ``Y`` (E): `Event Y position (SKY frame)` [pixel]
* ``DETPHI`` (E): `Photolectron emission angle (GPD frame)` [rad]
* ``PHI`` (E): `Photoelectron emission angle (SKY frame)` [rad]
* ``Q`` (E): `Corrected event q Stokes parameter` [None]
* ``U`` (E): `Corrected event u Stokes parameter` [None]
* ``W_MOM`` (E): `Event weight` [None]

GTI extension
-------------

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

* ``START`` (D): `GTI start time` [s]
* ``STOP`` (D): `GTI stop time` [s]

MONTE_CARLO extension
---------------------

Header keywords
~~~~~~~~~~~~~~~

* ``IRFNAME``: `name of the IRFs used for the simulation`
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

* ``SRC_ID`` (I): `Monte Carlo source identifier` [None]
* ``MC_ENERGY`` (E): `Monte Carlo event energy` [keV]
* ``MC_PHA`` (J): `Monte Carlo pulse height` [None]
* ``MC_PI`` (E): `Monte Carlo pulse invariant` [None]
* ``MC_RA`` (E): `Monte Carlo right ascension` [degrees]
* ``MC_DEC`` (E): `Monte Carlo declination` [degrees]
* ``MC_X`` (I): `Monte Carlo event X position (SKY frame)` [degrees]
* ``MC_Y`` (I): `Monte Carlo event Y position (SKY frame)` [degrees]
* ``MC_GAIN`` (E): `Relative GEM gain used for the event` [None]

ROITABLE extension
------------------

Header keywords
~~~~~~~~~~~~~~~

* ``ROIRA``: `right ascension of the ROI center`
* ``ROIDEC``: `declination of the ROI center`
* ``EQUINOX``: `equinox for RA and DEC`
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

* ``SRCID`` (I): `source identifier` [None]
* ``SRCNAME`` (A20): `source name` [None]

TIMELINE extension
------------------

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

* ``START`` (D): `Epoch start time` [s]
* ``STOP`` (D): `Epoch stop time` [s]
* ``IN_SAA`` (I): `SAA flag` [None]
* ``TARGET_OCCULT`` (I): `Earth occultation flag` [None]

SC_DATA extension
-----------------

Header keywords
~~~~~~~~~~~~~~~

* ``ROLL``: `spacecraft roll angle`
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

* ``MET`` (D): `Mission elapsed time` [s]
* ``RA_PNT`` (E): `Pointing RA` [degrees]
* ``DEC_PNT`` (E): `Pointing DEC` [degrees]
* ``LAT_GEO`` (E): `Spacecraft latitude` [degrees]
* ``LON_GEO`` (E): `Spacecraft longitude` [degrees]
* ``ALT_GEO`` (E): `Spacecraft elevation` [km]
* ``SUN_ANGLE`` (E): `Angle between the Sun and the target` [degrees]
* ``IN_SAA`` (I): `SAA flag` [None]
* ``TARGET_OCCULT`` (I): `Earth occultation flag` [None]

OCTI extension
--------------

Header keywords
~~~~~~~~~~~~~~~

* ``CALRUNS``: `Number of on-orbit calibration runs`
* ``CALTIME``: `Total time for on-orbit calibration in s`
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

* ``START`` (D): `OCTI start time` [s]
* ``STOP`` (D): `OCTI stop time` [s]

CHRG_MAP extension
------------------

Header keywords
~~~~~~~~~~~~~~~

* ``VERSION``: `Extension version number`
* ``CVSD0001``: `Date when this file should first be used`
* ``CVST0001``: `Time of day when this file should first be used`
* ``NUM_BINS``: `Number of bins per side of the map`

Columns
~~~~~~~

* ``BINX`` (I): `Index for the x coordinate` [None]
* ``BINY`` (I): `Index for the y coordinate` [None]
* ``SLOW`` (D): `Parameters for the charging slow component` [None]
* ``FAST`` (D): `Parameters for the charging fast component` [None]

