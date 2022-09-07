.. _implementation:

Implementation details
======================

This section describes the technical details of the simulation flow for a 
:ref:`standard model component <model-components>` and for the conversion of a 
:ref:`Chandra observation <chandra-observation>`. In order to take advantage of 
the efficient array manipulation capabilities provided by numpy, the entire 
implementation is vectorized, i.e. we don't have an explicit event loop in 
python.

Standard model component
------------------------

The basic flow of the simulation for a single standard model component is coded 
in :py:meth:`ixpeobssim.srcmodel.roi.xModelComponentBase.rvs_event_list()`.
Mathematically speaking, the simulation algorithm can be spelled out in the form
of the following basic sequence:

1. Given the source spectrum :math:`\mathcal{S}(E, t)` and the effective area
   :math:`A_{\rm eff}(E)`, we calculate the count spectrum as a function of
   the energy and time:

   .. math::
      \mathcal{C}(E, t) = \mathcal{S}(E, t) \times A_{\rm eff}(E)
      \quad [\text{s}^{-1}~\text{keV}^{-1}].
   
   In case of a source with cosmological redshift :math:`Z` or hydrogen column 
   density :math:`n_H` different from zero, these two quantities are used to 
   correctly re-scale :math:`\mathcal{C}(E, t)`:
   
    .. math::
      \mathcal{C}'(E, t) = \mathcal{C}(E(1 + E), t) \exp{(-n_H \sigma(E))}

2. We calculate the light curve of the model component (in counts space) by
   integrating over the energy:

   .. math::
      \mathcal{L}(t) = \int_{E_{\rm min}}^{E_{\rm max}} \mathcal{C}'(E, t) dE
      \quad [\text{Hz}].

3. We calculate the total number of expected events :math:`N_{\rm exp}` by
   integrating the count rate over the observation time:

   .. math::
      N_{\rm exp} = \int_{t_{\rm min}}^{t_{\rm max}} \mathcal{L}(t) dt

   and we extract the number of *observed* events :math:`N_{\rm obs}` according
   to a Poisson distribution with mean :math:`N_{\rm exp}`.

4. We treat the count rate as a one-dimensional probability density function
   in the random variable :math:`t`, we extract a vector :math:`\hat{t}`
   of :math:`N_{\rm obs}` values of :math:`t` according to this pdf---and we
   sort the vector itself. (Here and in the following we shall use the hat
   to indicate vectors of lenght :math:`N_{\rm obs}`.)

5. We treat the array of count spectra :math:`\hat{\mathcal{C}}(E, \hat{t})`,
   evaluated at the time array :math:`\hat{t}`, as an array of one-dimensional
   pdf objects, from which we extract a corresponding array :math:`\hat{E}` of
   (true) energy values. (In an event-driven formulation this would
   be equivalent to loop over the values :math:`t_i` of the array
   :math:`\hat{t}`, calculate the corresponding count spectrum

   .. math::
      \mathcal{C}_i(E, t_i)

   and treat that as a one-dimensional pdf from which we extract a (true) energy
   value :math:`E_i`, but the vectorized description is more germane to what
   the code is actually doing internally.)

6. We treat the energy dispersion :math:`\hat{\mathcal{D}}(\epsilon; \hat{E})`
   as an array of one-dimensional pdf objects that we use to extract
   the measured energies :math:`\hat{\epsilon}` and the corresponding
   PHA values.

7. We extract suitable arrays of (true) :math:`\hat{\text{RA}}`,
   :math:`\hat{\text{Dec}}` values and, similarly to what we do with the energy
   dispersion, we smear them with the PSF in order to get the correponding
   measured quantities.

8. We use the polarization degree :math:`P` and angle :math:`\alpha` of the
   model component to calculate the visibility :math:`\hat{\xi}` and the phase
   :math:`\hat{\phi}_0` of the azimuthal distribution modulation, given the
   modulation factor :math:`\mu(E)` of the polarimeter

   .. math::
      \hat{\xi} =
      \hat{P}(\hat{E}, \hat{t}, \hat{\text{RA}}, \hat{\text{Dec}}) \times 
      \mu(\hat{E})
      
      \hat{\phi}_0 =
      \hat{\alpha}(\hat{E}, \hat{t}, \hat{\text{RA}}, \hat{\text{Dec}}),

   and we use these values to extract the directions of emission of the
   photoelectron.

(For periodic sources all of the above is done in phase, rather than in time,
and the latter is recovered at the very end using the source ephemeris, but
other than that there is no real difference.)

For source models involving more than one component, this is done for each
component separately, and the different resulting phothon lists are then
merged and ordered in time at the end of the process.

Conversion of a Chandra observation
-----------------------------------

In this case the starting point is a Chandra photon list, publicly available on
the on-line `database <http://cda.harvard.edu/pop/>`_, which includes
temporal, energetic and spatial information for each detected event. 
Being the Chandra angular and energetic resolution (respectively HPD 
:math:`\sim` 1 arcsec and FWHM :math:`\sim` 5% at 5.9 keV) much better
than those expected for the IXPE mission), energies :math:`E_C` and positions
:math:`RA_C`, :math:`Dec_C` measured by Chandra can be assumed as Monte 
Carlo truth. 

From the implementation standpoint, the conversion is coded 
in :py:meth:`ixpeobssim.srcmodel.roi.xChandraObservation.rvs_event_list()` and 
can be summarized in the following steps:

1. We compute the ratio :math:`T_R` between the Chandra observation time 
   :math:`T_C` and the provided IXPE simulation time :math:`T_X`:
   
   .. math::
       T_R = T_C / T_X

2. We scale Chandra detected events by first concatenating the measured arrays 
   (:math:`E_C`, :math:`RA_C` and :math:`Dec_C`) as many times as the integer 
   part of :math:`T_R` (denoted :math:`\lfloor T_R \rfloor`) and then appending 
   to the end of the sequence a portion whose relative length is equal to the 
   fractional part of :math:`T_R` (called :math:`\lbrace T_R \rbrace`). 
   For instance, in case of the energy:
   
   .. math::
       \hat{E} = \underbrace{E_C \: E_C \cdots E_C}_{\lfloor T_R \rfloor} 
       \underbrace{E_C}_{\lbrace T_R \rbrace}.

3. We calculate the effective area ratio :math:`A_R` between one IXPE telescope 
   and Chandra as a function of energy and off-axis angle :math:`\theta`:
   
   .. math::
       A_R(E, \theta) = \dfrac{A_{\mathrm{eff},X} (E,\theta)} 
       {A_{\mathrm{eff},C}(E,\theta)}

4. This ratio allows to conveniently down-sample the Chandra time-scaled photon 
   list. In practice, this is done by throwing an array of random numbers 
   :math:`\hat{r}` between 0 and 1, and accepting as IXPE events only those 
   that meet the condition:
   
   .. math::
       \hat{r} < A_R(\hat{E},\hat{\theta})

5. We randomly extract the event times :math:`\hat{t}` uniformly between the 
   starting and ending observation time.
   
6. From this point, once the arrays :math:`\hat{E}`, :math:`\hat{t}`, 
   :math:`\hat{RA}` and :math:`\hat{Dec}` of true values are selected, 
   the simulation goes ahead as in the standard model component case. 
   In particular, we smear these quantities with the IXPE response functions 
   and we extract the emission angle :math:`\hat{\phi}` based on the input 
   polarization model.
   
In case of definition of multiple sub-regions of the ROI, the conversion is done
separately for each of them (skipping those flagged as to be removed), and the
resulting photon lists are merged and ordered in time at the end of the process.
