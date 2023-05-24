#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""Polarization analysis as described in https://arxiv.org/pdf/1409.6214.pdf
"""

from __future__ import print_function, division

import numbers

import numpy
import astropy
import scipy

from ixpeobssim.irf.ebounds import PI_ENERGY_MIN, PI_ENERGY_MAX
from ixpeobssim.utils.misc import pairwise
from ixpeobssim.utils.logging_ import logger, abort


# pylint: disable=invalid-name, too-many-arguments, too-many-locals, too-few-public-methods


class xModulationAnalysis:

    """Small class implements the polarization analysis based on the
    event-by-event Stokes parameters described in Kislat et al. 2015,
    see https://arxiv.org/pdf/1409.6214.pdf

    Arguments
    ---------
    phi : array_like
        The array of photoelectron azimuthal directions.
    """

    def __init__(self, phi, weights=None):
        """Constructor.
        """
        # Initialize the weights to one, if no weights are passed.
        if weights is None:
            weights = numpy.full(phi.shape, 1.)
        self._w = weights
        self._q = xStokesAnalysis.stokes_q(phi, weights)
        self._u = xStokesAnalysis.stokes_u(phi, weights)

    def calculate_modulation(self, degrees=False):
        """Calculate the modulation.
        """
        I = numpy.sum(self._w)
        Q = numpy.sum(self._q)
        U = numpy.sum(self._u)
        mu = numpy.full(I.shape, 1.)
        W2 = numpy.sum(self._w**2.)
        return xStokesAnalysis.calculate_polarization(I, Q, U, mu, W2, degrees)



class xStokesAnalysis:

    """This class implements the polarization analysis based on the
    event-by-event Stokes parameters described in Kislat et al. 2015,
    see https://arxiv.org/pdf/1409.6214.pdf

    The basic idea is that we pass to the constructor a list of photoelectron
    direction and energy arrays, along with the necessary response functions
    (and optional weights), and the proper functional combination are
    turned into the correponding (weighted) event-by-event Stokes parameters,
    that can then be easily summed in the proper energy range and turned into
    polarization degree and error.

    Note that (modulo a factor of 2) all the book-keeping and processing is
    done in terms of the reconstructed Stokes parameters defined at the end of
    section 3 of the paper, i.e., the Stokes parameters are divided by the
    modulation factor on an event-by-event basis in the constructor.

    Arguments
    ---------
    q : array_like
        The array of event-by-event Q Stokes parameters.

    u : array_like
        The array of event-by-event U Stokes parameters.

    energy : array_like
        The array of the event energies (in KeV)---must have the same shape as q and u.

    modf : xModulationFactor instance
        The modulation factor---this is evaluated at the event energies and used
        to normalize q and u.

    aeff : xEffectiveArea instance
        The effective area---this is used for the energy flux calculation, and
        and to calculate the acceptance correction, if necessary.

    livetime : float
        The livetime (in s) for the energy flux calculation.

    weights : array_like
        Additional (optional) multiplicative event weights---must have the same
        shape as q and u.

    acceptcorr : bool
        If True, the Stokes parameters are weighted by the inverse of the acceptance.
    """

    def __init__(self, q, u, energy, modf, aeff, livetime, weights=None, acceptcorr=True):
        """Constructor.
        """
        # After https://bitbucket.org/ixpesw/ixpeobssim/issues/539 we need to
        # cache a mask for the events that we keep, in order to be able to
        # do the same filter on "friend" columns down the road, e.g., when
        # doing polarization maps using X and Y.
        self._filter_mask = numpy.ones(energy.shape, dtype=bool)

        # Preliminary check for unphysical Q or U values. Note this is done right
        # at the beginning, prior to the energy filtering, so that we can
        # potentially give useful indications on the row number of the corrupted
        # event(s) on file.
        # This was added in response to
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/539
        mask = numpy.logical_or(numpy.isnan(q), numpy.isnan(u))
        num_events = mask.sum()
        if num_events > 0:
            rows, = numpy.where(mask)
            logger.error('Unphysical Stokes parameters at row(s): %s', rows + 1)
            logger.error('Q = %s, U = %s', q[rows], u[rows])
            logger.error('Filtering out %d corrupted event(s)...', num_events)
            logger.error('This might get you up and running, but should not happen!!!')
            self._filter_mask *= numpy.logical_not(mask)

        # Since the modulation factor is only defined between 1 and 12 keV, and
        # potentially we have events much above 12 keV we need this filtering
        # stage---this kind of polarization analysis makes only sense where the
        # modulation factor is defined or can be realiably extrapolated.
        # This was added in response to
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/539
        # Note we are fairly liberal an take everything from 0 to 15 keV.
        mask = numpy.logical_and(energy >= PI_ENERGY_MIN, energy <= PI_ENERGY_MAX)
        num_events = numpy.logical_not(mask).sum()
        if num_events > 0:
            logger.warning('Filtering out %d event(s) outside the %.2f--%.2f keV energy range',
                num_events, PI_ENERGY_MIN, PI_ENERGY_MAX)
            self._filter_mask *= mask

        # And now we should be good to go!
        self._energy = energy[self._filter_mask]
        self._mu = modf(self._energy)
        self._aeff = aeff(self._energy)
        self._livetime = livetime
        # Initialize the weights to one, if no weights are passed.
        if weights is None:
            weights = numpy.full(self._energy.shape, 1.)
        else:
            weights = weights[self._filter_mask]
        # If the effective area is available, correct for the acceptance.
        if acceptcorr:
            weights /= self._aeff
        # Cache the weights for later use---note the copy() call, here, as we
        # shall modify the other array later.
        self._w = weights.copy()
        # For the Q and U Stokes parameters, we do want to divide by the
        # modulation factor evaluated on an event-by-event basis.
        weights /= self._mu
        self._q = q[self._filter_mask] * weights
        self._u = u[self._filter_mask] * weights

    @property
    def n(self):
        """Simple property function to make the accumulation of the counts easy
        to read.
        """
        return numpy.full(self._energy.shape, 1.)

    @staticmethod
    def stokes_i(phi, weights=None):
        """Convert the event azimuthal angle to the I Stokes parameter, see
        equations (9a) and (A.1a) in Kislat et al., 2015.

        Arguments
        ---------
        phi : array_like
            The array of azimuthal angles.

        weights : array_like
            Optional event weights.
        """
        i = numpy.full(phi.shape, 1.)
        if weights is not None:
            i *= weights
        return i

    @staticmethod
    def stokes_q(phi, weights=None):
        """Convert the event azimuthal angle to the Q Stokes parameter, see
        equations (9b) and (A.1b) in Kislat et al., 2015.

        Note that, compared to equation (9b), we have an extra factor of 2,
        here, that renders the calculations downstream more natural. The rest
        of the class is implemented consistently.

        This is factored out in a staticmethod so that it can be reused
        consistently in other places.

        Arguments
        ---------
        phi : array_like
            The array of azimuthal angles.

        weights : array_like
            Optional event weights.
        """
        q = 2. * numpy.cos(2. * phi)
        if weights is not None:
            q *= weights
        return q

    @staticmethod
    def stokes_u(phi, weights):
        """Convert the event azimuthal angle to the U Stokes parameter, see
        equations (9c) and (A.1c) in Kislat et al., 2015.

        Note that, compared to equation (9c), we have an extra factor of 2,
        here, that renders the calculations downstream more natural. The rest
        of the class is implemented consistently.

        This is factored out in a staticmethod so that it can be reused
        consistently in other places.

        Arguments
        ---------
        phi : array_like
            The array of azimuthal angles.

        weights : array_like
            Optional event weights.
        """
        u = 2. * numpy.sin(2. * phi)
        if weights is not None:
            u *= weights
        return u

    def _energy_mask(self, emin, emax):
        """Return the proper mask to select events in a given energy range.

        Warning
        -------
        I am not quite sure about the binary operands, here. The numpy
        documentation seems to suggest that we should use >= emin and
        < emax, but in that case the unit tests are failing due to numerical
        roundings.

        https://numpy.org/devdocs/reference/generated/numpy.histogram.html
        All but the last (righthand-most) bin is half-open. In other words,
        if bins is:

        [1, 2, 3, 4]

        then the first bin is [1, 2) (including 1, but excluding 2) and the
        second [2, 3). The last bin, however, is [3, 4], which includes 4.
        """
        return numpy.logical_and(self._energy > emin, self._energy <= emax)

    def _weighted_average(self, values, mask):
        """Calculate the weighted average of a given quantity over the input
        weights.
        """
        return numpy.sum(values[mask] * self._w[mask]) / numpy.sum(self._w[mask])

    def _average_energy(self, mask):
        """Return the average energy over a given event mask.
        """
        return self._weighted_average(self._energy, mask)

    def average_energy(self, emin, emax):
        """Return the average energy over a given energy range.
        """
        return self._average_energy(self._energy_mask(emin, emax))

    def _effective_mu(self, mask):
        """Return the effective modulation factor weighted over the input events.
        """
        return self._weighted_average(self._mu, mask)

    def effective_mu(self, emin, emax):
        """Return the effective modulation factor weighted over the input events.
        """
        return self._effective_mu(self._energy_mask(emin, emax))

    def _sum_stokes_parameters(self, mask):
        """Sum the reconstructed Stokes parameters over an event mask.
        """
        I = numpy.sum(self._w[mask])
        Q = numpy.sum(self._q[mask])
        U = numpy.sum(self._u[mask])
        return I, Q, U

    def sum_stokes_parameters(self, emin, emax):
        """Sum the Stokes parameters into a given energy range.
        """
        return self._sum_stokes_parameters(self._energy_mask(emin, emax))

    def W2(self, mask):
        """Return the sum of the squares of weights.

        For an un-weighted analysis, and if the acceptance correction is not
        applied, the weights are all equal to unity, and this sum reduces to
        the number of events passing the cut, which is in turn equal to the
        sum of the I Stokes parameter.
        """
        return numpy.sum(self._w[mask]**2.)

    @staticmethod
    def _check_polarization_input(I, Q, U):
        """Miminal check on the inputs to the calculate_polarization() hook.

        This was added while debugging
        https://bitbucket.org/ixpesw/ixpeobssim/issues/539/
        """
        error_flag = False
        msg = 'Passing invalid %s Stokes parameter to polarization analysis: %s'
        if not numpy.isfinite(I).all():
            logger.error(msg, ('I', I))
            error_flag = True
        if not numpy.isfinite(Q).all():
            logger.error(msg, ('Q', Q))
            error_flag = True
        if not numpy.isfinite(U).all():
            logger.error(msg, ('U', U))
            error_flag = True
        return error_flag

    @staticmethod
    def calculate_polarization(I, Q, U, mu, W2=None, degrees=False):
        """Calculate the polarization degree and angle, with the associated
        uncertainties, for a given q and u.

        This implements equations (21), (36), (22) and (37) in the paper,
        respectively.

        Note that the Stokes parameters passed as the input arguments are assumed
        to be normalized to the modulation factor (for Q and U) on an
        event-by-event basis and summed over the proper energy range.

        Great part of the logic is meant to avoid runtime zero-division errors.
        """
        if xStokesAnalysis._check_polarization_input(I, Q, U):
            abort('Invalid input to xStokesAnalysis.calculate_polarization()')
        # If W2 is not, i.e, we are not passing the sum of weights, we assume
        # that the analysis is un-weighted, and the acceptance correction is
        # not applied, in which case W2 = I and the scale for the errors is 1.
        if W2 is None:
            W2 = I
        # Initialize the output arrays.
        err_scale = numpy.full(I.shape, 1.)
        pd = numpy.full(I.shape, 0.)
        pd_err = numpy.full(I.shape, 0.)
        pa = numpy.full(I.shape, 0.)
        pa_err = numpy.full(I.shape, 0.)
        # Define the basic mask---we are only overriding the values for the array
        # elements that pass the underlying selection.
        # Note we need I > 1., and not simply I > 0., to avoid any possible
        # zero-division runtime error in the calculations, including the error
        # propagation.
        mask = I > 1.
        # First pass at the polarization degree, which is needed to compute the
        # modulation, which is in turn one of the ingredients of the error
        # propagation (remember that Q and U are the reconstructed quantities,
        # i.e., already divided by the modulation factor).
        pd[mask] = numpy.sqrt(Q[mask]**2. + U[mask]**2.) / I[mask]
        # Convert the polarization to modulation---this is needed later for the
        # error propagation.
        m = pd * mu
        # We want the bins to satify the relation (m^2 < 2), since (2 - m^2)
        # is one of the factors of the errors on the polarization.
        mask = numpy.logical_and(mask, m**2. < 2.)
        # We also want to make sure that the modulation factor is nonzero--see
        # formula for the polarization error.
        # It's not entirely clear to me why that would happen, but I assume that
        # if you have a bin with a couple of very-low energy events it is maybe
        # possible?
        mask = numpy.logical_and(mask, mu > 0.)
        # Create a masked version of the necessary arrays.
        _I = I[mask]
        _Q = Q[mask]
        _U = U[mask]
        _W2 = W2[mask]
        _mu = mu[mask]
        _m = m[mask]
        # Second pass on the polarization with the final mask.
        pd[mask] = numpy.sqrt(_Q**2. + _U**2.) / _I
        # See equations (A.4a) and (A.4b), and compare with equations (17a) and
        # (17b) for the origin of the factor sqrt(W2 / I). Also note that a
        # square root is missing in (A.4a) and (A.4b).
        err_scale[mask] = numpy.sqrt(_W2 / _I)
        # Calculate the errors on the polarization degree
        pd_err[mask] = err_scale[mask] * numpy.sqrt((2. - _m**2.) / ((_I - 1.) * _mu**2.))
        assert numpy.isfinite(pd).all()
        assert numpy.isfinite(pd_err).all()
        # And, finally, the polarization angle and fellow uncertainty.
        pa[mask] = 0.5 * numpy.arctan2(_U, _Q)
        pa_err[mask] = err_scale[mask] / (_m * numpy.sqrt(2. * (_I - 1.)))
        assert numpy.isfinite(pa).all()
        assert numpy.isfinite(pa_err).all()
        # Convert to degrees, if needed.
        if degrees:
            pa = numpy.degrees(pa)
            pa_err = numpy.degrees(pa_err)
        return pd, pd_err, pa, pa_err

    def _polarization(self, mask, degrees=False):
        """Return the average polarization degree and angle, along with the
        corresponding statistical uncertainties, over a given event mask.
        """
        I, Q, U = self._sum_stokes_parameters(mask)
        mu = self._effective_mu(mask)
        W2 = self.W2(mask)
        return self.calculate_polarization(I, Q, U, mu, W2, degrees)

    def polarization(self, emin, emax, degrees=False):
        """Return the average polarization degree and angle, along with the
        corresponding statistical uncertainties, into a given energy range.
        """
        return self._polarization(self._energy_mask(emin, emax), degrees)

    @staticmethod
    def calculate_mdp99(mu, I, W2, clip=True):
        """Calculate the MDP based on equation (A.8) in the paper.

        Arguments
        ---------
        mu : array_like
            The effective modulation factor.

        I : array_like
            The I Stokes parameter.

        W2 : array_like
            The sum of the weights squared.

        clip : bool
            If true, the MDP is clipped within the physical bounds 0--100%.
        """
        mdp = numpy.full(I.shape, 1.)
        mask = numpy.logical_and(I > 0., mu > 0.)
        mdp[mask] = 4.29 * numpy.sqrt(W2[mask]) / (mu * I)[mask]
        if clip:
            mdp = numpy.clip(mdp, 0., 1.)
        return mdp

    @staticmethod
    def calculate_n_eff(counts, I, W2):
        """Calculate the effective number of events.
        """
        n_eff = numpy.zeros(I.shape)
        frac_w = numpy.zeros(I.shape)
        mask = I > 0.
        n_eff[mask] = I[mask]**2. / W2[mask]
        if isinstance(counts, numbers.Number):
            frac_w = n_eff / counts
        else:
            frac_w[mask] = n_eff[mask] / counts[mask]
        return n_eff, frac_w

    @staticmethod
    def normalized_stokes_parameters(I, Q, U, I_ERR, Q_ERR, U_ERR):
        """Return the normalized Stokes parameters QN = Q / I and UN = U /I
        properly handling possible zero-division error issues.
        """
        # Initialize the proper arrays with zeros.
        QN = numpy.zeros(I.shape)
        QN_ERR = numpy.zeros(I.shape)
        UN = numpy.zeros(I.shape)
        UN_ERR = numpy.zeros(I.shape)
        # Calculate the central values.
        mask = I > 0.
        QN[mask] = Q[mask] / I[mask]
        UN[mask] = U[mask] / I[mask]
        # Propagate the uncertainties.
        I_REL_ERR = I_ERR / I
        Q_REL_ERR = Q_ERR / Q
        U_REL_ERR = U_ERR / U
        QN_ERR[mask] = abs(QN[mask]) * numpy.sqrt(Q_REL_ERR**2. + I_REL_ERR**2.)
        UN_ERR[mask] = abs(UN[mask]) * numpy.sqrt(U_REL_ERR**2. + I_REL_ERR**2.)
        return QN, UN, QN_ERR, UN_ERR

    @staticmethod
    def stokes_covariance(I, QN, UN, W2):
        """Calculate the covariance between Q and U.
        """
        covariance = numpy.zeros(I.shape)
        mask = I > 0.
        covariance[mask] = W2[mask] / I[mask]**2. * QN * UN
        return covariance

    @staticmethod
    def significance(Q, U, Q_ERR, U_ERR):
        """Calculate the significance of the polarization detection.

        Here we take advantage of the fact that Q^2 / Var(Q) + U^2 / Var(U) is a
        chisquare with two degrees of freedom, a.k.a. an exponential, and we're
        using the corresponding cumulative function.

        The significance in gaussian equivalent sigma is then calculated with
        the ppf of a normal distribution.
        """
        pval = numpy.full(Q.shape, -1.)
        conf = numpy.full(Q.shape, -1.)
        sig = numpy.full(Q.shape ,-1.)
        mask = numpy.logical_and(Q_ERR > 0., U_ERR > 0.)
        chi2 = Q**2. / Q_ERR**2. + U**2. / U_ERR**2.
        pval[mask] = numpy.exp(-0.5 * (chi2[mask]))
        conf[mask] = 1. - pval[mask]
        sig[mask] =  scipy.stats.norm.ppf(conf[mask])
        # And this is obviously ill-conditioned when we're many sigmas from zero,
        # as the pvalue is too small, the detection confidence is 1, and the
        # significance becomes infinite. Under this conditions we might as well
        # take the number of sigma calculated via error propagation on the polarization
        # degree.
        _mask = numpy.isinf(sig)
        sig[_mask] = (Q[_mask]**2. + U[_mask]**2.) /\
            numpy.sqrt(Q[_mask]**2. * Q_ERR[_mask]**2. + U[_mask]**2. * U_ERR[_mask]**2.)
        return pval, conf, sig

    @staticmethod
    def calculate_polarization_sub(I, Q, U, I_ERR, Q_ERR, U_ERR, degrees=False):
        """Calculate the polarization degree and angle, along with the fellow
        associated uncertainties.

        This is using the linear error propagation and is currently only used in
        the polarization cube subtraction.
        """
        pd = numpy.full(I.shape, 0.)
        pd_err = numpy.full(I.shape, 0.)
        pa = numpy.full(I.shape, 0.)
        pa_err = numpy.full(I.shape, 0.)
        mask = I > 0.
        N = (Q[mask]**2. + U[mask]**2)**2.
        pd[mask] = numpy.sqrt(Q[mask]**2. + U[mask]**2) / I
        pd_err[mask] = pd[mask] * numpy.sqrt(
            Q[mask]**2. / N * Q_ERR[mask]**2 + \
            U[mask]**2. / N * U_ERR[mask]**2 + \
            (I_ERR[mask] / I[mask])**2.)
        pa[mask] = 0.5 * numpy.arctan2(U[mask], Q[mask])
        pa_err[mask] = 0.5 * numpy.sqrt(
            (Q[mask]**2. * U_ERR**2. + U[mask]**2. * Q_ERR**2.) / N
        )
        if degrees:
            pa = numpy.degrees(pa)
            pa_err = numpy.degrees(pa_err)
        return pd, pd_err, pa, pa_err

    @staticmethod
    def calculate_stokes_errors(I, Q, U, mu, W2):
        """Calculation of the errors on the Stokes parameters.
        """
        # Initialize the output arrays.
        QN = numpy.zeros(I.shape)
        UN = numpy.zeros(I.shape)
        dQN = numpy.zeros(I.shape)
        dUN = numpy.zeros(I.shape)
        cov = numpy.zeros(I.shape)
        pval = numpy.full(I.shape, -1.)
        conf = numpy.full(I.shape, -1.)
        sig = numpy.full(I.shape ,-1.)
        # Calculate the error on I---this is easy.
        dI = numpy.sqrt(W2)
        # Calculate the normalized Stokes parameters---note the mask to protect
        # from zero division.
        mask = I > 0.
        QN[mask] = Q[mask] / I[mask]
        UN[mask] = U[mask] / I[mask]
        # From this point on, for the errors to be defined, we need a sligthly
        # different mask, and we nee to cast it in such a way we're not dividing
        # by the effective modulation factor, which is zero wherevere there are
        # no events. We define the new mask and cache the values, so that
        # the formulae downstream look easier.
        mask = numpy.logical_and(mask, (QN * mu)**2. <= 2.)
        mask = numpy.logical_and(mask, (UN * mu)**2. <= 2.)
        _W2N = W2[mask] / I[mask]**2.
        _mu = mu[mask]
        _Q = Q[mask]
        _U = U[mask]
        _QN = QN[mask]
        _UN = UN[mask]
        # Propagate the errors on the normalized Stokes parameters. These are
        # equations A.9a and A.9b from Kislat et al. 2015.
        dQN[mask] = numpy.sqrt(_W2N * (2. / _mu**2. - _QN**2.))
        dUN[mask] = numpy.sqrt(_W2N * (2. / _mu**2. - _UN**2.))
        # Estimated covariance between QN and UN, see equation A.10 in Kislat et al. 2015.
        cov[mask] = -_W2N * _QN * _UN
        # Calculate the errors on Q and U
        dQ = I * dQN
        dU = I * dUN
        # Calculate the confidence of the significance. Here we take advantage of
        # the fact that Q^2 / Var(Q) + U^2 / Var(U) is a chisquare with two
        # degrees of freedom, a.k.a. an exponential, and we're using the corresponding
        # cumulative function. The significance in gaussian equivalent sigma
        # is then calculated with the ppf of a normal distribution.
        pval[mask] = numpy.exp(-0.5 * (_Q**2. / dQ[mask]**2. + _U**2. / dU[mask]**2.))
        conf[mask] = 1. - pval[mask]
        sig[mask] =  scipy.stats.norm.ppf(conf[mask])
        # And this is obviously ill-conditioned when we're many sigmas from zero,
        # as the pvalue is too small, the detection confidence is 1, and the
        # significance becomes infinite. Under this conditions we might as well
        # take the number of sigma calculated via error propagation on the polarization
        # degree.
        _mask = numpy.isinf(sig)
        sig[_mask] = (Q[_mask]**2. + U[_mask]**2.) /\
            numpy.sqrt(Q[_mask]**2. * dQ[_mask]**2. + U[_mask]**2. * dU[_mask]**2.)
        return QN, UN, dI, dQ, dU, dQN, dUN, cov, pval, conf, sig

    def polarization_table(self, ebinning, degrees=True):
        """Return a table with all the relevant polarization parameters.

        Note the column names, here, are taken from the definition of the
        FITS file holding the data structure, and we have to be careful in passing
        out the actual values in the right order.

        .. warning::

           It would be nice to be able to grab the column definition from the
           proper file (binning.fmt) but for some reason I do get a circular
           import error, when I try and do that. It might be a sensible thing to
           understand why and do the proper refactoring.
        """
        col_names = \
            ['ENERG_LO', 'ENERG_HI', 'E_MEAN', 'COUNTS', 'MU', 'W2', 'N_EFF', 'FRAC_W', 'MDP_99',
            'I', 'I_ERR', 'Q', 'Q_ERR', 'U', 'U_ERR', 'QN', 'QN_ERR', 'UN', 'UN_ERR', 'QUN_COV',
            'PD', 'PD_ERR', 'PA', 'PA_ERR', 'P_VALUE', 'CONFID', 'SIGNIF']
        table = astropy.table.Table(names=col_names)
        for emin, emax in pairwise(ebinning):
            mask = self._energy_mask(emin, emax)
            # Average energy, effective modulation factor and sum of weights squared.
            emean = self._average_energy(mask)
            mu = self._effective_mu(mask)
            counts = numpy.count_nonzero(mask)
            W2 = self.W2(mask)
            # Stokes parameters and associated uncertainties.
            I, Q, U = self._sum_stokes_parameters(mask)
            QN, UN, dI, dQ, dU, dQN, dUN, cov, pval, conf, sig = \
                self.calculate_stokes_errors(I, Q, U, mu, W2)
            # Effective number of counts, and MDP.
            n_eff, frac_w = self.calculate_n_eff(counts, I, W2)
            mdp = self.calculate_mdp99(mu, I, W2)
            # Polarization analysis.
            pd, pd_err, pa, pa_err = self.calculate_polarization(I, Q, U, mu, W2, degrees)
            # And now assemble each row in the correct order.
            row = (emin, emax, emean, counts, mu, W2, n_eff, frac_w, mdp, I, dI,\
                Q, dQ, U, dU, QN, dQN, UN, dUN, cov, pd, pd_err, pa, pa_err, pval, conf, sig)
            table.add_row(row)
        return table

    def _binned_cube(self, x, y, binning, keys, values):
        """Generic function to compute binned cubes of a given list of
        quantities.

        This is meant to create three-dimensional histograms (in sky-coordinates
        and energy) of a given set of quantities, and returns a dictionary
        with the three-dimensional accumulated arrays, indexed by the keys passed
        as the third argument.

        This is acting as a convenience method to pass suitable data structure
        to the binning routines---to be stored in FITS files.

        Arguments
        ---------
        x, y : array_like
            Pixelized sky-coordinates of the event list (matching the energy array).

        binning : 3-element iterable
            The cube binning in sky-coordinates and energy.

        keys : list or tuple of strings
            The keys for the output dictionary.

        values : list or tuple of arrays (typically class members)
            The actual quantities to be histogrammed.
        """
        x = x[self._filter_mask]
        y = y[self._filter_mask]
        assert x.shape == y.shape == self._energy.shape
        assert len(binning) == 3
        data = numpy.column_stack((self._energy, x, y))
        binned_data = {}
        for key, val in zip(keys, values):
            binned_data[key], _ = numpy.histogramdd(data, weights=val, bins=binning)
        return binned_data

    def mdp_map_cube(self, x, y, binning):
        """Return the necessary cube to create a MDPMAPCUBE binned object.

        Note the extensions names, here, are taken from the definition of the
        FITS file holding the data structure, and we have to be careful in passing
        out the actual values in the right order.

        For reference, here is a snapshot the field, in the appropriate order:
        ['E_MEAN', 'COUNTS', 'MU', 'W2', 'N_EFF', 'FRAC_W', 'MDP_99', 'I']
        """
        keys = ('E_MEAN', 'COUNTS', 'MU', 'I', 'W2')
        values = (self._energy * self._w, self.n, self._mu * self._w, self._w, self._w**2.)
        data = self._binned_cube(x, y, binning, keys, values)
        # Normalize the effective modulation factor to the summed I.
        # Note we have to protect the division for the bins where I == 0.
        counts = data['COUNTS']
        I = data['I']
        W2 = data['W2']
        mask = I > 0.
        data['E_MEAN'][mask] /= I[mask]
        data['MU'][mask] /= I[mask]
        args = [data[key] for key in ('MU', 'I', 'W2')]
        data['MDP_99'] = self.calculate_mdp99(*args)
        data['N_EFF'], data['FRAC_W'] = self.calculate_n_eff(counts, I, W2)
        return data

    def polarization_map_cube(self, x, y, binning):
        """Return the necessary cube to create a PMAPCUBE binned object.
        """
        data = self.mdp_map_cube(x, y, binning)
        keys = ('Q', 'U')
        values = (self._q, self._u)
        data.update(self._binned_cube(x, y, binning, keys, values))
        args = [data[key] for key in ('I', 'Q', 'U', 'MU', 'W2')]
        data['QN'], data['UN'], data['I_ERR'], data['Q_ERR'], data['U_ERR'], data['QN_ERR'], \
            data['UN_ERR'], data['QUN_COV'], data['P_VALUE'], data['CONFID'], data['SIGNIF'] = \
            self.calculate_stokes_errors(*args)
        data['PD'], data['PD_ERR'], data['PA'], data['PA_ERR'] = \
            xStokesAnalysis.calculate_polarization(*args, degrees=True)
        return data
