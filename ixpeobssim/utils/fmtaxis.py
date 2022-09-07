#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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

from __future__ import print_function, division


def axis_label(name, units=None):
    """
    """
    if units is None:
        return name
    return '%s [%s]' % (name, units)


class units:

    """Dumb container class for units.
    """

    flux = 'cm$^{-2}$ s$^{-1}$ keV$^{-1}$'
    norm_cnts = 's$^{-1}$ keV$^{-1}$'



class label:

    """Dumb containter class for axis labels
    """

    energy = axis_label('Energy', 'keV')
    flux = axis_label('Flux', units.flux)
    norm_cnts = axis_label('Normalized counts', units.norm_cnts)
    norm_cnts_I = axis_label('I Normalized counts', units.norm_cnts)
    norm_cnts_U = axis_label('U Normalized counts', units.norm_cnts)
    norm_cnts_Q = axis_label('Q Normalized counts', units.norm_cnts)
    mdp = axis_label('MDP 99\\%')
    phase = axis_label('Pulse phase')
    pl_index = axis_label('Power-law index')
    pl_norm = axis_label('Power-law normalization', units.flux)
    pol_deg = axis_label('Polarization degree')
    pol_ang_deg = axis_label('Polarization angle', 'deg')
    pol_ang_rad = axis_label('Polarization angle', 'rad')
    rate = axis_label('Rate', 'Hz')
    res_sigma = axis_label('Residuals', '$\\sigma$')
    pulls_sigma = axis_label('Pulls', '$\\sigma$')


class fmtaxis:

    """Dumb containter class for axis formats.
    """

    spec = dict(xlabel=label.energy, ylabel=label.flux)
    ene_pol_deg = dict(xlabel=label.energy, ylabel=label.pol_deg)
    ene_pol_ang_rad = dict(xlabel=label.energy, ylabel=label.pol_ang_rad)
    ene_pol_ang_deg = dict(xlabel=label.energy, ylabel=label.pol_ang_deg)
    ene_pol_ang = ene_pol_ang_deg
    norm_cnt_spec = dict(xlabel=label.energy, ylabel=label.norm_cnts)
    norm_cnt_spec_I = dict(xlabel=label.energy, ylabel=label.norm_cnts_I)
    norm_cnt_spec_U = dict(xlabel=label.energy, ylabel=label.norm_cnts_U)
    norm_cnt_spec_Q = dict(xlabel=label.energy, ylabel=label.norm_cnts_Q)
    ene_res_sigma = dict(xlabel=label.energy, ylabel=label.res_sigma)
    ene_pulls_sigma = dict(xlabel=label.energy, ylabel=label.pulls_sigma)
    pp_flux = dict(xlabel=label.phase, ylabel=label.flux)
    pp_pl_index = dict(xlabel=label.phase, ylabel=label.pl_index)
    pp_pl_norm = dict(xlabel=label.phase, ylabel=label.pl_norm)
    pp_pol_deg = dict(xlabel=label.phase, ylabel=label.pol_deg)
    pp_pol_ang_rad = dict(xlabel=label.phase, ylabel=label.pol_ang_rad)
    pp_pol_ang_deg = dict(xlabel=label.phase, ylabel=label.pol_ang_deg)
    pp_pol_ang = pp_pol_ang_deg
