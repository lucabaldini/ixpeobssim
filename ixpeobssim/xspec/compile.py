#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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


import subprocess

from ixpeobssim import IXPEOBSSIM_XSPEC


def cmd(expr):
    """Exec a shell command.
    """
    subprocess.call(expr, shell=True)

def compile_():
    """Compile the local models.

    See more information about local models at
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html

    In addition to the subroutine, XSPEC requires a text file describing the
    model and its parameters. The standard models are specified in the model.dat
    file so we usually refer to this text file by that name. A sample model.dat
    entry has the following form:

    modelentry        4  0.    1.e20      modelfunc    add  0 0
    lowT    keV     0.1   0.0808  0.0808 79.9      79.9       0.001
    highT   keV     4.    0.0808  0.0808 79.9      79.9       0.001
    Abundanc " "    1.    0.      0.      5.        5.        0.01
    *redshift " "   0.0

    The first line for each model gives the model name, the number of parameters,
    the low and high energies for which the model is valid, the name of the
    subroutine to be called and the type of model (add, mul, mix, or con, or acn).
    The final two arguments are flags: the first should be set to 1 if model
    variances are calculated by modelfunc and the second should be set to 1 if
    the model should be forced to perform a calculation for each spectrum. This
    final flag is necessary because if multiple spectra have the same energy bins,
    the default behavior is to perform the model calculation for just one
    spectrum and copy the results for each of the others. However, if a model
    depends on information about the spectrum in addition to its energy ranges,
    it must be forced to perform a calculation for each spectrum.

    The remaining lines in the text file specify each parameter in the model.
    For regular model parameters the first two fields are the parameter name
    followed by an optional units label. If there is no units label, then there
    must be a quoted blank (`` '') placeholder. The remaining 6 numerical entries
    are the default parameter value, hard min, soft min, soft max, hard max, and
    fit delta, which are described in the newpar command section.

    There are three special types of parameter which can be used. If the name of
    the parameter is prefixed with a ``*'' the parameter is a ``scale'' parameter
    and cannot be made variable or linked to any kind of parameter other than
    another scale parameter. Since the parameter value can never vary only the
    initial value need be given. If the name of the parameter is prefixed with a
    ``$'' the parameter is a ``switch'' parameter which is not used directly as
    part of the calculation, but switches the model component function's mode of
    operation (i.e. calculate or interpolate). Switch parameters only have 2
    fields: the parameter name and an integer value.

    Finally, if a P is added at the end of the line for a parameter then the
    parameter is defined to be periodic. During a fit, a periodic parameter will
    not be pegged if it tries to exceed its hard limits. Instead it will be
    assigned a value within its limits:
    f(max + delta) = f(min + delta), f(min-delta) = f(max-delta).
    The soft min and max settings are irrelevant for period parameters and will
    be ignored.
    """
    cmd('initpackage ixpeobssim ixpeobssim_model.dat %s' % IXPEOBSSIM_XSPEC)
    cmd('hmake')

def cleanup():
    """Cleanup.
    """
    cmd('hmake clean')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--cleanup', action='store_true',
                        help='cleanup to compiled files')
    args = parser.parse_args()
    if args.cleanup:
        cleanup()
    else:
        compile_()
