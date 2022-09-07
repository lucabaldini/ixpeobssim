# Preamble

Synch-ed with https://github.com/HEASARC/xspec_localmodels/tree/master/polmodels
on December 17, 2021.

See https://bitbucket.org/ixpesw/ixpeobssim/issues/489/update-the-xspec-polarization-models


# Simple multiplicative models for spectropolarimetry

<p>Simple multiplicative models for spectropolarimetry. These models assume
that the input spectra include XFLT keywords (likely XFLT0001) with
the values 'Stokes:0', 'Stokes:1', 'Stokes:2' for the I, Q, and U Stokes
parameters, respectively. These models are provided as templates which
can be modified to perform more complicated operations.</p>

## Models


### polconst: an energy-independent polarization fraction and angle

The parameters are:
<table>
<tr> <td>1</td> <td>A</td>    <td>polarization fraction</td>        </tr>
<tr> <td>2</td> <td>psi</td>  <td>polarization angle (degrees)</td> </tr>
</table>


### pollin: a linear dependence on energy for both polarization fraction and energy

    A(E)   = A1 + (E-1)*Aslope
    psi(E) = psi1 + (E-1)*psislope

The parameters are:
<table>
<tr> <td>1</td> <td>A1</td>         <td>polarization fraction at 1 keV</td>         </tr>
<tr> <td>2</td> <td>Aslope</td>     <td>polarization fraction slope</td>            </tr>
<tr> <td>3</td> <td>psi1</td>       <td>polarization angle at 1 keV (degrees)</td>  </tr>
<tr> <td>4</td> <td>psislope</td>   <td>polarization angle slope</td>               </tr>
</table>


### polpow: a power-law dependence on energy for both polarization fraction and energy

    A(E)   = Anorm * E^(-Aindex)
    psi(E) = psinorm * E^(-psiindex)

The parameters are:
<table>
<tr> <td>1</td> <td>Anorm</td>    <td>polarization fraction at 1 keV</td>         </tr>
<tr> <td>2</td> <td>Aindex</td>   <td>polarization fraction index</td>            </tr>
<tr> <td>3</td> <td>psinorm</td>  <td>polarization angle at 1 keV (degrees)</td>  </tr>
<tr> <td>4</td> <td>psiindex</td> <td>polarization angle index</td>               </tr>
</table>

