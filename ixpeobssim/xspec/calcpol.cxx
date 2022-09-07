// This calculates the polarization multiplicative correction for the three Stokes
// parameter spectra for general A(E) and psirad(E). psirad is the polarization
// angle in radians.

// This requires that the input spectral files have an XFLT keyword set
// (likely XFLT0001) to 'Stokes:0', 'Stokes:1', 'Stokes:2' for the I, Q, U cases,
// respectively.

#include "xspec_headers.h"

void calcpol(const RealArray& energyArray, int spectrumNumber, const RealArray& A,
	     const RealArray& psirad, RealArray& fluxArray)
{
  // find out which Stokes parameter this spectrum is from
  // integer code is 0 = I, 1 = Q, 2 = U

  string xname("Stokes");
  Real xvalue = FunctionUtility::getXFLT(spectrumNumber, xname);

  int Stokes = (int)xvalue;
  if ( xvalue == BADVAL ) {
    Stokes = 0;
    FunctionUtility::xsWrite("Failed to read Stokes parameter from XFLTnnnn keyword - applying no polarization correction.", 5);
  }

  fluxArray = 1.0;
  if ( Stokes == 1 ) {
    fluxArray = A * cos (2.0*psirad);
  } else if ( Stokes == 2 ) {
    fluxArray = A * sin (2.0*psirad);
  }

  return;
}
