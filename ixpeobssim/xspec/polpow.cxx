// Multiplicative model for polarization modification
// Assumes a polarization fraction and angle with a power-law dependence on energy
// parameters:
//   0       Anorm: polarization fraction at 1 keV
//   1       Aindex: polarization fraction index
//   2       psinorm: polarization angle at 1 keV (degrees)
//   3       psiindex: polarization angle index

//  A(E) = Anorm * E^(-Aindex)
//  psi(E) = psinorm * E^(-psiindex)

#include "xspec_headers.h"

// calcpol is found in calcpol.cxx
void calcpol(const RealArray& energyArray, int spectrumNumber, const RealArray& A,
	     const RealArray& psirad, RealArray& fluxArray);


extern "C"
void polpow(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	    const string& initString)
{
  size_t nF = energyArray.size() - 1;
  fluxArray.resize(nF);
  fluxErrArray.resize(0);

  RealArray A(nF);
  RealArray psirad(nF);

  A = params[0] * pow(energyArray, -params[1]);
  psirad = (params[2] * M_PI / 180.) * pow(energyArray, -params[3]);

  calcpol(energyArray, spectrumNumber, A, psirad, fluxArray);

  return;
}
