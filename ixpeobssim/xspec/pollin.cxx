// Multiplicative model for polarization modification
// Assumes a polarization fraction and angle with a linear dependence on energy
// parameters:
//   0       A1: polarization fraction at 1 keV
//   1       Aslope: polarization fraction slope
//   2       psi1: polarization angle at 1 keV (degrees)
//   3       psislope: polarization angle slope

//  A(E) = A1 + (E-1)*Aslope
//  psi(E) = psi1 + (E-1)*psislope

#include "xspec_headers.h"

// calcpol is found in calcpol.cxx
void calcpol(const RealArray& energyArray, int spectrumNumber, const RealArray& A,
	     const RealArray& psirad, RealArray& fluxArray);


extern "C"
void pollin(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	    const string& initString)
{
  size_t nF = energyArray.size() - 1;
  fluxArray.resize(nF);
  fluxErrArray.resize(0);

  RealArray A(nF);
  RealArray psirad(nF);

  A = params[0] + (energyArray-1.0)*params[1];
  psirad = (M_PI / 180.) * (params[2] + (energyArray-1.0)*params[3]);

  calcpol(energyArray, spectrumNumber, A, psirad, fluxArray);

  return;
}
