// Multiplicative model for polarization modification
// Assumes an energy independent polarization fraction and angle
// parameters:
//   0       A: polarization fraction
//   1       psi: polarization angle (degrees)

#include "xspec_headers.h"

// calcpol is found in calcpol.cxx
void calcpol(const RealArray& energyArray, int spectrumNumber, const RealArray& A,
	     const RealArray& psirad, RealArray& fluxArray);


extern "C"
void polconst(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	      const string& initString)
{
  size_t nF = energyArray.size() - 1;
  fluxArray.resize(nF);
  fluxErrArray.resize(0);

  RealArray A(params[0], nF);
  RealArray psirad(params[1] * M_PI / 180.0, nF);

  calcpol(energyArray, spectrumNumber, A, psirad, fluxArray);

  return;
}
