/*
  XSPEC headers.

  This file was factored out the other source files to help supporting different
  XSPEC versions, since many header files where moved in XSPEC 12.12.0. See, e.g.,
  https://gcc.gnu.org/onlinedocs/cpp/_005f_005fhas_005finclude.html
  for the documentation about the __has_include directive.

  As noted in the relevant issue,
  https://bitbucket.org/ixpesw/ixpeobssim/issues/472
  selecting the header on the basis of the XSPEC version is less than trivial,
  as XSPEC doesn't seem to export a version in any easily accessible way.
*/


#if __has_include(<FunctionUtility.h>)
  #include <FunctionUtility.h>
#endif
#include <XSFunctions/Utilities/FunctionUtility.h>
