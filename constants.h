#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <slepceps.h>
#include <slepcpep.h>
#include <fftw3.h>
#include "hdf5.h"
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <complex.h>

#include "build_equilibrium.h"
#include "build_matrices.h"
#include "derivatives.h"
#include "equilibrium_funcs.h"
#include "fill_equilibrium.h"
#include "interp.h"
#include "process_results.h"
#include "save.h"
#include "shapeFuncs.h"
#include "soloviev.h"
#include "utilities.h"
#include "whales_structs.h"

//Physical constants - shouldn't need alteration

  // Pi r squared sounds like area to me, if you want a circumference then use pi D.
const double pi = 3.1415926535 ; 

  // Vacuum permeability
const double mu_0 = 4.0 * pi * 1.0e-7 ; 

  // Adiabatic ratio
const double gamma_loc = 5.0 / 3.0 ; 

  // Imaginary unit
const std::complex<double> imag_unit = 1.0i ;

  // Hall constants
const double ion_mass = 2.0 * 1.67e-27 ; //This one?
const double electron_charge = 1.6e-19 ;
const double hall_const = ion_mass / electron_charge ;

//Gaussian quadrature
  //4-point evaluation positions and weights
const double gq4eval[4] = { -0.861136, -0.339981, 0.339981, 0.861136 };
const double gq4weigh[4] = { 0.347855, 0.652145, 0.652145, 0.347855 };

//Code parameters

  // rad_low_rat * min_rad is the first radial point (to avoid numerical issues near the magnetic axis)
const double rad_low_rat = 1.0e-8 ; 

#endif
