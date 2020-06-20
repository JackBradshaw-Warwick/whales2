#ifndef EQUILIBRIUM_FUNCS_H
#define EQUILIBRIUM_FUNCS_H

#include "constants.h"

double pol_flux_theta(double r);
double pol_flux_screw(double r);
double safety(double r);
double maj_rad(double psi,double theta);
double height(double psi,double theta);
double f_psi(double psi, double theta);
double dens(double psi, double psi_max, geom_shape geom);
double pres(double psi, double theta);


#endif
