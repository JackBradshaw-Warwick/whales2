#ifndef INTERP_H
#define INTERP_H

#include "constants.h"

void interpolate_radial(PetscScalar grid_values[],PetscScalar halfway_values[],int gridpoints,int nconv,std::string order);
void interpolate_bilinear(double r_grid[],double r_grid_polate[],int N_r,double theta_grid[],double theta_grid_polate[],int N_theta,PetscScalar orig_vals[],PetscScalar polate_vals[],int nconv);
std::complex<double> interp_lagrange(std::complex<double> grid_func[],double grid[],double int_point,int count_lower,int num_order);
std::complex<double> interpolate_1d(std::complex<double> grid_func[],double grid[],int N_psi,double int_point,std::string order);

#endif
