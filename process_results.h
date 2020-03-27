#ifndef PROCESS_RESULTS_H
#define PROCESS_RESULTS_H

#include "constants.h"

void convert_to_mag(std::complex<double> b_par[],std::complex<double> b_perp[],std::complex<double> b_wedge[],std::complex<double> xi_perp[],std::complex<double> xi_wedge[],equil_fields eq,geom_shape geom,int num_sols,std::string shape_order);
int find_pol_mode(PetscReal polmode_norm[],int m_range,int m_min,Vec eigenvectors[],int particular_val,int dim_locs[]);
int find_rad_mode(double xi_perp[],int length,int eigenmode_number);
int find_mode_number(double xi_perp[],int length,int eigenmode_number);

#endif
