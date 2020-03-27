#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "constants.h"

void deriv_1d(double deriv_grid[],double dep_var[],double ind_var[],int num_points,std::string order);
void deriv_1d(double deriv_grid[],double dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first,std::string order);
void deriv_1d(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_points,std::string order);
void deriv_1d(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first,std::string order);
double ffd_cb(double dep_var[],double diff[],int deriv);
double fd_cb(double dep_var[],double diff[],int deriv);
double cd_cb(double dep_var[],double diff[],int deriv);
double bd_cb(double dep_var[],double diff[],int deriv);
double bbd_cb(double dep_var[],double diff[],int deriv);
void deriv_ang(double deriv_grid[],double dep_var[],double ind_var[],int num_points);
void deriv_ang(double deriv_grid[],double dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first);
void deriv_ang(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_points);
void deriv_ang(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first);
std::complex<double> dpol_xi(std::complex<double> *funcs_fourier[],bool fg_sym[],int main_mod,int sec_mod,int N_psi,int N_theta,int psi_index);
std::complex<double> dpol_sq_xi(std::complex<double> *funcs_fourier[], int main_mod,int sec_mod,int N_psi,int N_theta,int psi_index);
void dpar_xi(double *input_funcs[],std::complex<double> result[],bool fg_sym[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod);
void dpar_sq_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod);
void dwedge_xi(double *input_funcs[],std::complex<double> result[],bool fg_sym[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod);
void dwedge_sq_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod);
void dwedge_dpol_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod);
void dpol_dwedge_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod);

#endif
