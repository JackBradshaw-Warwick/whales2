#ifndef UTILITIES_H
#define UTILITIES_H

#include "constants.h"

void init_equil(equil_fields *equil,int N_psi,int N_interp,int N_theta);
void init_coeffs(matrix_coeffs* coeffs,geom_shape geom);
void init_gq_mod(gq_mod* gq_mod, geom_shape geom, equil_fields eq);
void delete_gq_mod(gq_mod* gq_mod, int N_psi);
void delete_full_geom_shape(geom_shape* geom);
void delete_full_equil_fields(equil_fields* eq);
void delete_full_matrix_coeffs(matrix_coeffs* coeffs);
void assign_gq( const double*& gq_eval , const double*& gq_weigh , int num_quad);
void udsym(equil_fields* eq, geom_shape geom);
int calc_hermitian_diff(Mat M,int row_dim,int col_dim);
void rotate(PetscScalar perp_values[],PetscScalar wedge_values[],int N_r,int m_range);
void sort_eigenvecs(Vec eigenvecs[],PetscScalar xi_perp[],PetscScalar xi_wedge[],geom_shape geom,int keep_indices[], int index);
double newt_raph(double (*func)(double),double guess=0.0,int N_max=100,double tol=1.2e-16);
double expand_bracket(double (*func)(double),double a,double b,double brac_size,double ratio,bool a_fixed=true);
double bisect(double (*func)(double),double a,double b,double target=0.0,double tol=1.2e-16,int N_max=200);
void trim(std::string & s);
void read_in(std::string token_name,std::string & value);
void fill_dim(int &dim_loc,int pol_mode,int N_psi,std::string shape_order);
void fill_indices(int *perp_indices[],int *wedge_indices[],int pol_mode,int N_psi,std::string shape_order);
void fill_2dgrid(double (&foo)(double,double),double grid[],double psi_vals[],int N_psi,double theta_vals[],int N_theta);
void fill_2dgrid(std::complex<double> (&foo)(double,double),std::complex<double> grid[],double psi_vals[],int N_psi,double theta_vals[],int N_theta);
void clean_grid(std::complex<double> grid[],int size);
void clean_grid(double grid[],int size);
void print_grid(double grid[],int N_psi,int N_theta);
void print_grid(std::complex<double> grid[],int N_psi,int N_theta);

#endif
