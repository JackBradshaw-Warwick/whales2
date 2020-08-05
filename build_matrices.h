#ifndef BUILD_MATRICES_H
#define BUILD_MATRICES_H

#include "constants.h"

int build_matrices(Mat mats[] , geom_shape geom , equil_fields equil , matrix_coeffs *coeffs);

void build_gq_mod(gq_mod *gq_mod,geom_shape geom, equil_fields eq);
int solve_weights(double weights[],double upp, double low, int num_quad, int num_div);

void build_submatrices(Mat mats[],equil_fields eq,geom_shape* geom,matrix_coeffs* coeffs,gq_mod* gq_mod,int main_mod,int sec_mod,int row_offset,int col_offset);
void build_coeffs(equil_fields eq,geom_shape geom,matrix_coeffs* coeffs,int main_mod,int sec_mod, std::string which_mat);
void fill_submatrix(Mat M,geom_shape *geom, gq_mod* gq_mod,double rad_interp[],std::complex<double> matrix_coeffs[],int main_indices[],int sec_indices[],int dim_sub,int row_offset,int col_offset,bool dw_dh[],InsertMode insert_mode);
std::complex<double> integrate1d(std::complex<double> func_grid[],double rad_interp[],double (*shape[])(double,double,double),int N_psi,int num_quad,int upper_point,int lower_point);
std::complex<double> integrate1d(std::complex<double> func_grid[],gq_mod* gq_mod,double rad_interp[],void (*shape[])(std::vector<double>*,double,double),int N_psi,int num_quad,int upper_point,int lower_point);
std::vector<double> combine_pols(std::vector<double> shape1,std::vector<double> shape2);
void select_shape(std::string* out_shape,std::string shape_order,int shape_num,bool is_perp);
void fill_shapes(void (*shape[])(std::vector<double>*,double,double),bool deriv,std::string shape_order);
void fill_shapes(double (*shape[])(double,double,double),bool deriv,std::string shape_order);
void calc_fourier_sym(double grid_func[],std::complex<double> fourier_func[],int N_psi,int N_theta);
void calc_fourier_full(std::complex<double> grid_func[],std::complex<double> fourier_func[],int N_psi,int N_theta);

#endif
