#ifndef WHALES_STRUCTS_H
#define WHALES_STRUCTS_H

#include "constants.h"

//Structs
struct geom_shape
{
  // Read-in values

  int N_psi;
  int N_theta;
  int num_quad;
  std::string quad_type;
  
  int m_min;
  int m_max;
  int m_coup;
  double tor_mod;

  std::string shape_order;
  std::string deriv_order;
  std::string interp_order;

  std::string fill_type;
  std::string analytical_type;
  std::string numerical_type;
  std::string input_filepath;

  double min_rad;
  double beta_0;
  double B_0;
  double R_0;
  double elong;
  double triang;
  double alpha_sol;
  double flux_max;
  double solov_A;
  double solov_C;

  int dens_form;
  double A_dens;
  double B_dens;
  double mu;
  double nu;

  std::string output_dir;
  std::string results_filename;

  int write_equil;
  std::string equil_filename;

  int load_mats;
  int write_mats;
  std::string mats_dir;
  std::string mats_filename;

  int shear_on;
  int Hall_on;
  
  // Calculated/program-defined values

  int N_interp;
  int *pol_pos;

  int m_range;

  int dim; //Overall dimensions
  int dims_loc[3]; //Dimension of local block for m=0,|m|=1,|m|>1

  int fourier_size_full;
  int fourier_size_sym;

  int num_var_perp;
  int num_var_wedge;

  int **ind_perp_main; //[num_var_perp][N_psi]
  int **ind_perp_sec; //[num_var_perp][N_psi]
  int **ind_wedge_main; //[num_var_wedge][N_psi]
  int **ind_wedge_sec; //[num_var_wedge][N_psi]

  std::string main_shape; //Temporary shape values
  std::string sec_shape; //Temporary shape values
};

struct equil_fields
{
  double* psi_grid; //[N_psi]
  double* theta_grid; //[N_theta]
  double* rad_var; //[N_psi]
  double* rad_interp; //[N_interp]
  double* psi_interp; //[N_interp]
  double* drad_dpsi; //[N_interp]
  double* maj_rad; //[N_interp*N_theta]
  double* height; //[N_interp*N_theta]

  double* f_psi; //[N_interp*N_theta]
  double* dens; //[N_interp*N_theta]
  double* pres; //[N_interp*N_theta]

  double* dpsi_dR; //[N_interp*N_theta]
  double* dpsi_dZ; //[N_interp*N_theta]
  double* dth_dR; //[N_interp*N_theta]
  double* dth_dZ; //[N_interp*N_theta]

  double* g_pp; //[N_interp*N_theta]
  double* g_pt; //[N_interp*N_theta]
  double* g_tt; //[N_interp*N_theta]
  double* g_phph; //[N_interp*N_theta]	
  double* jacob; //[N_interp*N_theta]
  double* mag_sq; //[N_interp*N_theta]
  double* ones; //[N_interp*N_theta]

  double* j_dot_b; //[N_interp*N_theta]
  double* curv_psi; //[N_interp*N_theta]
  double* neg_shear; //[N_interp*N_theta]
  
  double* cc_0; //[N_interp*N_theta]
  double* cc_1; //[N_interp*N_theta]
  double* cc_2; //[N_interp*N_theta]
  double* cc_3;	//[N_interp*N_theta]
  double* cc_4;	//[N_interp*N_theta]
};

struct matrix_coeffs
{
  std::complex<double>* f_dpdp; //[m_range*m_coup*N_interp]
  std::complex<double>* f_dpp; //[m_range*m_coup*N_interp]
  std::complex<double>* f_pdp; //[m_range*m_coup*N_interp]
  std::complex<double>* f_pp; //[m_range*m_coup*N_interp]
  std::complex<double>* f_dpdw; //[m_range*m_coup*N_interp]
  std::complex<double>* f_dpw; //[m_range*m_coup*N_interp]
  std::complex<double>* f_pw; //[m_range*m_coup*N_interp]
  std::complex<double>* f_wdp; //[m_range*m_coup*N_interp]
  std::complex<double>* f_wp; //[m_range*m_coup*N_interp]
  std::complex<double>* f_wdw; //[m_range*m_coup*N_interp]	
  std::complex<double>* f_ww; //[m_range*m_coup*N_interp]
};

struct gq_mod
{  
  double **weights ; //[N_psi][num_quad]

  int num_div ;
  int wt_gq ;
  double tol ;
};

#endif
