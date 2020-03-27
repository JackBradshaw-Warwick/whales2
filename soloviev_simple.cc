#include "constants.h"

double s_elong, s_triang, s_R_0, s_min_rad, s_inv_asp, s_B_0, s_beta_0, s_alpha;



double sol_radius(double psi,double theta)
{
  return (1.0/s_inv_asp)*(sqrt(1.0+2.0*s_inv_asp*sqrt(psi)*cos(theta)+s_inv_asp*s_inv_asp)-1.0);
}

double sol_height(double psi,double theta)
{
  return s_elong*sqrt(psi)*sin(theta)/sqrt((1.0-0.25*s_inv_asp*s_inv_asp)*(1.0+2.0*s_inv_asp*s_triang*sqrt(psi)*cos(theta)+s_inv_asp*s_inv_asp*s_triang));
}


void fill_sol_simp(equil_fields *equil,geom_shape geom)
{
  std::string value;
  read_in("elong",value);  s_elong=std::stod( value );
  read_in("triang",value);  s_triang=std::stod( value );
  read_in("R_0",value);  s_R_0=std::stod( value );
  read_in("min_rad",value);  s_min_rad=std::stod( value );
  read_in("B_0",value);  s_B_0=std::stod( value );
  read_in("alpha_sol",value);  s_alpha=std::stod( value );
  s_inv_asp=s_min_rad/s_R_0;

  //Ensure f_psi is always sqrt( val > 0.0 )
  if(s_triang < 1.0){ assert(s_alpha>(2.0*s_inv_asp/s_elong)*sqrt((1.0-0.25*s_inv_asp*s_inv_asp)*(1.0-s_triang))); }
  else{ assert(s_alpha>0.0); }

  double sol_A,sol_B,sol_E,sol_F;

  sol_A=2.0*(1.0+(1.0-0.25*s_inv_asp*s_inv_asp)/(s_elong*s_elong));
  sol_B=4.0*s_inv_asp*(1.0+s_triang*(1.0-0.25*s_inv_asp*s_inv_asp)/(s_elong*s_elong));

  s_beta_0=sol_B*s_inv_asp*mu_0/(s_alpha*s_alpha); //s_beta_0 is constrained by condition of zero pressure boundary
  //s_alpha=sqrt(sol_B*s_inv_asp*mu_0*s_B_0/s_beta_0); //S_Alpha is constrained by condition of zero pressure boundary
  double flux_0=s_min_rad*s_min_rad*s_B_0/s_alpha;


  fill_rad( equil, geom, flux_0, flux_0 * rad_low_rat ) ;
  rad_to_psi( equil, geom.N_psi, geom.N_interp );
  
  /*
  //Create 'equally spaced' radial grid (normalised to 1)
  double *psi_norm=new double[geom.N_interp];
  double sqrt_diff=(sqrt(1.0)-sqrt(rad_low_rat))/(geom.N_interp-1);
  for(int iii=0;iii<geom.N_interp;iii++){psi_norm[iii]=(sqrt(rad_low_rat)+iii*sqrt_diff)*(sqrt(rad_low_rat)+iii*sqrt_diff);}
  */


  equil->maj_rad=new double[geom.N_interp*geom.N_theta];
  equil->height=new double[geom.N_interp*geom.N_theta];

  for(int iii=0;iii<geom.N_interp;iii++){
    for(int jjj=0;jjj<geom.N_theta;jjj++){
      equil->maj_rad[iii*geom.N_theta+jjj]=s_R_0+s_min_rad*sol_radius(equil->psi_interp[iii] / flux_0,equil->theta_grid[jjj]);
      equil->height[iii*geom.N_theta+jjj]=s_min_rad*sol_height(equil->psi_interp[iii] / flux_0,equil->theta_grid[jjj]);
    }}



  equil->dens=new double[geom.N_interp*geom.N_theta];
  fill_2dgrid(dens,equil->dens,equil->psi_interp,geom.N_interp,equil->theta_grid,geom.N_theta);



  sol_E=flux_0*(s_R_0/(s_min_rad*s_min_rad*s_min_rad))*(s_inv_asp*sol_A-0.5*sol_B);
  sol_F=flux_0*sol_B/(2.0*s_R_0*s_min_rad*s_min_rad*s_min_rad);
  
  
  equil->pres=new double[geom.N_interp*geom.N_theta];
  equil->f_psi=new double[geom.N_interp*geom.N_theta];

  for(int iii=0;iii<geom.N_interp;iii++){
    for(int jjj=0;jjj<geom.N_theta;jjj++){
      equil->pres[iii*geom.N_theta+jjj]=(0.5*s_beta_0*s_B_0*s_B_0/mu_0)-sol_F*equil->psi_interp[iii];
      equil->f_psi[iii*geom.N_theta+jjj]=sqrt(s_R_0*s_R_0*s_B_0*s_B_0-2.0*sol_E*equil->psi_interp[iii]);
    }}
}
