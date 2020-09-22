#include "constants.h"


//double con_tmp=0.60;
//double rec=(0.5-con_tmp)+2.0;
double con=2.0;

double g_B_0,g_min_rad,g_alpha;
int g_counter=0;

double pol_flux_theta(double r)
{
  //Never needs to be touched
  
  std::string value;
  double B_0; 
  read_in("B_0",value); B_0=std::stod( value );
  
  return 0.5*B_0*r*r;
}

double pol_flux_screw(double r)
{
  //std::string value;
  //double B_0=1.0; 
  //read_in("B_0",value); B_0=std::stod( value );
  
  return 0.5*con*r*r;
}

double safety(double r)
{
  return 1.0/con;
}

double dens(double psi, double psi_max, geom_shape geom) //psi not r!
{
  if(geom.dens_form == 0){//Default
    return geom.A_dens * pow( 1.0 - geom.B_dens * pow( psi / psi_max , geom.mu ) , geom.nu );
  }
  else if(geom.dens_form == 1){//User defined function (requires compilation)
    return 1.0 ;
  }
  else{
    return 0.0 ;
  }
}

double pres(double r, double theta)
{
  std::string value;
  double B_0=1.0; 
  read_in("B_0",value); B_0=std::stod( value );

  double beta_0=1.0; 
  read_in("beta_0",value); beta_0=std::stod( value );

  return  ( 0.5*B_0*B_0*beta_0 - B_0*B_0*r - con*con*r*r - 0.5*B_0*B_0*r*r ) / mu_0 ; 
}

double f_psi(double r, double theta)
{
  std::string value;
  double B_0=1.0; 
  read_in("B_0",value); B_0=std::stod( value );
  
  theta=1.0/0.0; //In case I accidentally put theta dependance in return again
  return -B_0 * ( 1.0 + r ) ;
}

//Are these even used? Check!
double maj_rad(double psi,double theta)
{
  std::string value;
  double R_0,B_0,min_rad; 
  read_in("R_0",value); R_0=std::stod( value );
  read_in("B_0",value); B_0=std::stod( value );
  read_in("min_rad",value); min_rad=std::stod( value );
  
  return R_0+sqrt(2.0*psi/B_0)*cos(theta);  //(9.0/(8.0*R_0))*(min_rad*min_rad-(2.0*psi/B_0))+sqrt(2.0*psi/B_0)*cos(theta);
}

double height(double psi,double theta)
{
  std::string value;
  double B_0; 
  read_in("B_0",value); B_0=std::stod( value);
  
  return sqrt(2.0*psi/B_0)*sin(theta);
}
