#include "constants.h"


double con_tmp=0.60;
double rec=(0.5-con_tmp)+2.0;
double con=1.0/rec;

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

double dens(double psi, double theta)
{
  if(g_counter==0)
    {
  std::string value;
  read_in("B_0",value); g_B_0=std::stod( value);
  read_in("alpha_sol",value); g_alpha=std::stod( value);
  read_in("min_rad",value); g_min_rad=std::stod( value);
  g_counter=1;
    }
  //double flux_max=g_min_rad*g_min_rad*g_B_0/g_alpha;
  double psi_0=1.0;
  double psi_d=0.2;
  double prop=0.5;

  return ion_mass * 1.0e19 ;//* ( 1.0 + prop * exp ( - pow( ( psi - psi_0 ) / psi_d , 2 ) ) );  //pow(1.0 - (psi/1.47)*(psi/1.47),0.4); //flux_max*flux_max*flux_max*flux_max*flux_max*flux_max-0.8*psi*psi*psi*psi*psi*psi;
}

double pres(double psi, double theta)
{
  //std::string value;
  double R_0,B_0,min_rad=1.0; 
  //read_in("R_0",value); R_0=std::stod( value);
  //read_in("B_0",value); B_0=std::stod( value);
  //read_in("min_rad",value); min_rad=std::stod( value);

  return con*con*(min_rad*min_rad-psi*psi)/mu_0;  //2.0*B_0*B_0*(min_rad*min_rad-(2.0*psi/B_0))/(mu_0*R_0*R_0); 
}

double f_psi(double psi, double theta)
{
  std::string value;
  double R_0,B_0; 
  read_in("R_0",value); R_0=std::stod( value);
  read_in("B_0",value); B_0=std::stod( value);

  theta=1.0/0.0; //In case I accidentally put theta dependance in return again
  return -1.0; //B_0*R_0; 
}

/*
double maj_rad(double psi,double theta)
{
  return R_0+(9.0/(8.0*R_0))*(min_rad*min_rad-(2.0*psi/B_phi_0))+sqrt(2.0*psi/B_phi_0)*cos(theta);
}

double height(double psi,double theta)
{
  return sqrt(2.0*psi/B_phi_0)*sin(theta);
}

double dens(double psi, double theta)
{
  double alpha=1.0/3.0;

  return 1.0; //-alpha*sqrt(2.0)*sqrt(psi);
}

double pres(double psi, double theta)
{

  return 2.0*B_phi_0*B_phi_0*(min_rad*min_rad-(2.0*psi/B_phi_0))/(mu_0*R_0*R_0); 
}

double f_psi(double psi, double theta)
{

  theta=1.0/0.0; //In case I accidentally put theta dependance in return again
  return B_phi_0*R_0; 
}
*/
