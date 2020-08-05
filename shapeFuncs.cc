/****************************************************************************************/
#include "constants.h"


//Implement shape functions in terms of psi,psi_upp,psi_low. Saves doing the whole coordinate transform


//QUARTIC
double N_A_nhqr(double psi,double upp,double low){return (psi - low) * (psi - low) * (low * low - 4.0 * low * upp - 5.0 * upp * upp + 2.0 * low * psi + 14.0 * upp * psi - 8.0 * psi * psi ) / ( (upp-low)*(upp-low)*(upp-low)*(upp-low) );}
double N_B_nhqr(double psi,double upp,double low){return (psi - upp) * (psi - upp) * ( - 5.0 * low * low - 4.0 * low * upp + upp * upp + 14.0 * low * psi + 2.0 * upp * psi - 8.0 * psi * psi ) / ( (upp-low)*(upp-low)*(upp-low)*(upp-low) );}
double N_C_nhqr(double psi,double upp,double low){return 16.0 * ( psi - low ) * ( psi - low ) * ( psi - upp ) * ( psi - upp ) / ( (upp-low)*(upp-low)*(upp-low)*(upp-low) );}
double N_D_nhqr(double psi,double upp,double low){return ( low + upp - 2.0 * psi ) * ( upp - psi ) * ( low - psi) * ( low - psi) / ( (upp-low)*(upp-low)*(upp-low) );}
double N_E_nhqr(double psi,double upp,double low){return - ( low + upp - 2.0 * psi ) * ( upp - psi ) * ( upp - psi) * ( low - psi) / ( (upp-low)*(upp-low)*(upp-low) );}
double dN_A_nhqr(double psi,double upp,double low){return 2.0 * ( 11.0 * low + 5.0 * upp - 16.0 * psi ) * ( low - psi) * ( upp - psi ) / ( (upp-low)*(upp-low)*(upp-low)*(upp-low) );}
double dN_B_nhqr(double psi,double upp,double low){return 2.0 * ( 5.0 * low + 11.0 * upp - 16.0 * psi ) * ( low - psi) * ( upp - psi ) / ( (upp-low)*(upp-low)*(upp-low)*(upp-low) );}
double dN_C_nhqr(double psi,double upp,double low){return - 32.0 * ( low + upp - 2.0 * psi ) * ( low - psi) * ( upp - psi ) / ( (upp-low)*(upp-low)*(upp-low)*(upp-low) );}
double dN_D_nhqr(double psi,double upp,double low){return (low * low + 5.0 * low * upp + 2.0 * upp * upp - 7.0 * low * psi - 9.0 * upp * psi + 8.0 * psi * psi ) * ( psi - low) / ( (upp-low)*(upp-low)*(upp-low) );}
double dN_E_nhqr(double psi,double upp,double low){return ( 2.0 * low * low + 5.0 * low * upp + upp * upp - 9.0 * low * psi - 7.0 * upp * psi + 8.0 * psi * psi ) * ( upp - psi ) / ( (upp-low)*(upp-low)*(upp-low) );}


//CUBIC
double N_A_nhcb(double psi,double upp,double low){return ((psi-low)*(psi-low)/((upp-low)*(upp-low)))*(3.0-2.0*(psi-low)/(upp-low));}
double N_B_nhcb(double psi,double upp,double low){return ((psi-upp)*(psi-upp)/((upp-low)*(upp-low)))*(3.0-2.0*(upp-psi)/(upp-low));}
double N_C_nhcb(double psi,double upp,double low){return (psi-upp)*((psi-low)*(psi-low)/((upp-low)*(upp-low)));}
double N_D_nhcb(double psi,double upp,double low){return (psi-low)*((psi-upp)*(psi-upp)/((upp-low)*(upp-low)));}
double dN_A_nhcb(double psi,double upp,double low){return 6.0*(psi-low)*(upp-psi)/((upp-low)*(upp-low)*(upp-low));}
double dN_B_nhcb(double psi,double upp,double low){return 6.0*(low-psi)*(upp-psi)/((upp-low)*(upp-low)*(upp-low));}
double dN_C_nhcb(double psi,double upp,double low){return (psi-low)*(3.0*psi-2.0*upp-low)/((upp-low)*(upp-low));}
double dN_D_nhcb(double psi,double upp,double low){return (upp-psi)*(2.0*low+upp-3.0*psi)/((upp-low)*(upp-low));}

//QUADRATIC

double N_A_nhqd(double psi,double upp,double low){return (2.0*psi-upp-low)*(psi-low)/((upp-low)*(upp-low));}
double N_B_nhqd(double psi,double upp,double low){return (2.0*psi-upp-low)*(psi-upp)/((upp-low)*(upp-low));}
double N_C_nhqd(double psi,double upp,double low){return 4.0*(psi-low)*(upp-psi)/((upp-low)*(upp-low));}
double dN_A_nhqd(double psi,double upp,double low){return (4.0*psi-upp-3.0*low)/((upp-low)*(upp-low));}
double dN_B_nhqd(double psi,double upp,double low){return (4.0*psi-3.0*upp-low)/((upp-low)*(upp-low));}
double dN_C_nhqd(double psi,double upp,double low){return 4.0*(upp+low-2.0*psi)/((upp-low)*(upp-low));}

double N_A_hqd(double psi,double upp,double low){return (3.0*psi-upp-2.0*low) / ( 3.0 * (upp-low) ) ;}
double N_B_hqd(double psi,double upp,double low){return (2.0*upp+low-3.0*psi) / ( 3.0 * (upp-low) ) ;}
double N_C_hqd(double psi,double upp,double low){return 2.0 / 3.0 ;}
double dN_A_hqd(double psi,double upp,double low){return (4.0*psi-upp-3.0*low)/((upp-low)*(upp-low));}
double dN_B_hqd(double psi,double upp,double low){return (4.0*psi-3.0*upp-low)/((upp-low)*(upp-low));}
double dN_C_hqd(double psi,double upp,double low){return 4.0*(upp+low-2.0*psi)/((upp-low)*(upp-low));}

/****************************************************************************************/
//LINEAR

double N_A_nhln(double psi,double upp,double low){return (psi-low)/(upp-low);}
double N_B_nhln(double psi,double upp,double low){return (upp-psi)/(upp-low);}
double dN_A_nhln(double psi,double upp,double low){return 1.0/(upp-low);}
double dN_B_nhln(double psi,double upp,double low){return -1.0/(upp-low);}

double N_A_hln(double psi,double upp,double low){return 0.5;}
double N_B_hln(double psi,double upp,double low){return 0.5;}
double dN_A_hln(double psi,double upp,double low){return 1.0/(upp-low);}
double dN_B_hln(double psi,double upp,double low){return -1.0/(upp-low);}

/****************************************************************************************/
//CONSTANT

double N_A_nhcn(double psi,double upp,double low){return 1.0;}
double N_NULL(double psi,double upp,double low){return 0.0;}



///MOD functions for gq_mod integration

//NHQR
void N_A_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { low * low * (upp + low) * (low - 5.0*upp) * rec , (22.0*low*low*upp + 10.0*upp*upp*low) * rec , - (11.0*low*low + 32.0*upp*low + 5.0*upp*upp) * rec, (18.0*low + 14.0*upp) * rec , -8.0 * rec } ;
  *result = temp ;
}

void N_B_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { upp * upp * (low + upp) * (upp - 5.0*low) * rec , (22.0*upp*upp*low + 10.0*low*low*upp) * rec , - (11.0*upp*upp + 32.0*low*upp + 5.0*low*low) * rec, (18.0*upp + 14.0*low) * rec , -8.0 * rec } ;
  *result = temp ;
}

void N_C_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 16.0 / ( (upp-low) * (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { upp * upp * low * low * rec , -2.0*upp*low*(low + upp) * rec , (upp*upp + 4.0*low*upp + low*low) * rec, -2.0*(upp + low) * rec , rec } ;
  *result = temp ;
}

void N_D_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = {  low*low*upp*(low+upp) * rec , -low*(low*low + 5.0*upp*low + 2.0*upp*upp) * rec , (upp*upp + 7.0*low*upp + 4.0*low*low) * rec, -(3.0*upp + 5.0*low) * rec , 2.0 * rec } ;
  *result = temp ;
}

void N_E_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = {  - upp*upp*low*(upp+low) * rec , upp*(upp*upp + 5.0*low*upp + 2.0*low*low) * rec , - (low*low + 7.0*upp*low + 4.0*upp*upp) * rec, (3.0*low + 5.0*upp) * rec , -2.0 * rec } ;
  *result = temp ;
}

void dN_A_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 2.0 / ( (upp-low) * (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { upp * low * (11.0*low + 5.0*upp) * rec , -(5.0*upp*upp + 32.0*upp*low + 11.0*low*low) * rec , (21.0*upp + 27.0*low) * rec, -16.0 * rec } ;
  *result = temp ;
}

void dN_B_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 2.0 / ( (upp-low) * (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { low * upp * (11.0*upp + 5.0*low) * rec , -(5.0*low*low + 32.0*low*upp + 11.0*upp*upp) * rec , (21.0*low + 27.0*upp) * rec, -16.0 * rec } ;
  *result = temp ;
}

void dN_C_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 32.0 / ( (upp-low) * (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { - low * upp * (upp + low) * rec , (low*low + 4.0*low*upp + upp*upp) * rec , -3.0*(low + upp) * rec, 2.0 * rec } ;
  *result = temp ;
}

void dN_D_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = {  -low*(low*low + 5.0*low*upp + 2.0*upp*upp) * rec , 2.0*(upp*upp + 7.0*low*upp + 4.0*low*low) * rec, -3.0*(3.0*upp + 5.0*low) * rec , 8.0 * rec } ;
  *result = temp ;
}

void dN_E_nhqr_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = {  upp*(upp*upp + 5.0*upp*low + 2.0*low*low) * rec , -2.0*(low*low + 7.0*upp*low + 4.0*upp*upp) * rec, 3.0*(3.0*low + 5.0*upp) * rec , -8.0 * rec } ;
  *result = temp ;
}

//NHCB
void N_A_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { low * low * (3.0*upp - low ) * rec , - 6.0 * upp * low * rec , 3.0 * (upp + low) * rec, -2.0 * rec } ;
  *result = temp ;
}

void N_B_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { upp * upp * (upp - 3.0*low ) * rec , 6.0 * upp * low * rec , -3.0 * (upp + low) * rec, 2.0 * rec } ;
  *result = temp ;
}

void N_C_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { - low * low * upp * rec , ( low*low + 2.0*low*upp ) * rec , - ( upp + 2.0*low ) * rec , rec } ;
  *result = temp ;
}

void N_D_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { - low * upp * upp * rec , ( upp*upp + 2.0*low*upp ) * rec , - ( 2.0*upp + low ) * rec , rec } ;
  *result = temp ;
}

void dN_A_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 6.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { - low * upp * rec , ( upp + low ) * rec , -rec } ;
  *result = temp ;
}

void dN_B_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 6.0 / ( (upp-low) * (upp-low) * (upp-low) ) ;
  std::vector<double> temp = { low * upp * rec , - ( upp + low ) * rec , rec } ;
  *result = temp ;
}

void dN_C_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { low * (2.0*upp + low) * rec , -2.0 * ( upp + 2.0*low ) * rec , 3.0 * rec } ;
  *result = temp ;
}

void dN_D_nhcb_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { upp * (upp + 2.0*low) * rec , -2.0 * ( 2.0*upp + low ) * rec , 3.0 * rec } ;
  *result = temp ;
}

//NHQD
void N_A_nhqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { low * (upp + low) * rec , ( - upp - 3.0 * low ) * rec , 2.0 * rec } ;
  *result = temp ;
}

void N_B_nhqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { upp * (upp + low) * rec , ( - 3.0 * upp - low ) * rec , 2.0 * rec } ;
  *result = temp ;
}

void N_C_nhqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 4.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { - upp * low * rec , (upp + low) * rec , -rec } ;
  *result = temp ;
}

void dN_A_nhqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { ( - upp - 3.0 * low ) * rec , 4.0 * rec } ;
  *result = temp ;
}

void dN_B_nhqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { ( - 3.0 * upp - low ) * rec , 4.0 * rec } ;
  *result = temp ;
}

void dN_C_nhqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { 4.0 * (upp + low) * rec , -8.0 * rec } ;
  *result = temp ;
}

//HQD
void N_A_hqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( 3.0 * (upp-low)) ;
  std::vector<double> temp = { - ( upp + 2.0 * low ) * rec , 3.0 * rec } ;
  *result = temp ;
}

void N_B_hqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( 3.0 * (upp-low)) ;
  std::vector<double> temp = { ( 2.0 * upp + low ) * rec , -3.0 * rec } ;
  *result = temp ;
}

void N_C_hqd_mod(std::vector<double>* result,double upp,double low)
{
  std::vector<double> temp = { 2.0 / 3.0 } ;
  *result = temp ;
}

void dN_A_hqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { ( - upp - 3.0 * low ) * rec , 4.0 * rec } ;
  *result = temp ;
}

void dN_B_hqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { ( - 3.0 * upp - low ) * rec , 4.0 * rec } ;
  *result = temp ;
}

void dN_C_hqd_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / ( (upp-low) * (upp-low)) ;
  std::vector<double> temp = { 4.0 * (upp + low) * rec , -8.0 * rec } ;
  *result = temp ;
}


//NHLN
void N_A_nhln_mod(std::vector<double>* result,double upp,double low)
{
double rec = 1.0 / (upp-low) ;
 std::vector<double> temp = { -low * rec , rec } ;
  *result = temp ;
}
void N_B_nhln_mod(std::vector<double>* result,double upp,double low)
{
double rec = 1.0 / (upp-low) ;
 std::vector<double> temp = { upp * rec , -rec } ;
  *result = temp ;
}
void dN_A_nhln_mod(std::vector<double>* result,double upp,double low)
{
double rec = 1.0 / (upp-low) ;
  std::vector<double> temp = { rec } ;
  *result = temp ;
}
void dN_B_nhln_mod(std::vector<double>* result,double upp,double low)
{
double rec = 1.0 / (upp-low) ;
  std::vector<double> temp = { -rec } ;
  *result = temp ;
}

void N_A_hln_mod(std::vector<double>* result,double upp,double low)
{
  std::vector<double> temp = {0.5} ;
  *result = temp ;
}
void N_B_hln_mod(std::vector<double>* result,double upp,double low)
{
  std::vector<double> temp = {0.5} ;
  *result = temp ;
}
void dN_A_hln_mod(std::vector<double>* result,double upp,double low)
{
 double rec = 1.0 / (upp-low) ;
  std::vector<double> temp = { rec } ;
  *result = temp ;
}
void dN_B_hln_mod(std::vector<double>* result,double upp,double low)
{
  double rec = 1.0 / (upp-low) ;
  std::vector<double> temp = { -rec } ;
  *result = temp ;
}

//NHCN
void N_A_nhcn_mod(std::vector<double>* result,double upp,double low)
{
  std::vector<double> temp = {1.0} ;
  *result = temp ;
}
void N_NULL_mod(std::vector<double>* result,double upp,double low)
{
  std::vector<double> temp = {0.0} ;
  *result = temp ;
}
