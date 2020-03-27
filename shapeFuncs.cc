/****************************************************************************************/
#include "constants.h"


//Implement shape functions in terms of psi,psi_upp,psi_low. Saves doing the whole coordinate transform


//QUARTIC
double N_A_nhqr(double psi,double psi_upp,double psi_low){return (psi - psi_low) * (psi - psi_low) * (psi_low * psi_low - 4.0 * psi_low * psi_upp - 5.0 * psi_upp * psi_upp + 2.0 * psi_low * psi + 14.0 * psi_upp * psi - 8.0 * psi * psi ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double N_B_nhqr(double psi,double psi_upp,double psi_low){return (psi - psi_upp) * (psi - psi_upp) * ( - 5.0 * psi_low * psi_low - 4.0 * psi_low * psi_upp + psi_upp * psi_upp + 14.0 * psi_low * psi + 2.0 * psi_upp * psi - 8.0 * psi * psi ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double N_C_nhqr(double psi,double psi_upp,double psi_low){return 16.0 * ( psi - psi_low ) * ( psi - psi_low ) * ( psi - psi_upp ) * ( psi - psi_upp ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double N_D_nhqr(double psi,double psi_upp,double psi_low){return ( psi_low + psi_upp - 2.0 * psi ) * ( psi_upp - psi ) * ( psi_low - psi) * ( psi_low - psi) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double N_E_nhqr(double psi,double psi_upp,double psi_low){return - ( psi_low + psi_upp - 2.0 * psi ) * ( psi_upp - psi ) * ( psi_upp - psi) * ( psi_low - psi) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double dN_A_nhqr(double psi,double psi_upp,double psi_low){return 2.0 * ( 11.0 * psi_low + 5.0 * psi_upp - 16.0 * psi ) * ( psi_low - psi) * ( psi_upp - psi ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double dN_B_nhqr(double psi,double psi_upp,double psi_low){return 2.0 * ( 5.0 * psi_low + 11.0 * psi_upp - 16.0 * psi ) * ( psi_low - psi) * ( psi_upp - psi ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double dN_C_nhqr(double psi,double psi_upp,double psi_low){return - 32.0 * ( psi_low + psi_upp - 2.0 * psi ) * ( psi_low - psi) * ( psi_upp - psi ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double dN_D_nhqr(double psi,double psi_upp,double psi_low){return (psi_low * psi_low + 5.0 * psi_low * psi_upp + 2.0 * psi_upp * psi_upp - 7.0 * psi_low * psi - 9.0 * psi_upp * psi + 8.0 * psi * psi ) * ( psi - psi_low) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}
double dN_E_nhqr(double psi,double psi_upp,double psi_low){return ( 2.0 * psi_low * psi_low + 5.0 * psi_low * psi_upp + psi_upp * psi_upp - 9.0 * psi_low * psi - 7.0 * psi_upp * psi + 8.0 * psi * psi ) * ( psi_upp - psi ) / ( (psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low) );}


//CUBIC
double N_A_nhcb(double psi,double psi_upp,double psi_low){return ((psi-psi_low)*(psi-psi_low)/((psi_upp-psi_low)*(psi_upp-psi_low)))*(3.0-2.0*(psi-psi_low)/(psi_upp-psi_low));}
double N_B_nhcb(double psi,double psi_upp,double psi_low){return ((psi-psi_upp)*(psi-psi_upp)/((psi_upp-psi_low)*(psi_upp-psi_low)))*(3.0-2.0*(psi_upp-psi)/(psi_upp-psi_low));}
double N_C_nhcb(double psi,double psi_upp,double psi_low){return (psi-psi_upp)*((psi-psi_low)*(psi-psi_low)/((psi_upp-psi_low)*(psi_upp-psi_low)));}
double N_D_nhcb(double psi,double psi_upp,double psi_low){return (psi-psi_low)*((psi-psi_upp)*(psi-psi_upp)/((psi_upp-psi_low)*(psi_upp-psi_low)));}
double dN_A_nhcb(double psi,double psi_upp,double psi_low){return 6.0*(psi-psi_low)*(psi_upp-psi)/((psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low));}
double dN_B_nhcb(double psi,double psi_upp,double psi_low){return 6.0*(psi_low-psi)*(psi_upp-psi)/((psi_upp-psi_low)*(psi_upp-psi_low)*(psi_upp-psi_low));}
double dN_C_nhcb(double psi,double psi_upp,double psi_low){return (psi-psi_low)*(3.0*psi-2.0*psi_upp-psi_low)/((psi_upp-psi_low)*(psi_upp-psi_low));}
double dN_D_nhcb(double psi,double psi_upp,double psi_low){return (psi_upp-psi)*(2.0*psi_low+psi_upp-3.0*psi)/((psi_upp-psi_low)*(psi_upp-psi_low));}

//QUADRATIC

double N_A_nhqd(double psi,double psi_upp,double psi_low){return (2.0*psi-psi_upp-psi_low)*(psi-psi_low)/((psi_upp-psi_low)*(psi_upp-psi_low));}
double N_B_nhqd(double psi,double psi_upp,double psi_low){return (2.0*psi-psi_upp-psi_low)*(psi-psi_upp)/((psi_upp-psi_low)*(psi_upp-psi_low));}
double N_C_nhqd(double psi,double psi_upp,double psi_low){return 4.0*(psi-psi_low)*(psi_upp-psi)/((psi_upp-psi_low)*(psi_upp-psi_low));}
double dN_A_nhqd(double psi,double psi_upp,double psi_low){return (4.0*psi-psi_upp-3.0*psi_low)/((psi_upp-psi_low)*(psi_upp-psi_low));}
double dN_B_nhqd(double psi,double psi_upp,double psi_low){return (4.0*psi-3.0*psi_upp-psi_low)/((psi_upp-psi_low)*(psi_upp-psi_low));}
double dN_C_nhqd(double psi,double psi_upp,double psi_low){return 4.0*(psi_upp+psi_low-2.0*psi)/((psi_upp-psi_low)*(psi_upp-psi_low));}

double N_A_hqd(double psi,double psi_upp,double psi_low){return 0.0;}
double N_B_hqd(double psi,double psi_upp,double psi_low){return 0.0;}
double N_C_hqd(double psi,double psi_upp,double psi_low){return 0.0;}
double dN_A_hqd(double psi,double psi_upp,double psi_low){return 0.0;}
double dN_B_hqd(double psi,double psi_upp,double psi_low){return 0.0;}
double dN_C_hqd(double psi,double psi_upp,double psi_low){return 0.0;}

/****************************************************************************************/
//LINEAR

double N_A_nhln(double psi,double psi_upp,double psi_low){return (psi-psi_low)/(psi_upp-psi_low);}
double N_B_nhln(double psi,double psi_upp,double psi_low){return (psi_upp-psi)/(psi_upp-psi_low);}
double dN_A_nhln(double psi,double psi_upp,double psi_low){return 1.0/(psi_upp-psi_low);}
double dN_B_nhln(double psi,double psi_upp,double psi_low){return -1.0/(psi_upp-psi_low);}

double N_A_hln(double psi,double psi_upp,double psi_low){return 0.5;}
double N_B_hln(double psi,double psi_upp,double psi_low){return 0.5;}
double dN_A_hln(double psi,double psi_upp,double psi_low){return 1.0/(psi_upp-psi_low);}
double dN_B_hln(double psi,double psi_upp,double psi_low){return -1.0/(psi_upp-psi_low);}

/****************************************************************************************/
//CONSTANT

double N_A_nhcn(double psi,double psi_upp,double psi_low){return 1.0;}
double N_NULL(double psi,double psi_upp,double psi_low){return 0.0;}
