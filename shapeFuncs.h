#ifndef SPACFUNCS_H
#define SPACFUNCS_H

#include "constants.h"

double N_A_nhqr(double psi,double psi_upp,double psi_low);
double N_B_nhqr(double psi,double psi_upp,double psi_low);
double N_C_nhqr(double psi,double psi_upp,double psi_low);
double N_D_nhqr(double psi,double psi_upp,double psi_low);
double N_E_nhqr(double psi,double psi_upp,double psi_low);
double dN_A_nhqr(double psi,double psi_upp,double psi_low);
double dN_B_nhqr(double psi,double psi_upp,double psi_low);
double dN_C_nhqr(double psi,double psi_upp,double psi_low);
double dN_D_nhqr(double psi,double psi_upp,double psi_low);
double dN_E_nhqr(double psi,double psi_upp,double psi_low);

double N_A_nhcb(double psi,double psi_upp,double psi_low);
double N_B_nhcb(double psi,double psi_upp,double psi_low);
double N_C_nhcb(double psi,double psi_upp,double psi_low);
double N_D_nhcb(double psi,double psi_upp,double psi_low);
double dN_A_nhcb(double psi,double psi_upp,double psi_low);
double dN_B_nhcb(double psi,double psi_upp,double psi_low);
double dN_C_nhcb(double psi,double psi_upp,double psi_low);
double dN_D_nhcb(double psi,double psi_upp,double psi_low);

double N_A_nhqd(double psi,double psi_upp,double psi_low);
double N_B_nhqd(double psi,double psi_upp,double psi_low);
double N_C_nhqd(double psi,double psi_upp,double psi_low);
double dN_A_nhqd(double psi,double psi_upp,double psi_low);
double dN_B_nhqd(double psi,double psi_upp,double psi_low);
double dN_C_nhqd(double psi,double psi_upp,double psi_low);

double N_A_hqd(double psi,double psi_upp,double psi_low);
double N_B_hqd(double psi,double psi_upp,double psi_low);
double N_C_hqd(double psi,double psi_upp,double psi_low);
double dN_A_hqd(double psi,double psi_upp,double psi_low);
double dN_B_hqd(double psi,double psi_upp,double psi_low);
double dN_C_hqd(double psi,double psi_upp,double psi_low);

double N_A_nhln(double psi,double psi_upp,double psi_low);
double N_B_nhln(double psi,double psi_upp,double psi_low);
double dN_A_nhln(double psi,double psi_upp,double psi_low);
double dN_B_nhln(double psi,double psi_upp,double psi_low);

double N_A_hln(double psi,double psi_upp,double psi_low);
double N_B_hln(double psi,double psi_upp,double psi_low);
double dN_A_hln(double psi,double psi_upp,double psi_low);
double dN_B_hln(double psi,double psi_upp,double psi_low);

double N_A_nhcn(double psi,double psi_upp,double psi_low);
double N_NULL(double psi,double psi_upp,double psi_low);

void N_A_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_B_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_C_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_D_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_E_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_A_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_B_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_C_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_D_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_E_nhqr_mod(std::vector<double>* result,double psi_upp,double psi_low);

void N_A_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_B_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_C_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_D_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_A_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_B_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_C_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_D_nhcb_mod(std::vector<double>* result,double psi_upp,double psi_low);

void N_A_nhqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_B_nhqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_C_nhqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_A_nhqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_B_nhqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_C_nhqd_mod(std::vector<double>* result,double psi_upp,double psi_low);

void N_A_hqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_B_hqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_C_hqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_A_hqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_B_hqd_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_C_hqd_mod(std::vector<double>* result,double psi_upp,double psi_low);

void N_A_nhln_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_B_nhln_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_A_nhln_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_B_nhln_mod(std::vector<double>* result,double psi_upp,double psi_low);

void N_A_hln_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_B_hln_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_A_hln_mod(std::vector<double>* result,double psi_upp,double psi_low);
void dN_B_hln_mod(std::vector<double>* result,double psi_upp,double psi_low);

void N_A_nhcn_mod(std::vector<double>* result,double psi_upp,double psi_low);
void N_NULL_mod(std::vector<double>* result,double psi_upp,double psi_low);

#endif
