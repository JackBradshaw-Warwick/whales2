#ifndef BUILD_EQUILIBRIUM_H
#define BUILD_EQUILIBRIUM_H

#include "constants.h"
#include "whales_structs.h"

void read_geom(geom_shape *geom);
void calc_geom(geom_shape *geom);
void build_geom(geom_shape *geom);
void build_equil(equil_fields *equil,geom_shape *geom);
void build_jacobian(equil_fields *equil,geom_shape geom);
void build_equil_funcs(equil_fields *equil,geom_shape geom);
void rad_to_psi(equil_fields *equil,int N_psi,int N_interp);
void fill_cc(equil_fields *equil,int N_interp,int N_theta);
void fill_rad(equil_fields *equil,geom_shape geom, double flux_max, double flux_min);

#endif
