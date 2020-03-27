#ifndef FILL_EQUILIBRIUM_H
#define FILL_EQUILIBRIUM_H

#include "constants.h"

void fill_theta(equil_fields *equil,geom_shape geom);
void fill_screw(equil_fields *equil,geom_shape geom);
void fill_full(equil_fields *equil,geom_shape geom);
int fill_whales(equil_fields *equil,geom_shape *geom);

#endif
