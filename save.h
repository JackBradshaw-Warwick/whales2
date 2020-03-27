#ifndef SAVE_H
#define SAVE_H

#include "constants.h"

int save_equil_full(geom_shape geom,equil_fields equil);
int save_equil_part(geom_shape geom,equil_fields equil);
int load_mats(geom_shape geom , Mat mats[] );
int save_mats(geom_shape geom , Mat mats[] );

#endif
