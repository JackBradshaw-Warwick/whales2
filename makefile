ALL:
CC=gcc
CFLAGS=-I. -I${BOOST_DIR} ${SLEPC_CC_INCLUDES}

DEPS = build_equilibrium.h build_matrices.h constants.h derivatives.h equilibrium_funcs.h fill_equilibrium.h interp.h process_results.h save.h shapeFuncs.h soloviev.h utilities.h 
OBJ = build_equilibrium.o build_matrices.o derivatives.o equilibrium_funcs.o fill_equilibrium.o interp.o process_results.o save.o shapeFuncs.o soloviev_equil.o soloviev_simple.o utilities.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

#--------------------------------------------------------------------------

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

whales2: whales2.o $(OBJ) chkopts
	-${CLINKER} $(OBJ) -o ../whales2 whales2.o ${SLEPC_LIB} -L/usr/local/lib -L${BOOST_DIR}/stage/lib -lfftw3 -lm
	${RM} whales2.o $(OBJ)


