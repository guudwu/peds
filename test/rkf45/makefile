# Test Runge-Kutta-Fehlberg algorithm

# Last modified on 2013-10-22 by GuuD WU.

# 8x8 with 32 non-zero entries (Non-sparse input).
# Dense intercept.
# Time-Point from 0 to 10 with step 1.
# Observation 8x11.

all: rkf45_test.out

PATH_ROOT=$(PEDS)
PATH_ENV=$(PATH_ROOT)
PATH_RKF45=$(PATH_ROOT)/solver/rkf45
PATH_LINEAR_ODE=$(PATH_ROOT)/model/linear_ode
PATH_TEST=$(PATH_ROOT)/test/rkf45

COMPILER=gfortran -fdefault-real-8
STATIC_LIB_COMPILER=$(COMPILER) -c

rkf45_test.out: ENV.o rkf45.o d1mach.o Linear_Ode_mod.o rkf45_test.o
	$(COMPILER) ENV.o rkf45.o d1mach.o \
	  Linear_Ode_mod.o rkf45_test.o \
	  -o rkf45_test.out

ENV.o: $(PATH_ENV)/ENV.f90
	$(STATIC_LIB_COMPILER) $(PATH_ENV)/ENV.f90 \
	  -o ENV.o
rkf45.o: $(PATH_RKF45)/rkf45.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/rkf45.f \
	  -o rkf45.o
d1mach.o: $(PATH_RKF45)/d1mach.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/d1mach.f \
	  -o d1mach.o
Linear_Ode_mod.o: $(PATH_LINEAR_ODE)/Linear_Ode_mod.f90
	$(STATIC_LIB_COMPILER) $(PATH_LINEAR_ODE)/Linear_Ode_mod.f90 \
	  -o Linear_Ode_mod.o
rkf45_test.o: $(PATH_TEST)/rkf45_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/rkf45_test.f90 \
	  -o rkf45_test.o

clear:
	rm *.o *.mod *.out
