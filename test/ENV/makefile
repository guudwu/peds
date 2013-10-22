# Test module ENV

# Last modified on 2013-10-22 by GuuD WU.

all: ENV_test.out

PATH_ROOT=$(PEDS)
PATH_ENV=$(PATH_ROOT)
PATH_TEST=$(PATH_ROOT)/test/ENV

COMPILER=gfortran -fdefault-real-8
STATIC_LIB_COMPILER=$(COMPILER) -c

ENV_test.out: ENV.o ENV_test.o
	$(COMPILER) ENV.o ENV_test.o \
	  -o ENV_test.out

ENV.o: $(PATH_ENV)/ENV.f90
	$(STATIC_LIB_COMPILER) $(PATH_ENV)/ENV.f90 \
	  -o ENV.o

ENV_test.o: $(PATH_TEST)/ENV_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/ENV_test.f90 \
	  -o ENV_test.o

clear:
	rm *.o *.mod
