# Makefile for KMC. T.Pope, 2020
FC          = mpifort 
FLAGS_DEBUG = -g -fcheck=all -Wall
FLAGS_OPT   = -O3
FLAGS_MPI   = 
FLAGS_FORM  = -ffixed-line-length-none -ffree-form
LIBRARY     = 
MODULES     = ./modules
FLAGS       = $(FLAGS_MPI) $(FLAGS_DEBUG) $(FLAGS_OPT) $(FLAGS_FORM) 
PROGRAM     = CharTED-KMC.x
PRECOMP     = CharTED-KMC.o
SOURCE      = CharTED-KMC.f90

MODOBJS= $(MODULES)/kinds.o \
         $(MODULES)/param.o \
         $(MODULES)/general_functions.o \
         $(MODULES)/type_indexes.o \
         $(MODULES)/type_coor.o \
         $(MODULES)/type_charge.o \
         $(MODULES)/mpi.o \
         $(MODULES)/readings.o \
         $(MODULES)/printings.o \
         $(MODULES)/checks.o \
         $(MODULES)/initialization.o  \
         $(MODULES)/timings.o \
         $(MODULES)/pbc.o \
         $(MODULES)/coulomb.o \
         $(MODULES)/random.o \
         $(MODULES)/event.o

.PHONY: all $(MODULES) $(PROGRAM)

$(MODULES):
	$(MAKE) --directory=$@

all: $(MODULES)
	$(FC) -I$(MODULES) $(FLAGS) $(LIBRARY) -c $(SOURCE) -o $(PRECOMP)
	$(FC) $(MODOBJS) $(PRECOMP) $(FLAGS) $(LIBRARY) -o $(PROGRAM)
clean:
	rm $(MODULES)/*.o $(MODULES)/*.mod $(PRECOMP)
