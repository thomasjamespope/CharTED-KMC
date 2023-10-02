# Makefile for KMC. T.Pope, 2020
FC          = mpifort 
FLAGS_DEBUG = -g -fcheck=all -Wall
FLAGS_OPT   = -O3
FLAGS_MPI   =
FLAGS_FORM  = -ffixed-line-length-none -ffree-form 

FLAGS       =  $(FLAGS_MPI) $(FLAGS_OPT) $(FLAGS_FORM) $(FLAGS_DEBUG)

.f.o:
	$(FC) -c $(FLAGS) $< 

OBJS = kinds.o \
       param.o \
       general_functions.o \
       type_indexes.o \
       type_coor.o \
       type_charge.o \
       mpi.o \
       readings.o \
       printings.o \
       checks.o \
       initialization.o  \
       timings.o \
       pbc.o \
       coulomb.o \
       random.o \
       event.o

ai: $(OBJS)
	

