FM = fermi_bath.o sovler_fermi.o solver_multi_op.o
DENSITY =  density_bath.o  
DIAG = diag_flip_bath.o solver_diag_flip.o 
DIAGFM = solver_diag_flip_fermi_bath.o
OFFDIAG = off_diag_flip_bath.o solver_off_diag_flip.o
OFFDIAGFM = solver_off_diag_flip_fermi_bath.o
TRACE = local_config.o trace.o  
MEAS = accumulator.o 
AUX = fourier_transformation.o
MAIN = solver.o main.o
OBJS = $(FM) $(DENSITY) $(DIAG) $(DIAGFM) $(OFFDIAG) $(OFFDIAGFM) $(TRACE) $(MEAS) $(AUX) $(MAIN)
CC = mpic++ -std=c++11
#DEBUG = -traceback -g -debug all
CFLAGS = -c -Wall -O3 $(DEBUG) -I /home/hh25/boost_1_70_0

main : $(OBJS)
	$(CC) -o main $(OBJS)

%.o : %.cpp
	$(CC) $(CFLAGS) $<
clean:
	\rm *.o *~main

