# Compiler
FC = gfortran

# The CFLAGS variable sets compile flags for gfortran:
#  -g                            compile with debug information
#  -Wall                         give verbose compiler warnings
#  -O0                           do not optimize generated code
#  -O2 -ftree-vectorize          enable auto-vectorization only
#  -O3                           apply level 3 optimizations
#  -std=gnu99                    use the GNU99 standard language definition
#  -fopt-info-all=opt_report.all save optimization report to file
# -fopenmp                       enable openmp
# -fcheck=all                    enable array bounds checking

CFLAGS = -O2 -ftree-vectorize -fopenmp #-fcheck=all 

OBJS = constants_mod.o data_mod.o init_mod.o \
       fieldSolver_mod.o particleMover_mod.o driver.o 

em_3d : $(OBJS)
	$(FC) $(CFLAGS) -o em_3d $(OBJS)
	
constants_mod.o : src/constants_mod.f90
	$(FC) $(CFLAGS) -c src/constants_mod.f90
		
data_mod.o : src/data_mod.f90 
	$(FC) $(CFLAGS) -c src/data_mod.f90	
	
fieldSolver_mod.o : src/fieldSolver_mod.f90  constants_mod.o data_mod.o  
	$(FC) $(CFLAGS) -c src/fieldSolver_mod.f90

particleMover_mod.o : src/particleMover_mod.f90  constants_mod.o data_mod.o  
	$(FC) $(CFLAGS) -c src/particleMover_mod.f90
	
init_mod.o : src/init_mod.f90 constants_mod.o data_mod.o particleMover_mod.o
	$(FC) $(CFLAGS) -c src/init_mod.f90
	
driver.o : src/driver.f90  constants_mod.o data_mod.o init_mod.o fieldSolver_mod.o particleMover_mod.o 
	$(FC) $(CFLAGS) -c src/driver.f90
	
clean:
	rm  -f *.o *.mod *.txt