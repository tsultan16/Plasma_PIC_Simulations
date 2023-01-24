OBJS = constants_mod.o data_mod.o init_mod.o \
       fieldSolver_mod.o particleMover_mod.o driver.o 

em_3d : $(OBJS)
	gfortran -fopenmp -fcheck=all -o em_3d $(OBJS)
	
constants_mod.o : src/constants_mod.f90
	gfortran -fopenmp -c src/constants_mod.f90
		
data_mod.o : src/data_mod.f90 
	gfortran -fopenmp -c src/data_mod.f90	
	
fieldSolver_mod.o : src/fieldSolver_mod.f90  constants_mod.o data_mod.o  
	gfortran -fopenmp -c src/fieldSolver_mod.f90

particleMover_mod.o : src/particleMover_mod.f90  constants_mod.o data_mod.o  
	gfortran -fopenmp -c src/particleMover_mod.f90
	
init_mod.o : src/init_mod.f90 constants_mod.o data_mod.o particleMover_mod.o
	gfortran -fopenmp -c src/init_mod.f90
	
driver.o : src/driver.f90  constants_mod.o data_mod.o init_mod.o fieldSolver_mod.o particleMover_mod.o 
	gfortran -fopenmp -c src/driver.f90
	
clean:
	rm  -f *.o *.mod *.txt