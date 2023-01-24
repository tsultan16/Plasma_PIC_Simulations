**********************3D Electromagnetic PIC Code************************

Some notes

(1) Grid:        dx,dy,dz and dt are set to 1 code units.
                 The speed of light is then set equal to
                 the maximum Courant number for stability
                 (which is 0.5 in 3D). The Bx, By,Bz arrays
                 are actually c*Bx, c*By, c*Bz. This 
                 makes E <-> B symmetry in Maxwell's equations.
                 We also set epsilon_0 = 1 which makes 
                 mu_0 = 1/c^2.


(2) Code stability: Three conditions need to be met to ensure 
                    that the time integration of particles are
                    fields remain stable. These are: 
                    dx <= lambda_debye, 
                    dt < dx/c*sqrt(3) and
                    dt < 2/ w_pe
                    i.e. debeye length needs to be spatially resolved,
                    electromagnetic waves must not propagate more than 
                    1/sqrt(3) cell lengths in a single time-step and
                    time resolutins needs to be high enough to resolve 
                    at least a few plasma oscillation periods.
************************************************************************

Features:

-> Current conservation algorithms: 

* Esirkepov algorithm (2001)
* Umeda Zig-zag algorithm (2003)

-> Boundary Conditions:
   
* Radiation outflow boundary conditions (lowest order Lindman method)
* Periodic boundaries

-> Quasi-particle shape factors:

* First Order (piecewise linear)


TO DO List:

-> Need to add conducting boundary conditions to neutralize outgoing charges


Parallelization:

-> OpenMP Thread Parallelization : Version 3 in-progress

-> MPI Domain Decomposition

