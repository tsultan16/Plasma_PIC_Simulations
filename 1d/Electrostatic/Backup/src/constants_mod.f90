MODULE constants_mod

IMPLICIT NONE


! simulation and grid parameters
INTEGER, PARAMETER :: nx       = 32    			   ! mesh size (has to be a power of 2)
INTEGER, PARAMETER :: ns       = 2     			   ! number of particle species
INTEGER, PARAMETER :: N        = 128  		       ! max number of particles per species (needs to be a multiple of nx)
INTEGER, PARAMETER :: maxsteps = 500              ! max. number of time steps
REAL*8,  PARAMETER :: L        = 100.d0            ! grid length
REAL*8,  PARAMETER :: twopi    = 8.d0*ATAN(1.d0)
INTEGER, PARAMETER :: interpolation_type = 3       ! 1: zeroth order (NGP), 2: first order (CIC), 3: first order for particles and zeroth order for fields
REAL*8, PARAMETER :: a1 = 0.5d0, a2 = 0.5d0 
LOGICAL, PARAMETER :: print_debug = .FALSE.

! init parameters
REAL*8, PARAMETER :: vmin = -1.d0  ! lower cut-off for intial velocity space distribution
REAL*8, PARAMETER :: vmax = 1.d0   ! upper cut-off for intial velocity space distribution 
REAL*8, PARAMETER :: v0   = 0.d0   ! intial drift velocity
REAL*8, PARAMETER :: vt   = 0.d0   ! R.M.S. thermal velocity
REAL*8, PARAMETER :: rho0 = 0.d0   ! neutralizing background charge distribution

! physics parameters
REAL*8, PARAMETER :: eps0 = 1.d0       ! vacuum permittivity
REAL*8, PARAMETER :: qm_e = -4.d0  ! charge to mass ratio - electrons
REAL*8, PARAMETER :: qm_i = 1.d0       ! charge to mass ratio - ions
REAL*8, PARAMETER :: wp_e = 2.d0     ! plasma frequency - electrons
REAL*8, PARAMETER :: wp_i = 1.d0       ! plasma frequency -  ions
REAL*8, PARAMETER :: wc_e = 0.d0  	   ! cyclotron frequency - electrons
REAL*8, PARAMETER :: wc_i = 0.d0       ! cyclotron frequency - ions

! File output parameters
LOGICAL, PARAMETER :: fields_out = .TRUE. 
LOGICAL, PARAMETER :: fieldsk_out = .TRUE. 
LOGICAL, PARAMETER :: particles_out = .TRUE. 
LOGICAL, PARAMETER :: modes_out = .TRUE.
LOGICAL, PARAMETER :: energy_out = .TRUE.

END MODULE constants_mod