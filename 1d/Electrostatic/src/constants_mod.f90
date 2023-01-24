MODULE constants_mod

IMPLICIT NONE


! simulation and grid parameters
INTEGER, PARAMETER :: nx       = 32    			   ! mesh size (has to be a power of 2)
INTEGER, PARAMETER :: ns       = 1     			   ! number of particle species
INTEGER, PARAMETER :: N        = 1024  		       ! max number of particles per species (needs to be a multiple of nx)
INTEGER, PARAMETER :: maxsteps = 300              ! max. number of time steps
REAL*8,  PARAMETER :: L        = 8.d0*ATAN(1.d0)            ! grid length
REAL*8,  PARAMETER :: twopi    = 8.d0*ATAN(1.d0)
INTEGER, PARAMETER :: interpolation_type = 2       ! 1: zeroth order (NGP), 2: first order (CIC), 3: first order for particles and zeroth order for fields
REAL*8, PARAMETER :: a1 = 0.d0, a2 = 1000.d0 
LOGICAL, PARAMETER :: print_debug = .FALSE.
INTEGER, PARAMETER :: tskip = 1
INTEGER, PARAMETER :: kmax = nx-1 

! init parameters
INTEGER, PARAMETER :: vbins = 500   ! velocity space bins 
REAL*8, PARAMETER :: v0   = 0.0d0  ! intial drift velocity
REAL*8, PARAMETER :: vt   = 0.5d0  ! R.M.S. thermal velocity
REAL*8, PARAMETER :: vmax = 3.0*vt ! Maxwellian velocity distribution cut-off
REAL*8 :: rho0 = 0.d0   ! neutralizing background charge distribution
LOGICAL, PARAMETER :: load_maxwellian = .TRUE. 
INTEGER, PARAMETER :: nlg = 1   ! number of velocity loading groups

! physics parameters
REAL*8, PARAMETER :: eps0 = 1.d0       ! vacuum permittivity
REAL*8, PARAMETER :: qm_e = -1.d0      ! charge to mass ratio - electrons
REAL*8, PARAMETER :: qm_i = 0.0001d0      ! charge to mass ratio - ions
REAL*8, PARAMETER :: wp_e = 1.d0       ! plasma frequency - electrons
REAL*8, PARAMETER :: wp_i = 0.01d0       ! plasma frequency -  ions
REAL*8, PARAMETER :: wc_e = 0.d0  	   ! cyclotron frequency - electrons
REAL*8, PARAMETER :: wc_i = 0.d0       ! cyclotron frequency - ions


! File output parameters
LOGICAL, PARAMETER :: fields_out = .TRUE. 
LOGICAL, PARAMETER :: fieldsk_out = .TRUE. 
LOGICAL, PARAMETER :: particles_out = .TRUE. 
LOGICAL, PARAMETER :: modes_out = .TRUE.
LOGICAL, PARAMETER :: energy_out = .TRUE.
LOGICAL, PARAMETER :: fv_out = .TRUE.

END MODULE constants_mod