MODULE constants_mod

IMPLICIT NONE


! simulation and grid parameters
INTEGER, PARAMETER :: nx       = 300    	     ! mesh size 
INTEGER, PARAMETER :: ny       = 150
INTEGER, PARAMETER :: nz       = 150

INTEGER, PARAMETER :: ns       = 2     			 ! number of particle species
INTEGER, PARAMETER :: Np_max   = 1500000    	 ! max number of particles per species

INTEGER, PARAMETER :: maxsteps = 800             ! max. number of time steps
REAL*8, PARAMETER  :: COUR     = 0.5             ! time step courant number
INTEGER, PARAMETER :: tskip    = 1               ! file output frequency


INTEGER, PARAMETER :: interpolation_type = 2     ! 1: zeroth order (NGP), 2: first order (CIC), 3: first order for particles and zeroth order for fields
INTEGER, PARAMETER :: bndry = 2                  ! 1: periodic (both x,y) , 2: open (i.e. radiation absorbing and particle outflow), 3:reflecting 
INTEGER, PARAMETER :: current_type = 1           ! 1: Zig-zag scheme 2: Esirkepov's scheme
INTEGER, PARAMETER :: N_filter = 1               ! how many times to repeat current filter (1 seems to be good enough for coarser grids)

LOGICAL, PARAMETER :: current_filter_on = .TRUE.
LOGICAL, PARAMETER :: print_debug       = .FALSE.
LOGICAL, PARAMETER :: override_checks   = .TRUE.
LOGICAL, PARAMETER :: deposit_charges_On  = .FALSE.

REAL*8,  PARAMETER :: twopi  = 8.d0*ATAN(1.d0)
REAL*8 , PARAMETER :: third  = 1.d0 / 3.d0 

! physics parameters
REAL*8, PARAMETER :: eps0 = 1.d0        ! vacuum permittivity
REAL*8, PARAMETER :: c    = 0.5d0       ! speed of light in vacuum 
REAL*8, PARAMETER :: ee   = 0.00000001!1.d0          ! electron charge
REAL*8, PARAMETER :: qm_e = -1.d0       ! charge to mass ratio - electrons
REAL*8, PARAMETER :: qm_i = 1.d0/16.d0  ! charge to mass ratio - ions
REAL*8, PARAMETER :: wp_e = 0.0001d0    ! plasma frequency - electrons
REAL*8, PARAMETER :: wp_i = 0.0001d0    ! plasma frequency -  ions
REAL*8, PARAMETER :: wc_e = 0.0d0  	    ! cyclotron frequency - electrons
REAL*8, PARAMETER :: wc_i = 0.d0        ! cyclotron frequency - ions
 
! boundary condition parameters 
REAL*8, PARAMETER :: cE = 0.5858
REAL*8, PARAMETER :: cB = 1- cE

! init parameters
REAL*8, PARAMETER  :: vmin  = -3.d0      ! lower cut-off for velocity space distribution
REAL*8, PARAMETER  :: vmax  = 3.d0       ! upper cut-off for velocity space distribution
INTEGER, PARAMETER :: vbins = 1000       ! velocity space bins 
REAL*8, PARAMETER  :: v0x   = 0.25d0     ! intial drift velocity
REAL*8, PARAMETER  :: v0y   = 0.d0 
REAL*8, PARAMETER  :: v0z   = 0.0d0 
REAL*8, PARAMETER  :: vth_e    = 0.1d0       ! R.M.S. thermal velocity
REAL*8, PARAMETER  :: vth_i    = 0.025d0     ! R.M.S. thermal velocity
REAL*8, PARAMETER  :: E0x   = 0.0d0      ! external electric field
REAL*8, PARAMETER  :: E0y   = 0.0d0
REAL*8, PARAMETER  :: E0z   = 0.0d0
REAL*8 :: rho0 = 0.d0                    ! neutralizing background charge density

LOGICAL, PARAMETER :: randPos = .FALSE.  ! setting for radomizing particle positions

! File output parameters
LOGICAL, PARAMETER :: fields_out    = .TRUE. 
LOGICAL, PARAMETER :: particles_out = .TRUE. 
LOGICAL, PARAMETER :: energy_out    = .FALSE.
LOGICAL, PARAMETER :: fv_out        = .FALSE.
LOGICAL, PARAMETER :: current_out   = .FALSE.

END MODULE constants_mod