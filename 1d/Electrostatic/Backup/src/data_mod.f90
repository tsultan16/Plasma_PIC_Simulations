 MODULE data_mod
 
 IMPLICIT NONE	
 
 TYPE particle
    ! type data
    REAL*8 :: q = 0.d0   ! charge
	REAL*8 :: m = 0.d0   ! mass
	REAL*8 :: x = 0.d0   ! position
	REAL*8 :: vx = 0.d0  ! velocity-x
	REAL*8 :: vy = 0.d0  ! velocity-y
	
 END TYPE particle
 
 
 REAL*8, ALLOCATABLE :: rho(:), phi(:), E(:), ne(:), ni(:)
 REAL*8, ALLOCATABLE :: rhok(:,:), phik(:,:), Ksqr(:,:), ESEk(:), SMk(:), ESEkt(:,:)
 TYPE(particle), ALLOCATABLE :: particles(:,:)
 REAL*8 :: n0  
 REAL*8 :: dx, dt
 REAL*8 :: KE, Px, Py, ESE, ESE2
 REAL*8 :: xmin, xmax
 
 END MODULE data_mod