MODULE data_mod
 
IMPLICIT NONE	
 
TYPE particle
    ! type data
    REAL*8 :: q = 0.d0   ! charge
    REAL*8 :: m = 0.d0   ! mass
    REAL*8 :: x = 0.d0   ! position-x
    REAL*8 :: y = 0.d0   ! position-y
    REAL*8 :: z = 0.d0   ! position-z
    REAL*8 :: ux = 0.d0  ! velocity-x
    REAL*8 :: uy = 0.d0  ! velocity-y
    REAL*8 :: uz = 0.d0  ! velocity-z
    LOGICAL :: oob = .FALSE. ! out-of-bounds flag to indicate whether particle has left the domain
END TYPE particle
 
 
REAL*8, ALLOCATABLE :: rho(:,:,:), phi(:,:,:),                          &
                       Ex(:,:,:), Ey(:,:,:), Ez(:,:,:),                 &
                       Bx(:,:,:), By(:,:,:), Bz(:,:,:),                 &
                       Ex_grid(:,:,:), Ey_grid(:,:,:), Ez_grid(:,:,:),  &
                       Bx_grid(:,:,:), By_grid(:,:,:), Bz_grid(:,:,:),  &
                       Jx(:,:,:), Jy(:,:,:), Jz(:,:,:),                 &
                       ne(:,:,:), ni(:,:,:), fv(:,:),                   &
                       buffer(:,:,:),                                   & ! general purpose buffer (used mainly for digital filtering)
                       pos_buffer(:,:,:)                                  ! particle positions buffer, needed during current deposition

TYPE(particle), ALLOCATABLE :: particles(:,:)

REAL*8 :: dx, dy, dz, dt
REAL*8 :: Lx, Ly, Lz         
REAL*8 :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER :: ix, iy, iz, mx, my, mz, lot

REAL*8 :: n0  
REAL*8 :: KE, Px, Py, Pz, ESE, ESE2
REAL*8 :: vd1, vd2, vspread1, vspread2

 
 END MODULE data_mod