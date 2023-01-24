MODULE data_mod

USE constants_mod
IMPLICIT NONE	
  
 
TYPE particle
    ! type data
    REAL*8 :: x = 0.d0   ! position-x
    REAL*8 :: y = 0.d0   ! position-y
    REAL*8 :: z = 0.d0   ! position-z
    REAL*8 :: ux = 0.d0  ! velocity-x
    REAL*8 :: uy = 0.d0  ! velocity-y
    REAL*8 :: uz = 0.d0  ! velocity-z
    LOGICAL :: oob = .TRUE. ! out-of-bounds flag to indicate whether particle has left the domain
END TYPE particle
 
 
REAL*8, ALLOCATABLE :: rho(:,:,:),                                          &
                       Ex(:,:,:), Ey(:,:,:), Ez(:,:,:),                     &
                       Bx(:,:,:), By(:,:,:), Bz(:,:,:),                     &
                       Ex_grid(:,:,:), Ey_grid(:,:,:), Ez_grid(:,:,:),      &
                       Bx_grid(:,:,:), By_grid(:,:,:), Bz_grid(:,:,:),      &
                       Jx(:,:,:), Jy(:,:,:), Jz(:,:,:),                     &
                       ne(:,:,:), ni(:,:,:), fv(:,:),                       &
                       buffer(:,:,:)                                        ! general purpose buffer (used mainly for digital filtering)

TYPE(particle), ALLOCATABLE :: particles(:,:)

INTEGER :: Np_in(2) = 0, Np_esc(2) = 0
INTEGER :: Ntot_inj(2) = 0, Ntot_esc(2) = 0
REAL*8  :: q(2), m(2)
REAL*8 :: dx, dy, dz, dt
REAL*8 :: Lx, Ly, Lz         
REAL*8 :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER :: i_sim, ix, iy, iz, mx, my, mz, lot

REAL*8 :: n0  
REAL*8 :: K_E, Px, Py, Pz, ESE, ESE2
REAL*8 :: vd1, vd2, vspread1, vspread2

REAL*8 :: o3, o2, o1, o, sm(27), table(64)
INTEGER :: ms(27), ie, je, ke
 
 
 
CONTAINS



SUBROUTINE create_grid_arrays()


	ALLOCATE(particles(ns,Np_max))
    ALLOCATE(rho(-1:nx+3,-1:ny+3,-1:nz+3))
    ALLOCATE(Ex(-1:nx+3,-1:ny+3,-1:nz+3), Ey(-1:nx+3,-1:ny+3,-1:nz+3),Ez(-1:nx+3,-1:ny+3,-1:nz+3))  ! -2,-1,0,n+1,n+1,n+3 are ghost cells
    ALLOCATE(Bx(-1:nx+3,-1:ny+3,-1:nz+3), By(-1:nx+3,-1:ny+3,-1:nz+3),Bz(-1:nx+3,-1:ny+3,-1:nz+3))  
    ALLOCATE(Ex_grid(-1:nx+1,-1:ny+1,-1:nz+1), Ey_grid(-1:nx+1,-1:ny+1,-1:nz+1), Ez_grid(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(Bx_grid(-1:nx+1,-1:ny+1,-1:nz+1), By_grid(-1:nx+1,-1:ny+1,-1:nz+1), Bz_grid(-1:nx+1,-1:ny+1,-1:nz+1))  
    ALLOCATE(Jx(-1:nx+3,-1:ny+3,-1:nz+3), Jy(-1:nx+3,-1:ny+3,-1:nz+3), Jz(-1:nx+3,-1:ny+3,-1:nz+3))
    ALLOCATE(ne(0:nx+1,0:ny+1,0:nz+1), ni(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(buffer(-1:nx+3,-1:ny+3,-1:nz+3))
    ALLOCATE(fv(1:ns,1:vbins))

    rho = 0.d0
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0
    Ex_grid = 0.d0
    Ey_grid = 0.d0
    Ez_grid = 0.d0
    Bx_grid = 0.d0
    By_grid = 0.d0
    Bz_grid = 0.d0
    Ex = 0.d0
    Ey = 0.d0
    Ez = 0.d0
    Bx = 0.d0
    By = 0.d0
    Bz = 0.d0
    buffer = 0.d0


END SUBROUTINE create_grid_arrays



SUBROUTINE destroy_grid_arrays()

    DEALLOCATE(particles)
    DEALLOCATE(rho, ne, ni)
    DEALLOCATE(Ex, Ey, Ez)
    DEALLOCATE(Bx, By, Bz)
    DEALLOCATE(Ex_grid, Ey_grid, Ez_grid)
    DEALLOCATE(Bx_grid, By_grid, Bz_grid)
    DEALLOCATE(Jx, Jy, Jz)
    DEALLOCATE(buffer)
    DEALLOCATE(fv)

END SUBROUTINE destroy_grid_arrays 


 
END MODULE data_mod