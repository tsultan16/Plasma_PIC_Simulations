MODULE data_mod

USE constants_mod
USE MPI

IMPLICIT NONE	
  
!INCLUDE 'mpif.h'

 
TYPE particle
    ! type data
    REAL*8 :: species    ! species type, 1=electron, 2=ion
    REAL*8 :: x = 0.d0   ! position-x
    REAL*8 :: y = 0.d0   ! position-y
    REAL*8 :: z = 0.d0   ! position-z
    REAL*8 :: vx = 0.d0  ! velocity-x
    REAL*8 :: vy = 0.d0  ! velocity-y
    REAL*8 :: vz = 0.d0  ! velocity-z
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


TYPE(particle), ALLOCATABLE :: particles(:,:)!, &
                               
                               
REAL*8, ALLOCATABLE :: bparticles_xp(:), bparticles_xm(:), &  ! arrays for out of bounds particles (separate array for each boundary face and edge)
                       bparticles_yp(:), bparticles_ym(:), &
                       bparticles_tr(:), bparticles_tl(:), &
                       bparticles_br(:), bparticles_bl(:)                       


INTEGER :: Np_in(2) = 0
INTEGER :: Np_esc_xm, Np_esc_xp, Np_esc_ym, Np_esc_yp, &
           Np_esc_bl, Np_esc_br, Np_esc_tl, Np_esc_tr 

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
 

! RGN seed
INTEGER :: s1, s2, s3


 
! MPI variables
INTEGER :: comm2d, ierr, ndim, myrank, numprocs(1), dims(2), mycoord(2), req(4), &
           coltype, status(MPI_STATUS_SIZE), source, tag, destination
INTEGER :: neighbor_rank(8), field_buffer_size_x, field_buffer_size_y, particle_buffer_size

LOGICAL :: isperiodic(2), reorder 
REAL*8, ALLOCATABLE :: particle_buffer_in(:),  &
                       field_buffer_in_xm(:), field_buffer_in_xp(:), &
                       field_buffer_in_ym(:), field_buffer_in_yp(:), &
                       field_buffer_out_xm(:), field_buffer_out_xp(:), &
                       field_buffer_out_ym(:), field_buffer_out_yp(:)
                       
                       !field_buffer_in_2(:), field_buffer_out_2(:), &
                       !field_buffer_in_3(:), field_buffer_out_3(:), &
                       !field_buffer_in_4(:), field_buffer_out_4(:)
INTEGER :: xlow, ylow, zlow, nb_xlow, nb_xhi, nb_ylow, nb_yhi 
INTEGER :: bxlow, bylow, bzlow, bxhi, byhi, bzhi, &
           exlow, eylow, ezlow, exhi, eyhi, ezhi  
 
CONTAINS



SUBROUTINE create_grid_arrays()


	ALLOCATE(particles(ns,Np_max))
    

	ALLOCATE(bparticles_xp(0:nvars_particles*max_particles), bparticles_yp(0:nvars_particles*max_particles))
	ALLOCATE(bparticles_xm(0:nvars_particles*max_particles), bparticles_ym(0:nvars_particles*max_particles))
    ALLOCATE(bparticles_tr(0:nvars_particles*max_particles), bparticles_br(0:nvars_particles*max_particles))
	ALLOCATE(bparticles_tl(0:nvars_particles*max_particles), bparticles_bl(0:nvars_particles*max_particles))
  
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
    DEALLOCATE(bparticles_xm, bparticles_xp)
    DEALLOCATE(bparticles_ym, bparticles_yp)
    DEALLOCATE(bparticles_tl, bparticles_tr)
    DEALLOCATE(bparticles_bl, bparticles_br)
    DEALLOCATE(rho, ne, ni)
    DEALLOCATE(Ex, Ey, Ez)
    DEALLOCATE(Bx, By, Bz)
    DEALLOCATE(Ex_grid, Ey_grid, Ez_grid)
    DEALLOCATE(Bx_grid, By_grid, Bz_grid)
    DEALLOCATE(Jx, Jy, Jz)
    DEALLOCATE(buffer)
    DEALLOCATE(fv)

END SUBROUTINE destroy_grid_arrays 



! Returns (double precision) random number between 0 and 1
! Adapted from Taus88 C-code in L'ECUYER(1996)
! (need to put this function in a better place...)
FUNCTION rand_taus() RESULT(r)


    REAL*8 :: r
    INTEGER :: b, c   ! Signed 64-bit INTEGER

    b = ISHFT(IEOR(ISHFT(s1,13),s1),-19)    
    c = 2147483647
    s1 = IEOR(ISHFT(IAND(s1,c),12),b)
    b = ISHFT(IEOR(ISHFT(s2,2),s2),-25)
    c = 2147483639
    s2 = IEOR(ISHFT(IAND(s2,c),4),b)
    b = ISHFT(IEOR(ISHFT(s3,3),s3),-11)
    c = 2147483631
    s3 = IEOR(ISHFT(IAND(s3,c),17),b)

    r = 0.5D0+IEOR(s1,IEOR(s2,s3))*2.3283064365D-10    ! 2.3283064365D-10 = 2^-32
    
    IF(r .GE. 1.D0 .OR. r .LT. 0.D0) THEN
        PRINT*,'Random number generator failed..r=',r
        STOP
    END IF

END FUNCTION rand_taus


! this subroutine returns an array of random numbers between [0,1.0]
SUBROUTINE RAND_NUM(p)

    REAL*8, INTENT(OUT) :: p(:) 
    INTEGER :: i
    
    DO i = 1, SIZE(p)
        p(i) = rand_taus() 
    END DO

END SUBROUTINE RAND_NUM
 
END MODULE data_mod