MODULE fieldSolver_mod

USE constants_mod
USE data_mod
USE field_bc_mod
USE domain_mod
IMPLICIT NONE


CONTAINS



! Electric and Magnetic fields are decentered across the grid. Electric field is define at cell faces while magnetic fields are defined at cell corners.
! (where each cell is a cube centered at x_i, y_j, z_k)
!
! Bx(i,j,k) = Bx(x_i, y_j+1/2, z_k+1/2)
! By(i,j,k) = By(x_i+1/2, y_j, z_k+1/2)
! Bz(i,j,k) = Bz(x_i+1/2, y_j+1/2, z_k)
! Ex(i,j,k) = Ex(x_i+1/2, y_j, z_k)
! Ey(i,j,k) = Ey(x_i, y_j+1/2, z_k)
! Ez(i,j,k) = Ez(x_i, y_j, z_k+1/2)
!
!
!   e.g. on the xy plane
!
!         Ey    Bz-
!    ----_*-----*(i+1/2,j+1/2) 
!   | (i,j+1/2) |
!   |           |
!   |     *     * Ex
!   |   (i,j)   |(i+1/2,j)
!   |           |
!    -----------
!
!     
!
! Solving Maxwell's equations in their surface-integral form (only hyperbolic ones with the time-derivatives) amounts to the following operation: 

! *********************************************************************************************************************
! * change in E (B) flux across cell-face = line integral of B (E) along boundary of that face + charge-flux (zero) *
! *********************************************************************************************************************
!
! Flux of Bx across x face:
! dy*dz*[(Bx(new)-Bx(old))_(i,j+1/2,k+1/2)] /dt = -dy*c*[ Ey_(i,j+1/2,k) -  Ey_(i,j+1/2,k+1) ] - dz*c*[ Ez_(i,j+1,k+1/2) -  Ez_(i,j,k+1/2) ]   
!
! Flux of By across y face:
! dx*dz*[(By(new)-By(old))_(i+1/2,j,k+1/2)] /dt = -dx*c*[ Ex_(i+1/2,j+1,k) - Ex_(i+1/2,j,k) ] - dz*c*[ Ez_(i,j,k+1/2) -  Ez_(i,j+1,k+1/2) ]   
!
! Flux of Bz across z face:
! dx*dy*[(Bz(new)-Bz(old))_(i+1/2,j+1/2,k)] /dt = -dx*c*[ Ex_(i+1/2,j,k) -  Ex_(i+1/2,j+1,k) ] - dy*c*[ Ey_(i+1,j+1/2,k) -  Ey_(i,j+1/2,k) ]   
! 
! Flux of Ex across x-face:
! dy*dz*[(Ex(new)-Ex(old))_(i+1/2,j,k)] /dt = dz*c*[ Bz_(i+1/2,j+1/2,k) -  Bz_(i+1/2,j-1/2,k) ] + dy*c*[ By_(i+1/2,j,k-1/2) - By_(i+1/2,j,k+1/2) ]
!                                             - dy*dz*Jx_(i+1/2,j,k)  
!
! Flux of Ey across y-face:
! dx*dz*[(Ey(new)-Ey(old))_(i,j+1/2,k)] /dt = dz*c*[ Bz_(i-1/2,j+1/2,k) -  Bz_(i+1/2,j+1/2,k) ] + dx*c*[ Bx_(i,j+1/2,k+1/2) - Bx_(i,j+1/2,k-1/2) ]
!                                            - dx*dz*Jy_(i,j+1/2,k) 
!
! Flux of Ez across z-face:
! dx*dy*[(Ez(new)-Ez(old))_(i,j,k+1/2)] /dt = dx*c*[ Bx_(i,j-1/2,k+1/2) -  Bx_(i,j+1/2,k+1/2) ] + dy*c*[ By_(i+1/2,j,k+1/2) - By_(i+1/2,j,k-1/2) ]
!                                            - dx*dz*Jz_(i,j,k+1/2) 
!
!
! Note: The quantities on the R.H.S. are time-centered, i.e. evaluated half-way betweem the "old" and "new" times on the L.H.S.
!       This means that the E- field components need to be time-decentered from the B-field and J components. 
!
!        We will use the following time-decentering:
!
!        t^n : x, E
!        t^n+1/2: v, J, B
!
!        The only issue that arises is that updating v(t^n-1/2)->v(t^n+1/2) requires B(t^n). To accomodate for this,
!        we only update B by a half time-step here. The remaining half-update needs to bve done separately.  



SUBROUTINE half_update_B()

	INTEGER :: i, j, k

    !************************************************
    ! Magnetic field half-update:  B (t-dt/2) -> B(t)
    !************************************************
    IF(bndry .EQ. 2) THEN
        ! preliminary partial update of boundary edges
        !CALL radiation_prel_B() 
    END IF    
    
    
    IF(nz .EQ. 1) THEN  

        !$OMP PARALLEL DO 
        DO k = bzlow, bzhi
            DO j = bylow, byhi
                DO i = bxlow, bxhi
                    Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * (Ez(i,j+1,k) - Ez(i,j,k))
                    By(i,j,k) = By(i,j,k) - 0.5 * c * (Ez(i,j,k) - Ez(i+1,j,k)) 
                    Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
    
    ELSE

    
        !$OMP PARALLEL DO 
        DO k = bzlow, bzhi
            DO j = bylow, byhi
                DO i = bxlow, bxhi
                    Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * ((Ey(i,j,k) - Ey(i,j,k+1)) + (Ez(i,j+1,k) - Ez(i,j,k)))
                    By(i,j,k) = By(i,j,k) - 0.5 * c * ((Ez(i,j,k) - Ez(i+1,j,k)) + (Ex(i,j,k+1) - Ex(i,j,k)))
                    Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
    
    END IF


    ! exchange with neighbors
    CALL exchange_B(1,0)
    

END SUBROUTINE half_update_B



SUBROUTINE full_update_EB()

	INTEGER :: i, j, k
   
   
    !*********************************************************
    ! Magnetic field remaining half-update: B (t) -> B(t+dt/2) 
    !*********************************************************
    
    IF(nz .EQ. 1) THEN  

        !$OMP PARALLEL DO 
        DO k = bzlow, bzhi
            DO j = bylow, byhi
                DO i = bxlow, bxhi
                    Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * (Ez(i,j+1,k) - Ez(i,j,k))
                    By(i,j,k) = By(i,j,k) - 0.5 * c * (Ez(i,j,k) - Ez(i+1,j,k)) 
                    Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
      
    ELSE

        !$OMP PARALLEL DO 
        DO k = bzlow, bzhi
            DO j = bylow, byhi
                DO i = bxlow, bxhi
                    Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * ((Ey(i,j,k) - Ey(i,j,k+1)) + (Ez(i,j+1,k) - Ez(i,j,k)))
                    By(i,j,k) = By(i,j,k) - 0.5 * c * ((Ez(i,j,k) - Ez(i+1,j,k)) + (Ex(i,j,k+1) - Ex(i,j,k)))
                    Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

    END IF    
       
    ! set boundary magnetic field 
    IF(bndry .EQ. 1) THEN
    
       CALL periodic_boundary_B()
       
    ELSE IF(bndry .EQ. 2) THEN
    
       CALL radiation_bc_B()
    
       !CALL radiation_post_B()
       
    END IF
    
    
    ! exchange with neighbors
    CALL exchange_B(1,0)
   
   
    !*********************************************************************
    ! Electric Field field update (w/o charge flux term) : E(t) -> E(t+dt)
    !*********************************************************************
 
    ! Note: E field is not updated at the boundary or the first interior layer (i,j,k = 0, n+1)  
        
    IF(bndry .EQ. 2) THEN
        ! preliminary partial update of boundary edges
        !CALL radiation_prel_E()  
    END IF 

    IF(nz .EQ. 1) THEN  

    
        !$OMP PARALLEL DO
        DO k = ezlow, ezhi
            DO j = eylow, eyhi
                DO i = exlow, exhi
                    Ex(i,j,k) = Ex(i,j,k) + c* (Bz(i,j,k) - Bz(i,j-1,k)) 
                    Ey(i,j,k) = Ey(i,j,k) + c* (Bz(i-1,j,k) - Bz(i,j,k))
                    Ez(i,j,k) = Ez(i,j,k) + c* ((Bx(i,j-1,k) - Bx(i,j,k)) + (By(i,j,k) - By(i-1,j,k))) 
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
     
    ELSE

        !$OMP PARALLEL DO
        DO k = ezlow, ezhi
            DO j = eylow, eyhi
                DO i = exlow, exhi
                    Ex(i,j,k) = Ex(i,j,k) + c* ((By(i,j,k-1) - By(i,j,k)) + (Bz(i,j,k) - Bz(i,j-1,k))) 
                    Ey(i,j,k) = Ey(i,j,k) + c* ((Bx(i,j,k) - Bx(i,j,k-1)) + (Bz(i-1,j,k) - Bz(i,j,k)))
                    Ez(i,j,k) = Ez(i,j,k) + c* ((Bx(i,j-1,k) - Bx(i,j,k)) + (By(i,j,k) - By(i-1,j,k))) 
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

    END IF    
     
    ! set boundary electric field 
    IF(bndry .EQ. 1) THEN
    
        CALL periodic_boundary_E()
        
    ELSE IF(bndry .EQ. 2) THEN
            
        CALL radiation_bc_E()
        
        !CALL radiation_post_E()
     
    END IF


    ! decrement current from electric field  
    CALL decrement_currents_E()


    ! exchange with neighbors
    CALL exchange_E(1,0)
      
   
	IF(print_debug) THEN
	PRINT*,''
	PRINT*,'Jy = '
	DO k = nz, 1, -1
        DO j = 1, ny
            WRITE(*,FMT='(1f7.1)',ADVANCE='NO') Jy(nx/2,j,k)
        END DO
        PRINT*,''
    END DO

    PRINT*,''
	PRINT*,'Ey = '
	DO k = nz, 1, -1
        DO j = 1, ny
            WRITE(*,FMT='(1f7.1)',ADVANCE='NO') Ey(nx/2,j,k)
        END DO
        PRINT*,''
    END DO

	END IF


END SUBROUTINE full_update_EB



! complete electric field update by adding in the charge flux term
SUBROUTINE decrement_currents_E()

	INTEGER :: i, j, k, ix1,ix2,iy1,iy2,iz1,iz2

    IF(mycoord(1)  .EQ. 0) THEN
        ix1 = 1
    ELSE
       ix1 = 0
    END IF

    IF(mycoord(1)  .EQ. nranks_x-1) THEN
        ix2 = 1
    ELSE
       ix2 = 0
    END IF
    
    IF(mycoord(2)  .EQ. 0) THEN
        iy1 = 1
    ELSE
       iy1 = 0
    END IF

    IF(mycoord(2)  .EQ. nranks_y-1) THEN
        iy2 = 1
    ELSE
       iy2 = 0
    END IF
    
    iz1 = 1
    iz2 = 1


    !$OMP PARALLEL DO
    DO k = 1-iz1, nz+iz2
        DO j = 1-iy1, ny+iy2
            DO i = 1-ix1, nx+ix2
                Ex(i,j,k) = Ex(i,j,k) - Jx(i,j,k)
                Ey(i,j,k) = Ey(i,j,k) - Jy(i,j,k)
                Ez(i,j,k) = Ez(i,j,k) - Jz(i,j,k)
            END DO
        END DO
     END DO
    !$OMP END PARALLEL DO


END SUBROUTINE decrement_currents_E


! neutralizes excess boundary charges by smoothing out surface currents
! (i.e. the boundaries are made slightly conductive)
SUBROUTINE neutralize_boundary_charge()

      ! x+ and x- faces
      !call boundary_conduct(ix, iy, iz, mx, my, mz, Ey, Ez, 0.25d0)
      ! y+ and y- faces
      !call boundary_conduct(iy, iz, ix, my, mz, mx, Ez, Ex, 0.25d0)
      ! z+ and z- faces
      !call boundary_conduct(iz, ix, iy, mz, mx, my, Ex, Ey, 0.25d0)
      
     IF(bndry .EQ. 2) THEN 
        CALL boundary_conduct_x(0.25d0)
        CALL boundary_conduct_y(0.25d0)
        CALL boundary_conduct_z(0.25d0)
   END IF
      
END SUBROUTINE neutralize_boundary_charge


SUBROUTINE grid_avg_fields()

    INTEGER :: i, j, k
   
    ! grid-point averages of em fields, i.e. Ex_j,k, Ey_j,k, Bz_j,k
    
    !$OMP PARALLEL DO
    DO k = 0, nz+1
        DO j = 0, ny+1
            DO i = 0, nx+1
                Bx_grid(i,j,k) = 0.25d0*( Bx(i,j-1,k-1) + Bx(i,j,k) + Bx(i,j-1,k) + Bx(i,j,k-1) )
                By_grid(i,j,k) = 0.25d0*( By(i-1,j,k-1) + By(i,j,k) + By(i-1,j,k) + By(i,j,k-1) )
                Bz_grid(i,j,k) = 0.25d0*( Bz(i-1,j-1,k) + Bz(i,j,k) + Bz(i-1,j,k) + Bz(i,j-1,k) )
                
                Ex_grid(i,j,k) = 0.5d0*( Ex(i-1,j,k) + Ex(i,j,k) )
                Ey_grid(i,j,k) = 0.5d0*( Ey(i,j-1,k) + Ey(i,j,k) )  
                Ez_grid(i,j,k) = 0.5d0*( Ez(i,j,k-1) + Ez(i,j,k) )                  
            END DO
        END DO    
    END DO    
    !$OMP END PARALLEL DO
  
END SUBROUTINE grid_avg_fields



! Note: For radiation boundary conditions, not to carefully treat 
! boundaries that are split among different domains. Currently
! edge bc routines not working correctly in 2d domain decomposition

SUBROUTINE radiation_bc_E()

    ! apply radiation b.c. for domains that touch the x and y boundaries 
    IF(mycoord(1)  .EQ. 0) CALL radiation_e_xm()
    
    IF(ndims .EQ. 2) THEN
        IF(mycoord(2)  .EQ. 0) CALL radiation_e_ym()
    ELSE
        CALL radiation_e_ym()
    END IF
    
    ! all domains contain part of the zm boundary
    CALL radiation_e_zm()
       

END SUBROUTINE radiation_bc_E



SUBROUTINE radiation_prel_E()

    ! apply radiation b.c. for domains that touch the x and y boundaries 
    IF(mycoord(1)  .EQ. 0)  THEN
        CALL radiation_boundary_prel_edge(Ey, Ez, Ex, By, Bz, Bx, -iy, -iz, -ix, my, mz, mx, lot)
    END IF
    
    IF(ndims .EQ. 2) THEN
        IF(mycoord(2)  .EQ. 0) THEN 
        CALL radiation_boundary_prel_edge(Ez, Ex, Ey, Bz, Bx, By, -iz, -ix, -iy, mz, mx, my, lot)
        END IF    
    ELSE
        CALL radiation_boundary_prel_edge(Ez, Ex, Ey, Bz, Bx, By, -iz, -ix, -iy, mz, mx, my, lot)
    END IF
    
    ! all domains contain part of the zm boundary
    CALL radiation_boundary_prel_edge(Ex, Ey, Ez, Bx, By, Bz, -ix, -iy, -iz, mx, my, mz, lot)         
       

END SUBROUTINE radiation_prel_E


SUBROUTINE radiation_post_E()

    ! apply radiation b.c. for domains that touch the x and y boundaries 
    IF(mycoord(1)  .EQ. 0)  THEN
        CALL radiation_boundary_post_edge(Ey, Ez, Ex, By, Bz, Bx, -iy, -iz, -ix, my, mz, mx, lot)
    END IF
    
    IF(ndims .EQ. 2) THEN
        IF(mycoord(2)  .EQ. 0) THEN 
        CALL radiation_boundary_post_edge(Ez, Ex, Ey, Bz, Bx, By, -iz, -ix, -iy, mz, mx, my, lot)
        END IF    
    ELSE
        CALL radiation_boundary_post_edge(Ez, Ex, Ey, Bz, Bx, By, -iz, -ix, -iy, mz, mx, my, lot)
    END IF
    
    ! all domains contain part of the zm boundary
    CALL radiation_boundary_post_edge(Ex, Ey, Ez, Bx, By, Bz, -ix, -iy, -iz, mx, my, mz, lot)  
       
      
END SUBROUTINE radiation_post_E


SUBROUTINE radiation_bc_B()

    ! apply radiation b.c. for domains that touch the x and y boundaries 
    IF(mycoord(1)  .EQ. nranks_x-1) CALL radiation_b_xp()
    
    IF(ndims .EQ. 2) THEN
        IF(mycoord(2)  .EQ. nranks_y-1) CALL radiation_b_yp()
    ELSE
        CALL radiation_b_yp()
    END IF

    ! all domains contain part of the zp boundary
    CALL radiation_b_zp()
       

END SUBROUTINE radiation_bc_B



!*****************************
! Field Boundary Conditions
!*****************************

SUBROUTINE radiation_prel_B()

    ! apply radiation b.c. for domains that touch the x and y boundaries 
    IF(mycoord(1)  .EQ. nranks_x-1) THEN 
        CALL radiation_boundary_prel_edge(By, Bz, Bx, Ey, Ez, Ex, iy, iz, ix, my, mz, mx, 1)
    END IF
    
    IF(ndims .EQ. 2) THEN
        IF(mycoord(2)  .EQ. nranks_y-1) THEN
        CALL radiation_boundary_prel_edge(Bz, Bx, By, Ez, Ex, Ey, iz, ix, iy, mz, mx, my, 1)
        END IF    
    ELSE
        CALL radiation_boundary_prel_edge(Bz, Bx, By, Ez, Ex, Ey, iz, ix, iy, mz, mx, my, 1)
    END IF

    
    ! all domains contain part of the zp boundary
    CALL radiation_boundary_prel_edge(Bx, By, Bz, Ex, Ey, Ez, ix, iy, iz, mx, my, mz, 1)    

       
        
END SUBROUTINE radiation_prel_B



SUBROUTINE radiation_post_B()

    ! apply radiation b.c. for domains that touch the x and y boundaries 
    IF(mycoord(1)  .EQ. nranks_x-1) THEN 
        CALL radiation_boundary_post_edge(By, Bz, Bx, Ey, Ez, Ex, iy, iz, ix, my, mz, mx, 1)
    END IF
    
    IF(ndims .EQ. 2) THEN
        IF(mycoord(2)  .EQ. nranks_y-1) THEN
        CALL radiation_boundary_post_edge(Bz, Bx, By, Ez, Ex, Ey, iz, ix, iy, mz, mx, my, 1)
        END IF    
    ELSE
        CALL radiation_boundary_post_edge(Bz, Bx, By, Ez, Ex, Ey, iz, ix, iy, mz, mx, my, 1)
    END IF

    
    ! all domains contain part of the zp boundary
    CALL radiation_boundary_post_edge(Bx, By, Bz, Ex, Ey, Ez, ix, iy, iz, mx, my, mz, 1)    

       
        
END SUBROUTINE radiation_post_B



! Radiation absorbing boundaries (Lindman method, outgoing plane waves, lowest order)

! This subroutine assumes radiation boundary conditions for the magnetic field surface components at the
! z+ boundary. It can be reused for the other boundaries (x+,y+, x-, etc.) and E-field components, by manipualting the ordering of
! the input arrays that are passed and passing the appropriate index strides. The inputs arrays will be assumed to be 1D.
! ( This subroutine also assumes dx = dy = dz = dt = 1 )
SUBROUTINE radiation_boundary_surface(bx, by, bz, ex, ey, ez, ix, iy, iz, mx, my, mz, m00)

    REAL*8, INTENT(INOUT) :: bx(1), by(1), bz(1)
    REAL*8, INTENT(IN) :: ex(1), ey(1), ez(1) 
    INTEGER, INTENT(IN) :: ix, iy, iz, mx, my, mz, m00
    REAL*8 :: rs, s, os
    INTEGER :: m, n
       
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
     
    
    DO m = m00 + iz * (mz-1), m00 + iz * (mz-1) + iy * (my-2), iy
        DO n = m, m + ix * (mx-2), ix
            bz(n) =    bz(n) + 0.5d0 * c * (ex(n+iy) - ex(n) - ey(n+ix) + ey(n))  
        END DO
        
        DO n = m + ix, m + ix * (mx-2), ix
            bx(n) = bx(n) + rs * (bx(n-iz) - bx(n) + s * (bz(n) - bz(n-ix))) - &     
                    os * (ez(n+iy) - ez(n)) - (os - c) * (ez(n+iy-iz) - ez(n-iz)) - &
                    c * (ey(n) - ey(n-iz))
        END DO
    END DO  

    DO m = m00 + iz * (mz-1), m00 + iz * (mz-1) + ix * (mx-2), ix
         DO n = m + iy, m + iy * (my-2), iy
            by(n) = by(n) + rs * (by(n-iz)  -by(n) + s * (bz(n) - bz(n-iy))) + &
                    os * (ez(n+ix) - ez(n)) + (os - c) *(ez(n+ix-iz) - ez(n-iz)) + &
                    c * (ex(n) - ex(n-iz))
        END DO
        
        DO n = m, m + iy * (my-2), iy
            bz(n) = bz(n) + 0.5d0 * c * (ex(n+iy) - ex(n) - ey(n+ix) + ey(n))
        END DO
        
    END DO


END SUBROUTINE radiation_boundary_surface



SUBROUTINE radiation_boundary_edge(bx, ix, iy, iz, mx, my, mz, m00)

    REAL*8, INTENT(INOUT) :: bx(1)
    INTEGER, INTENT(IN) :: ix, iy, iz, mx, my, mz, m00
    INTEGER :: lx, ly, lz, i, j, k, n
    
    lx = ix*(mx-1)
    ly = iy*(my-1)
    lz = iz*(mz-1)
    n = m00 + iy + lz
    
    
      DO j = 2,my-1
        bx(n) = bx(n+ix)+bx(n-iz)-bx(n+ix-iz)
        bx(n+lx) = bx(n+lx-ix)+bx(n+lx-iz)-bx(n+lx-ix-iz)
        n = n+iy
      END DO  
  
      n=m00+ly+iz

      DO k = 2, mz-1
          bx(n) = bx(n+ix)+bx(n-iy)-bx(n+ix-iy)
          bx(n+lx) = bx(n+lx-ix)+bx(n+lx-iy)-bx(n+lx-ix-iy)
          n = n+iz
      END DO  

      n=m00+ly+lz

       DO i = 1, mx
          bx(n) = bx(n-iy)+bx(n-iz)-bx(n-iy-iz)
          bx(n-ly) = bx(n-ly+iy)+bx(n-ly-iz)-bx(n-ly+iy-iz)
          bx(n-lz) = bx(n-lz-iy)+bx(n-lz+iz)-bx(n-lz-iy+iz)
          n = n+ix
      END DO
  
END SUBROUTINE radiation_boundary_edge  


  
SUBROUTINE radiation_boundary_prel_edge(bx, by, bz, ex, ey, ez, ix, iy, iz, mx, my, mz, m)

    REAL*8, INTENT(INOUT) :: bx(1), by(1), bz(1)
    REAL*8, INTENT(IN) :: ex(1), ey(1), ez(1) 
    INTEGER, INTENT(IN) :: ix, iy, iz, mx, my, mz, m
    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: n

    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
      
    DO n = m + iy * (my-1) + iz * (mz-1) + ix, m + iy * (my-1) + iz * (mz-1) + ix * (mx-2), ix
        bx(n) = bx(n-iy-iz) + (1.d0-4.d0*t)*bx(n) + (1.d0-2.d0*t)*(bx(n-iy) + bx(n-iz)) + &
                s * t * (by(n) - by(n-ix) + by(n-iz) - by(n-ix-iz) + &
                bz(n) - bz(n-ix) + bz(n-iy) - bz(n-ix-iy))
    END DO

    DO n = m + iz * (mz-1), m + iz * (mz-1) + iy * (my-2), iy
        bx(n) = (1.d0-c*r)*(bx(n)+bx(n+ix))+bx(n-iz)+bx(n+ix-iz)-r*(bz(n)+   &
                c*((1.d0-s)*(ex(n+iy)-ex(n))+(1.+s)*0.25d0*(ez(n+iy)-ez(n)+ez(n+ix+iy) - &
                ez(n+ix)+ez(n+iy-iz)-ez(n-iz)+ez(n+ix+iy-iz)-ez(n+ix-iz))))
    END DO
    
    DO n = m+iy*(my-1), m+iy*(my-1)+iz*(mz-2), iz
        bx(n) = (1.d0 - c*r) * (bx(n) + bx(n+ix)) + bx(n-iy) + bx(n+ix-iy) - & 
                r * (by(n) - c*((1.d0-s) * (ex(n+iz) - ex(n)) + &
                (1.d0+s)*0.25d0 * (ey(n+iz) - ey(n) + ey(n+ix+iz) - &
                ey(n+ix) + ey(n+iz-iy) - ey(n-iy) + ey(n+ix+iz-iy) - ey(n+ix-iy))))
    END DO
    
    DO n = m+iy*(my-1)+iz*(mz-1), m+ix*(mx-2)+iy*(my-1)+iz*(mz-1), ix
       temp = bz(n) - 0.5d0 * c * (1.d0 - s) * (ey(n+ix) - ey(n) + ey(n+ix-iy) - ey(n-iy))
       bz(n) = bz(n-iy) - bz(n) + p * temp + q * by(n)
       by(n) = by(n-iz) - by(n) + p * by(n) + q * temp
   END DO
    

END SUBROUTINE radiation_boundary_prel_edge


SUBROUTINE radiation_boundary_post_edge(bx, by, bz, ex, ey, ez, ix, iy, iz, mx, my, mz, m)

    REAL*8, INTENT(INOUT) :: bx(1), by(1), bz(1)
    REAL*8, INTENT(IN) :: ex(1), ey(1), ez(1) 
    INTEGER, INTENT(IN) :: ix, iy, iz, mx, my, mz, m
    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: n

      
    s = 0.4142136
    p = (1.d0 + c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    t = c / (2.d0 * c + 1.d0 + s)
    r = 4.d0 / (2.d0 * c + 2.d0 + s)

            
    DO n = m+iy*(my-1)+iz*(mz-1)+ix, m+iy*(my-1)+iz*(mz-1)+ix*(mx-2), ix
        bx(n) = bx(n) - (1.d0-4.d0*t)*bx(n-iy-iz) - (1.d0-2.d0*t)*(bx(n-iy)+bx(n-iz)) + &
                s * t * (by(n) - by(n-ix) + by(n-iz) - by(n-ix-iz) + &
                bz(n) - bz(n-ix) + bz(n-iy) - bz(n-ix-iy))
    END DO
    
    DO n = m+iz*(mz-1), m+iz*(mz-1)+iy*(my-2), iy
        bx(n) = bx(n) - bx(n+ix) - (1.d0 - c*r) * (bx(n-iz) + bx(n+ix-iz)) + r * bz(n)
    END DO
    
    DO n = m+iy*(my-1), m+iy*(my-1)+iz*(mz-2), iz
        bx(n) = bx(n) - bx(n+ix) - (1.d0 - c*r) * (bx(n-iy) + bx(n+ix-iy)) + r * by(n)
    END DO

    DO n = m+iy*(my-1)+iz*(mz-1), m+ix*(mx-2)+iy*(my-1)+iz*(mz-1), ix
        temp = by(n-iz)-.5*c*(1.-s)*(ez(n+ix)-ez(n)+ez(n+ix-iz)-ez(n-iz))
        bz(n) = bz(n-iy) + bz(n) - q * temp - p * bz(n-iy)
        by(n) = by(n-iz) + by(n) - q * bz(n-iy)- p * temp
    END DO
    
END SUBROUTINE radiation_boundary_post_edge


! conducting boundary at x+ face
SUBROUTINE boundary_conduct(ix, iy, iz, mx, my, mz, ey, ez, s)

    REAL*8, INTENT(INOUT) :: ey(1), ez(1)
    REAL*8, INTENT(IN) :: s
    INTEGER, INTENT(IN) :: ix, iy, iz, mx, my, mz
    INTEGER :: i, j, k
    REAL*8 :: temp

    DO i = 1 + 2*ix, 1 + ix*(mx-3), ix*(mx-5)
        DO j = i + 2*iy, i + iy*(my-4), iy
            DO k = j + 2*iz, j + iz*(mz-4), iz
            
                temp=.25*s*ey(k)
                ey(k-ix)=ey(k-ix)-temp
                ey(k+ix)=ey(k+ix)-temp
                ey(k)=ey(k)-2.*temp
      
                temp=.25*s*ez(k)
                ez(k-ix)=ez(k-ix)-temp
                ez(k+ix)=ez(k+ix)-temp
                ez(k)=ez(k)-2.*temp
   
            END DO
        END DO    
    END DO

END SUBROUTINE boundary_conduct



END MODULE fieldSolver_mod