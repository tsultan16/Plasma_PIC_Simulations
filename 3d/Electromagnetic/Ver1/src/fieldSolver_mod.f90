MODULE fieldSolver_mod

USE constants_mod
USE data_mod
IMPLICIT NONE

REAL*8 :: bound_bufferx(-1:ny+3), bound_buffery(-1:nx+3)

CONTAINS



! Electric and Magnetic fields are decentered across the grid:
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
    DO k = -1, nz+2
        DO j = -1, ny+2
            DO i = -1, nx+2
                 Bx(i,j,k) = Bx(i,j,k) - 0.5 * dt * c * ((Ey(i,j,k) - Ey(i,j,k+1)) / dz  + (Ez(i,j+1,k) - Ez(i,j,k)) / dy )
                 By(i,j,k) = By(i,j,k) - 0.5 * dt * c * ((Ez(i,j,k) - Ez(i+1,j,k)) / dx  + (Ex(i,j,k+1) - Ex(i,j,k)) / dz )
                 Bz(i,j,k) = Bz(i,j,k) - 0.5 * dt * c * ((Ex(i,j,k) - Ex(i,j+1,k)) / dy  + (Ey(i+1,j,k) - Ey(i,j,k)) / dx )
            END DO
        END DO
    END DO
   
      


END SUBROUTINE half_update_B



SUBROUTINE full_update_EB()

	INTEGER :: i, j, k
   
   
    !*********************************************************
    ! Magnetic field remaining half-update: B (t) -> B(t+dt/2) 
    !*********************************************************
    IF(bndry .EQ. 2) THEN
        ! store values of last interior layer, required later for boundary update
        !bound_bufferx(-1:ny+3) = Bz(nx+2,-1:ny+3)   
        !bound_buffery(-1:nx+3) = Bz(-1:nx+3,ny+2)
    END IF
    
    DO k = -1, nz+2
        DO j = -1, ny+2
            DO i = -1, nx+2
                 Bx(i,j,k) = Bx(i,j,k) - 0.5 * dt * c * ((Ey(i,j,k) - Ey(i,j,k+1)) / dz  + (Ez(i,j+1,k) - Ez(i,j,k)) / dy )
                 By(i,j,k) = By(i,j,k) - 0.5 * dt * c * ((Ez(i,j,k) - Ez(i+1,j,k)) / dx  + (Ex(i,j,k+1) - Ex(i,j,k)) / dz )
                 Bz(i,j,k) = Bz(i,j,k) - 0.5 * dt * c * ((Ex(i,j,k) - Ex(i,j+1,k)) / dy  + (Ey(i+1,j,k) - Ey(i,j,k)) / dx )
            END DO
        END DO
     END DO
   
   
    ! set boundary magnetic field 
    ! Bx,By,Bz on x+ (i=nx+2), y+ (j=ny+2) and z+ (k=nz+2) edges
    IF(bndry .EQ. 1) THEN
       !CALL periodic_boundary_B()
    ELSE IF(bndry .EQ. 2) THEN
       !CALL radiation_boundary_B_x()
       !CALL radiation_boundary_B_y()
       !CALL radiation_boundary_B_z()
    END IF
    
   
    !**************************************************
    ! Electric Field field full update: E(t) -> E(t+dt)
    !**************************************************
    DO k = 0, nz+3
        DO j = 0, ny+3
            DO i = 0, nx+3
                Ex(i,j,k) = Ex(i,j,k) + dt * c* ( (By(i,j,k-1) - By(i,j,k))  /dz  + (Bz(i,j,k) - Bz(i,j-1,k)) / dy  ) &
                            - dt * Jx(i,j,k)
                Ey(i,j,k) = Ey(i,j,k) + dt * c* ( (Bx(i,j,k) - Bx(i,j,k-1))  /dz  + (Bz(i-1,j,k) - Bz(i,j,k)) / dx  ) &
                            - dt * Jy(i,j,k)
                Ez(i,j,k) = Ez(i,j,k) + dt * c* ( (Bx(i,j-1,k) - Bx(i,j,k))  /dy  + (By(i,j,k) - By(i-1,j,k)) / dx  ) &
                            - dt * Jz(i,j,k)
            END DO
        END DO
     END DO

    
    ! set boundary electric field 
    ! Ex,Ey,Ez on x- (i=-1), y- (j=-1) and z- (k=-1) edges
    IF(bndry .EQ. 1) THEN
       !CALL periodic_boundary_E()
    ELSE IF(bndry .EQ. 2) THEN
       !CALL radiation_boundary_E_x()
       !CALL radiation_boundary_E_y()
       !CALL radiation_boundary_E_z()
    END IF
  
      
   
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



SUBROUTINE grid_avg_fields()

   INTEGER :: i, j, k
   
   ! grid-point averages of em fields, i.e. Ex_j,k, Ey_j,k, Bz_j,k
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
    
  
END SUBROUTINE grid_avg_fields


!*****************************
! Field Boundary Conditions
!*****************************

! Periodic boundaries
SUBROUTINE periodic_boundary_B()

INTEGER :: i
    
    DO i = 0, 1
        Bx(-i,:,:) = Bx(nx-i,:,:)
        Bx(:,-i,:) = Bx(:,ny-i,:)
        Bx(:,:,-i) = Bx(:,:,nz-i)
    
        Bx(nx+1+i,:,:) = Bx(1+i,:,:)
        Bx(:,ny+1+i,:) = Bx(:,1+i,:)
        Bx(:,:,nz+1+i) = Bx(:,:,1+i)
    
        By(-i,:,:) = By(nx-i,:,:)
        By(:,-i,:) = By(:,ny-i,:)
        By(:,:,-i) = By(:,:,nz-i)
    
        By(nx+1+i,:,:) = By(1+i,:,:)
        By(:,ny+1+i,:) = By(:,1+i,:)
        By(:,:,nz+1+i) = By(:,:,1+i)
    
        Bz(-i,:,:) = Bz(nx-i,:,:)
        Bz(:,-i,:) = Bz(:,ny-i,:)
        Bz(:,:,-i) = Bz(:,:,nz-i)
    
        Bz(nx+1+i,:,:) = Bz(1+i,:,:)
        Bz(:,ny+1+i,:) = Bz(:,1+i,:)
        Bz(:,:,nz+1+i) = Bz(:,:,1+i)      
    END DO
    
    

END SUBROUTINE periodic_boundary_B


SUBROUTINE periodic_boundary_E()

    INTEGER :: i
    
    DO i = 0, 1  
        Ex(-i,:,:) = Ex(nx-i,:,:)
        Ex(:,-i,:) = Ex(:,ny-i,:)
        Ex(:,:,-i) = Ex(:,:,nz-i)
    
        Ex(nx+1+i,:,:) = Ex(1+i,:,:)
        Ex(:,ny+1+i,:) = Ex(:,1+i,:)
        Ex(:,:,nz+1+i) = Ex(:,:,1+i)
    
        Ey(-i,:,:) = Ey(nx-i,:,:)
        Ey(:,-i,:) = Ey(:,ny-i,:)
        Ey(:,:,-i) = Ey(:,:,nz-i)
    
        Ey(nx+1+i,:,:) = Ey(1+i,:,:)
        Ey(:,ny+1+i,:) = Ey(:,1+i,:)
        Ey(:,:,nz+1+i) = Ey(:,:,1+i)
    
        Ez(-i,:,:) = Ez(nx-i,:,:)
        Ez(:,-i,:) = Ez(:,ny-i,:)
        Ez(:,:,-i) = Ez(:,:,nz-i)
    
        Ez(nx+1+i,:,:) = Ez(1+i,:,:)
        Ez(:,ny+1+i,:) = Ez(:,1+i,:)
        Ez(:,:,nz+1+i) = Ez(:,:,1+i)
    END DO
    
    
END SUBROUTINE periodic_boundary_E



! Radiation absorbing boundaries (Lindman method, outgoing plane waves, lowest order)

SUBROUTINE radiation_boundary_B_x()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, cdtdz, a, b
    
    cdtdx = c*dt/dx
    cdtdy = c*dt/dy
    cdtdz = c*dt/dz
     
    !  First, half-update the normal component of B at the x+ boundary
    DO k = -1, nz+2
        DO j = -1, ny+2
             Bx(nx+3,j,k) = Bx(nx+3,j,k) - 0.5 * dt * c * ((Ey(nx+3,j,k) - Ey(nx+3,j,k+1)) / dz + &
                            (Ez(nx+3,j+1,k) - Ez(nx+3,j,k)) / dy )
        END DO     
    END DO    

    a = 1.d0 / (1.d0 + cdtdx)
    b = (1.d0 - cdtdx)*a
    
    ! now update surface components according to the lowest order Lindman boundary conditions 
    ! (Note: The edges aren't included and have to be updated separately.)    
    DO k = 0, nz+2
        DO j = 0, ny+2
             By(nx+3,j,k) = By(nx+2,j,k) + b * (By(nx+3,j,k)-By(nx+2,j,k)) + &
                            2.d0 * cB * cdtdy * a * (Bx(nx+3,j,k) - Bx(nx+3,j-1,k)) - &
                            cE * cdtdz * a * (Ex(nx+3,j,k+1) + Ex(nx+2,j,k+1) - Ex(nx+3,j,k) - Ex(nx+2,j,k)) + & ! first three terms are from Lindman boundary conditions
                            dt * c * ((Ez(nx+2,j,k) - Ez(nx+3,j,k)) / dx  + (Ex(nx+2,j,k+1) - Ex(nx+2,j,k)) / dz )  ! this last term is to undo the regular full update
  
             
             Bz(nx+3,j,k) = Bz(nx+2,j,k) + b * (Bz(nx+3,j,k)-Bz(nx+2,j,k)) + &
                            2.d0 * cB * cdtdz * a * (Bx(nx+3,j,k) - Bx(nx+3,j,k-1)) + &
                            cE * cdtdy * a * (Ex(nx+3,j+1,k) + Ex(nx+2,j+1,k) - Ex(nx+3,j,k) - Ex(nx+2,j,k)) + &
                            dt * c * ((Ex(nx+2,j,k) - Ex(nx+2,j+1,k)) / dy  + (Ey(nx+3,j,k) - Ey(nx+2,j,k)) / dx )  ! this last term is to undo the regular full update
 
        END DO
    END DO    

    
     ! remaining half-update the normal component of B at the x+ boundary
    DO k = -1, nz+2
        DO j = -1, ny+2
             Bx(nx+3,j,k) = Bx(nx+3,j,k) - 0.5 * dt * c * ((Ey(nx+3,j,k) - Ey(nx+3,j,k+1)) / dz + &
                            (Ez(nx+3,j+1,k) - Ez(nx+3,j,k)) / dy )
        END DO     
    END DO    

    

END SUBROUTINE radiation_boundary_B_x



SUBROUTINE radiation_boundary_B_y()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, cdtdz, a, b
    
    cdtdx = c*dt/dx
    cdtdy = c*dt/dy
    cdtdz = c*dt/dz
     
    !  First, half-update the normal component of B at the y+ boundary
    DO k = -1, nz+2
        DO i = -1, nx+2
                 By(i,ny+3,k) = By(i,ny+3,k) - 0.5 * dt * c * ((Ez(i,ny+3,k) - Ez(i+1,ny+3,k)) / dx + &
                                (Ex(i,ny+3,k+1) - Ex(i,ny+3,k)) / dz )
        END DO     
    END DO    

    a = 1.d0 / (1.d0 + cdtdy)
    b = (1.d0 - cdtdy)*a
    
    ! now update surface components according to the lowest order Lindman boundary conditions 
    ! (Note: The edges aren't included and have to be updated separately.)    
    DO k = 0, nz+2
        DO i = 0, nx+2
             Bz(i,ny+3,k) = Bz(i,ny+2,k) + b * (Bz(i,ny+3,k)-Bz(i,ny+2,k)) + &
                            2.d0 * cB * cdtdz * a * (By(i,ny+3,k) - By(i,ny+3,k-1)) - &
                            cE * cdtdx * a * (Ey(i+1,ny+3,k) + Ey(i+1,ny+2,k) - Ey(i,ny+3,k) - Ey(i,ny+2,k)) + & ! first three terms are from Lindman boundary conditions
                            dt * c * ((Ex(i,ny+2,k) - Ex(i,ny+3,k)) / dy  + (Ey(i+1,ny+2,k) - Ey(i,ny+2,k)) / dx )  ! this last term is to undo the regular full update
  
             
             Bx(i,ny+3,k) = Bx(i,ny+2,k) + b * (Bx(i,ny+3,k)-Bx(i,ny+2,k)) + &
                            2.d0 * cB * cdtdx * a * (By(i,ny+3,k) - By(i-1,ny+3,k)) + &
                            cE * cdtdz * a * (Ey(i,ny+3,k+1) + Ey(i,ny+2,k+1) - Ey(i,ny+3,k) - Ey(i,ny+2,k)) + &
                            dt * c * ((Ey(i,ny+2,k) - Ey(i,ny+2,k+1)) / dz  + (Ez(i,ny+3,k) - Ez(i,ny+2,k)) / dy ) ! this last term is to undo the regular full update
 
        END DO
    END DO    

    
     ! remaining half-update the normal component of B at the y+ boundary
    DO k = -1, nz+2
        DO i = -1, nx+2
                 By(i,ny+3,k) = By(i,ny+3,k) - 0.5 * dt * c * ((Ez(i,ny+3,k) - Ez(i+1,ny+3,k)) / dx + &
                                (Ex(i,ny+3,k+1) - Ex(i,ny+3,k)) / dz )
        END DO     
    END DO   

    

END SUBROUTINE radiation_boundary_B_y



SUBROUTINE radiation_boundary_B_z()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, cdtdz, a, b
    
    cdtdx = c*dt/dx
    cdtdy = c*dt/dy
    cdtdz = c*dt/dz
     
    !  First, half-update the normal component of B at the z+ boundary
    DO j = -1, ny+2
        DO i = -1, nx+2
             Bz(i,j,nz+3) = Bz(i,j,nz+3) - 0.5 * dt * c * ((Ex(i,j,nz+3) - Ex(i,j+1,nz+3)) / dy + &
                            (Ey(i+1,j,nz+3) - Ey(i,j,nz+3)) / dx )
        END DO
    END DO    

    a = 1.d0 / (1.d0 + cdtdz)
    b = (1.d0 - cdtdz)*a
    
    ! now update surface components according to the lowest order Lindman boundary conditions 
    ! (Note: The edges aren't included and have to be updated separately.)    
    DO j = 0, ny+2
        DO i = 0, nx+2
             Bx(i,j,nz+3) = Bx(i,j,nz+2) + b * (Bx(i,j,nz+3)-Bx(i,j,nz+2)) + &
                            2.d0 * cB * cdtdx * a * (Bz(i,j,nz+3) - Bz(i-1,j,nz+3)) - &
                            cE * cdtdy * a * (Ez(i,j+1,nz+3) + Ez(i,j+1,nz+2) - Ez(i,j,nz+3) - Ez(i,j,nz+2)) + & ! first three terms are from Lindman boundary conditions
                            dt * c * ((Ey(i,j,nz+2) - Ey(i,j,nz+3)) / dz  + (Ez(i,j+1,nz+2) - Ez(i,j,nz+2)) / dy )  ! this last term is to undo the regular full update
  
             
             By(i,j,nz+3) = By(i,j,nz+2) + b * (By(i,j,nz+3)-By(i,j,nz+2)) + &
                            2.d0 * cB * cdtdy * a * (Bz(i,j,nz+3) - Bz(i-1,j,nz+3)) + &
                            cE * cdtdx * a * (Ez(i+1,j,nz+3) + Ez(i+1,j,nz+2) - Ez(i,j,nz+3) - Ez(i,j,nz+2)) + &
                            dt * c * ((Ez(i,j,nz+2) - Ez(i+1,j,nz+2)) / dx  + (Ex(i,j,nz+3) - Ex(i,j,nz+2)) / dz )  ! this last term is to undo the regular full update
 
        END DO
    END DO    

    
     ! remaining half-update the normal component of B at the z+ boundary
    DO j = -1, ny+2
        DO i = -1, nx+2
             Bz(i,j,nz+3) = Bz(i,j,nz+3) - 0.5 * dt * c * ((Ex(i,j,nz+3) - Ex(i,j+1,nz+3)) / dy + &
                            (Ey(i+1,j,nz+3) - Ey(i,j,nz+3)) / dx )
        END DO
    END DO    

    

END SUBROUTINE radiation_boundary_B_z



SUBROUTINE radiation_boundary_E_x()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, cdtdz, a, b
    
    cdtdx = c*dt/dx
    cdtdy = c*dt/dy
    cdtdz = c*dt/dz
     
    !  First, half-update the normal component of E at the x- boundary
    DO k = 0, nz+3
        DO j = 0, ny+3
            Ex(-1,j,k) = Ex(-1,j,k) + 0.5 * dt * c* ( (By(-1,j,k-1) - By(-1,j,k))  /dz  + (Bz(-1,j,k) - Bz(-1,j-1,k)) / dy  ) 
        END DO     
    END DO    

    a = 1.d0 / (1.d0 - cdtdx)
    b = (1.d0 + cdtdx)*a
    
    ! now update surface components according to the lowest order Lindman boundary conditions 
    ! (Note: The edges aren't included and have to be updated separately.)    
    DO k = 0, nz+2
        DO j = 0, ny+2
             Ey(-1,j,k) = -Ey(0,j,k) + b * (Ey(-1,j,k)-Ey(0,j,k)) - &
                            2.d0 * cB * cdtdy * a * (Ex(-1,j,k) - Ex(-1,j-1,k)) + &
                            cE * cdtdz * a * (Bx(-1,j,k+1) + Bx(0,j,k+1) - Bx(-1,j,k) - Bx(0,j,k)) - & ! first three terms are from Lindman boundary conditions
                            dt * c* ( (Bx(0,j,k) - Bx(0,j,k-1))  /dz  + (Bz(-1,j,k) - Bz(0,j,k)) / dx  )   ! this last term is to undo the regular full update           
            
             Ez(-1,j,k) = -Ez(0,j,k) + b * (Ez(-1,j,k)-Ez(0,j,k)) - &
                            2.d0 * cB * cdtdz * a * (Ex(-1,j,k) - Ex(-1,j,k-1)) - &
                            cE * cdtdy * a * (Bx(-1,j+1,k) + Bx(0,j+1,k) - Bx(-1,j,k) - Bx(0,j,k)) - &
                            dt * c* ( (Bx(0,j-1,k) - Bx(0,j,k))  /dy  + (By(0,j,k) - By(-1,j,k)) / dx  )  ! this last term is to undo the regular full update
 
        END DO
    END DO    

    
    ! remaining half-update the normal component of E at the x- boundary
    DO k = 0, nz+3
        DO j = 0, ny+3
            Ex(-1,j,k) = Ex(-1,j,k) + 0.5 * dt * c* ( (By(-1,j,k-1) - By(-1,j,k))  /dz  + (Bz(-1,j,k) - Bz(-1,j-1,k)) / dy  ) 
        END DO     
    END DO    
    
    
END SUBROUTINE radiation_boundary_E_x


SUBROUTINE radiation_boundary_E_y()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, cdtdz, a, b
    
    cdtdx = c*dt/dx
    cdtdy = c*dt/dy
    cdtdz = c*dt/dz
     
    !  First, half-update the normal component of E at the y- boundary
    DO k = 0, nz+3
        DO i = 0, nx+3
            Ey(i,-1,k) = Ey(i,-1,k) + 0.5 * dt * c* ( (Bx(i,-1,k) - Bx(i,-1,k-1))  /dz  + (Bz(i-1,-1,k) - Bz(i,-1,k)) / dx  )      
        END DO     
    END DO    

    a = 1.d0 / (1.d0 - cdtdy)
    b = (1.d0 + cdtdy)*a
    
    ! now update surface components according to the lowest order Lindman boundary conditions 
    ! (Note: The edges aren't included and have to be updated separately.)    
    DO k = 0, nz+2
        DO i = 0, nx+2
             Ez(i,-1,k) = -Ez(i,0,k) + b * (Ez(i,-1,k)-Ez(i,0,k)) - &
                            2.d0 * cB * cdtdz * a * (Ey(i,-1,k) - Ey(i,-1,k-1)) + &
                            cE * cdtdx * a * (By(i+1,-1,k) + By(i+1,0,k) - By(i,-1,k) - By(i,0,k)) - & ! first three terms are from Lindman boundary conditions
                            dt * c* ( (Bx(i,-1,k) - Bx(i,0,k))  /dy  + (By(i,0,k) - By(i-1,0,k)) / dx  ) ! this last term is to undo the regular full update
  
             
             Ex(i,-1,k) = -Ex(i,0,k) + b * (Ex(i,-1,k)-Ex(i,0,k)) - &
                            2.d0 * cB * cdtdx * a * (Ey(i,-1,k) - Ey(i-1,-1,k)) - &
                            cE * cdtdz * a * (By(i,-1,k+1) + By(i,0,k+1) - By(i,-1,k) - By(i,0,k)) - &
                            dt * c* ( (By(i,j,-1) - By(i,j,0))  /dz  + (Bz(i,j,0) - Bz(i,j-1,0)) / dy  ) ! this last term is to undo the regular full update
 
        END DO
    END DO    

    
    ! remaining half-update the normal component of E at the y- boundary
    DO k = 0, nz+3
        DO i = 0, nx+3
            Ey(i,-1,k) = Ey(i,-1,k) + 0.5 * dt * c* ( (Bx(i,-1,k) - Bx(i,-1,k-1))  /dz  + (Bz(i-1,-1,k) - Bz(i,-1,k)) / dx  )      
        END DO     
    END DO    


    

END SUBROUTINE radiation_boundary_E_y



SUBROUTINE radiation_boundary_E_z()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, cdtdz, a, b
    
    cdtdx = c*dt/dx
    cdtdy = c*dt/dy
    cdtdz = c*dt/dz
     
    !  First, half-update the normal component of E at the z- boundary
    DO j = 0, ny+3
        DO i = 0, nx+3
             Ez(i,j,-1) = Ez(i,j,-1) + 0.5 * dt * c * ( (Bx(i,j-1,-1) - Bx(i,j,-1))  /dy  + (By(i,j,-1) - By(i-1,j,-1)) / dx  )             
        END DO
    END DO    

    a = 1.d0 / (1.d0 - cdtdz)
    b = (1.d0 + cdtdz)*a
    
    ! now update surface components according to the lowest order Lindman boundary conditions 
    ! (Note: The edges aren't included and have to be updated separately.)    
    DO j = 0, ny+2
        DO i = 0, nx+2
             Ex(i,j,-1) = -Ex(i,j,0) + b * (Ex(i,j,-1)-Ex(i,j,0)) - &
                            2.d0 * cB * cdtdx * a * (Ez(i,j,-1) - Ez(i-1,j,-1)) + &
                            cE * cdtdy * a * (Bz(i,j+1,-1) + Bz(i,j+1,0) - Bz(i,j,-1) - Bz(i,j,0)) - & ! first three terms are from Lindman boundary conditions
                            dt * c* ( (By(i,j,-1) - By(i,j,0))  /dz  + (Bz(i,j,0) - Bz(i,j-1,0)) / dy  ) ! this last term is to undo the regular full update
  
             
             Ey(i,j,-1) = -Ey(i,j,0) + b * (Ey(i,j,-1)-Ey(i,j,0)) - &
                            2.d0 * cB * cdtdy * a * (Ez(i,j,-1) - Ez(i-1,j,-1)) - &
                            cE * cdtdx * a * (Bz(i+1,j,-1) + Bz(i+1,j,0) - Bz(i,j,-1) - Bz(i,j,0)) - &
                            dt * c* ( (Bx(i,j,0) - Bx(i,j,-1))  /dz  + (Bz(i-1,j,0) - Bz(i,j,0)) / dx  )  ! this last term is to undo the regular full update
 
        END DO
    END DO    

    
     ! remaining half-update the normal component of E at the z- boundary
    DO j = 0, ny+3
        DO i = 0, nx+3
             Ez(i,j,-1) = Ez(i,j,-1) + 0.5 * dt * c * ( (Bx(i,j-1,-1) - Bx(i,j,-1))  /dy  + (By(i,j,-1) - By(i-1,j,-1)) / dx  )             
        END DO
    END DO    

    

END SUBROUTINE radiation_boundary_E_z



END MODULE fieldSolver_mod