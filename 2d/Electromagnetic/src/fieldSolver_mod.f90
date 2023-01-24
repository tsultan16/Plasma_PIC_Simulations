MODULE fieldSolver_mod

USE constants_mod
USE data_mod
IMPLICIT NONE

REAL*8 :: bound_bufferx(-1:ny+3), bound_buffery(-1:nx+3)

CONTAINS


! This subroutine computes E.M field components by solving Maxwell's equations.
!
! Electric and Magnetic fields are decentered across the grid:
! Bz(j,k) = Bz(x_j+1/2, y_k+1/2)
! Ex(j,k) = Ex(x_j+1/2, y_k)
! Ey(j,k) = Ey(x_j, y_k+1/2)
!
!         Ey    Bz-
!    ----_*-----*(j+1/2,k+1/2) 
!   | (j,k+1/2) |
!   |           |
!   |     *     * Ex
!   |   (j,k)   |(j+1/2,k)
!   |           |
!    -----------
!
! Solving Maxwell's equations for Ex, Ey and Bz amounts to the following operation: 

! *********************************************************************************************************************
! * change in E (B) flux through cell-face = line integral of B (E) across boundary of that face + charge-flux (zero) *
! *********************************************************************************************************************
!
! Flux of Bz across z face:
! dx*dy*[(Bz(new)-Bz(old))_(j+1/2,k+1/2)] /dt = -dx*c*[ Ex_(j+1/2,k) -  Ex_(j+1/2,k+1) ] - dy*c*[ Ey_(j+1,k+1/2) -  Ey_(j,k+1/2) ]   
! 
! Flux of Ex across x-face:
! dy*dz*[(Ex(new)-Ex(old))_(j+1/2,k)] /dt = dz*c*[ Bz_(j+1/2,k+1/2) -  Bz_(j+1/2,k-1/2) ] - dy*dz*Jx_(j+1/2,k) + line integral of By along y direction off-plane (does not exist in 2D)  
!
! Flux of Ey across y-face:
! dx*dz*[(Ey(new)-Ey(old))_(j,k+1/2)] /dt = dz*c*[ Bz_(j-1/2,k+1/2) -  Bz_(j+1/2,k+1/2) ] - dx*dz*Jy_(j,k+1/2) + line integral of Bx along x direction off-plane (does not exist in 2D)  
!
!Integration of Ex, Ey fluxes leap-frogs over integration of Bz flux due to required time-centered values on the line integration on the R.H.S.

SUBROUTINE compute_fields()

	INTEGER :: kx, ky, i, j, k
   
    !***************************************************************************
    ! Magnetic field half-update: Bz_j+1/2,k+1/2 (t-dt/2) -> Bz_j+1/2,k+1/2 (t)  
    !***************************************************************************
    IF(bndry .EQ. 2) THEN
        ! store values of last interior layer, required later for boundary update
        bound_bufferx(-1:ny+3) = Bz(nx+2,-1:ny+3)   
        bound_buffery(-1:nx+3) = Bz(-1:nx+3,ny+2)
    END IF
    
    DO j = -1, nx+2
        DO k = -1, ny+2
            Bz(j,k) = Bz(j,k) - 0.5*dt*c*((Ey(j+1,k)-Ey(j,k))/dx - (Ex(j,k+1)-Ex(j,k))/dy  )
        END DO
     END DO
   
   
    ! set boundary magnetic field 
    ! Bz on x+ (j=nx+3) and y+ edges (k=ny+3)
    IF(bndry .EQ. 1) THEN
       CALL periodic_boundary_B()
    ELSE IF(bndry .EQ. 2) THEN
       CALL radiation_boundary_B()
    END IF
    
   
    !***************************************************************************
    ! Electric Field field update: Ex_j,k+1/2 (t) -> Ex_j,k+1/2 (t+dt)
    !                              Ey_j+1/2,k (t) -> Ey_j+1/2,k (t+dt)
    !***************************************************************************
    DO j = -1, nx+3
        DO k = 0, ny+3
            Ex(j,k) = Ex(j,k) + dt * ( c * (Bz(j,k)-Bz(j,k-1))/dy - Jx(j,k) )
        END DO
    END DO
   
    DO j = 0, nx+3
        DO k = -1, ny+3
            Ey(j,k) = Ey(j,k) + dt * (-c * (Bz(j,k)-Bz(j-1,k))/dx - Jy(j,k) )             
        END DO
    END DO

    
    ! set boundary electric field 
    ! Ey on x- edge (j=-1) and Ex on y- edge (k=-1)
    IF(bndry .EQ. 1) THEN
       CALL periodic_boundary_E()
    ELSE IF(bndry .EQ. 2) THEN
       CALL radiation_boundary_E()
    END IF
    
    
    
   
    !***************************************************************************
    ! Remaining magnetic field half-update: Bz_j+1/2,k+1/2 (t-dt/2) -> Bz_j+1/2,k+1/2 (t)  
    !***************************************************************************
    DO j = -1, nx+2
        DO k = -1, ny+2
            Bz(j,k) = Bz(j,k) - 0.5*dt*c*((Ey(j+1,k)-Ey(j,k))/dx - (Ex(j,k+1)-Ex(j,k))/dy  )
        END DO
    END DO
      
    ! set boundary magnetic field
    ! ......

    
  
   
    ! compute grid average fields (will be used in particle mover during force interpolation)
    CALL grid_avg_fields()
   
   
   
	IF(print_debug) THEN
	PRINT*,''
	PRINT*,'Ex = '
	DO j = ny+1, 0, -1
        DO i = 0, nx+1
            WRITE(*,FMT='(1f5.1)',ADVANCE='NO') Ex(i,j)
        END DO
        PRINT*,''
    END DO

    PRINT*,''
	PRINT*,'Ey = '
	DO j = ny+1, 0, -1
        DO i = 0, nx+1
            WRITE(*,FMT='(1f5.1)',ADVANCE='NO') Ey(i,j)
        END DO
        PRINT*,''
    END DO
    
	END IF

END SUBROUTINE compute_fields


SUBROUTINE grid_avg_fields()

   INTEGER :: j, k
   
   ! grid-point averages of em fields, i.e. Ex_j,k, Ey_j,k, Bz_j,k
    DO j = 0, nx+1
        DO k = 0, ny+1
            Bz_grid(j,k) = 0.25d0*( Bz(j-1,k-1) + Bz(j,k) + Bz(j-1,k) + Bz(j,k-1) )
            Ex_grid(j,k) = 0.5d0*( Ex(j-1,k) + Ex(j,k) )
            Ey_grid(j,k) = 0.5d0*( Ey(j,k-1) + Ey(j,k) )  
        END DO
    END DO    
    
  
END SUBROUTINE grid_avg_fields


!*****************************
! Field Boundary Conditions
!*****************************

! Periodic boundaries
SUBROUTINE periodic_boundary_B()

    Bz(nx+3,:) = Bz(-1,:) 
    Bz(:,ny+3) = Bz(:,-1) 


END SUBROUTINE periodic_boundary_B


SUBROUTINE periodic_boundary_E()

    Ex(:,-1) = Ex(:,ny+3)
    Ey(-1,:) = Ex(nx+3,:)

END SUBROUTINE periodic_boundary_E



! Radiation absorbing boundaries

!At x+ boundary, outgoing Bz plane wave satisfies: (1/c d/dt + d/dx) Bz = c_E (d/dy) Ex  
!At y+ boundary, "         "        "     "      : (1/c d/dt + d/dy) Bz = -c_E (d/dx) Ey
! where c_E = 1-c_B = 0.5*(1+tan^2(theta)) ~ 0.5858 (with theta = 45 degrees, which seems to be a good approximation for arbitrary theta)
!
! finite difference version of the above equatons can be solved to obtain the required boundary values

SUBROUTINE radiation_boundary_B()

    INTEGER :: i,j,k
    REAL*8 :: cdtdx, cdtdy, a, b
    
    cdtdx = 0.5d0*c*dt/dx
    cdtdx = 0.5d0*c*dt/dy
    a = (1.d0-cdtdx)/(1.d0+cdtdx)
    b = cE*c*dt/(dy*(1.d0+cdtdx))
     
    !  x+ boundary update (top-right corner value not updated)
    DO k = -1, ny+2
        Bz(nx+3,k) = bound_bufferx(k) + a*(Bz(nx+3,k)-Bz(nx+2,k)) + b*(Ex(nx+3,k+1)-Ex(nx+3,k))
    END DO    

    a = (1.d0-cdtdy)/(1.d0+cdtdy)
    b = cE*c*dt/(dx*(1.d0+cdtdy))
     
    !  y+ boundary update (top-right corner value not updated)
    DO j = -1, nx+2
        Bz(j,ny+3) = bound_buffery(j) + a*(Bz(j,ny+3)-Bz(j,ny+2)) - b*(Ey(j+1,ny+3)-Ey(j,ny+3))
    END DO    


END SUBROUTINE radiation_boundary_B



!At x- boundary, outgoing Ey plane wave satisfies: (1/c d/dt - d/dx) Ey = (c_B/(1-2 c_B)) (d/dy) Ex
!At y- boundary, "        Ex        "     "      : (1/c d/dt - d/dy) Ex = (c_E/(1-2 c_E)) (d/dx) Ey
! where c_E = 1-c_B = 0.5*(1+tan^2(theta)) ~ 0.5858 

SUBROUTINE radiation_boundary_E()

    ! 
    
END SUBROUTINE radiation_boundary_E


!*************************
! POISSON SOLVER
!*************************

! can use this subroutine to set up initial Electrostatic field
SUBROUTINE poisson_solve()

    INTEGER :: kx, ky, i, j, k

    
    ! compute DFT of rho
	CALL slow_DFT(rho, rhok)

	! set k=0 mode to zero (for periodic boundaries, average charge density must be zero)
	rhok(0,0,:) = 0.d0


	! compute K squared
	DO kx = -nx/2, nx/2-1
        DO ky = -ny/2, ny/2-1
            !attentuation/compensation factor
            SMk(kx,ky) = EXP(a1*SIN(0.5*twopi*kx*dx/Lx)**2-a2*TAN(0.5*twopi*kx*dx/Lx)**2)
            SMk(kx,ky) = SMk(kx,ky)*EXP(a1*SIN(0.5*twopi*ky*dy/Ly)**2-a2*TAN(0.5*twopi*ky*dy/Ly)**2)
            ! real part
            Ksqr(kx,ky,1) = ((2.d0*SIN(0.5*twopi*kx*dx/Lx)/dx)**2 + &
                             (2.d0*SIN(0.5*twopi*ky*dy/Ly)/dy)**2)
            Ksqr(kx,ky,1) = Ksqr(kx,ky,1)/(Smk(kx,ky)**2)
            ! imaginary part
            Ksqr(kx,ky,2) = 0.d0
        END DO        
	END DO
   	
	! compute phik
    DO kx = -nx/2, nx/2-1
        DO ky = -ny/2, ny/2-1
            phik(kx,ky,1) = rhok(kx,ky,1)/(eps0*Ksqr(kx,ky,1))
            phik(kx,ky,2) = rhok(kx,ky,2)/(eps0*Ksqr(kx,ky,1))    
        END DO        
	END DO
    phik(0,0,:) = 0.d0 ! get rid of singular value

    ! inverse transform phik
    CALL slow_DFTI(phik,phi)
	
	! apply periodic boundary conditions
	phi(nx+1,:) = phi(1,:)
	phi(0,:) = phi(nx,:)
    phi(:,ny+1) = phi(:,1)
	phi(:,0) = phi(:,ny)
    
	
	! compute electric field
	DO j = 1, nx
        DO k = 1, ny
            Ex_grid(j,k) = -0.5d0*(phi(j+1,k)-phi(j-1,k))/dx
            Ey_grid(j,k) = -0.5d0*(phi(j,k+1)-phi(j,k-1))/dy
        END DO        
	END DO

	
	! apply boundary conditions
    IF(bndry .EQ. 1) THEN
        Ex_grid(nx+1,:) = Ex_grid(1,:)
        Ex_grid(0,:) = Ex_grid(nx,:)
        Ex_grid(:,ny+1) = Ex_grid(:,1)
        Ex_grid(:,0) = Ex_grid(:,ny)
	
        Ey_grid(nx+1,:) = Ey_grid(1,:)
        Ey_grid(0,:) = Ey_grid(nx,:)
        Ey_grid(:,ny+1) = Ey_grid(:,1)
        Ey_grid(:,0) = Ey_grid(:,ny)
    END IF


    ! Interpolate to cell interfaces
	DO j = 1, nx
        DO k = 1, ny
            Ex(j,k) = 0.5d0*(Ex_grid(j+1,k)+Ex_grid(j,k))
            Ey(j,k) = 0.5d0*(Ey_grid(j+1,k)+Ey_grid(j,k))
        END DO        
	END DO           


END SUBROUTINE poisson_solve



! computes the discrete fourier transform of f and stores result in fk (both real and imaginary parts)
SUBROUTINE slow_DFT(f, fk)

	REAL*8, INTENT(IN) :: f(0:nx+1, 0:ny+1)
	REAL*8, INTENT(INOUT) :: fk(-nx/2:nx/2 - 1,-ny/2:ny/2-1 ,2) 

	INTEGER :: kx, ky, j, k
	
	fk = 0.d0
	
	
	DO kx = -nx/2, nx/2 -1
        DO ky = -ny/2, ny/2 -1
            DO j = 1, nx
                DO k = 1, ny
                    ! real part
                    fk(kx,ky,1) = fk(kx,ky,1) + f(j,k)*COS(twopi*( kx*j*dx/Lx + ky*k*dy/Ly )) 
                    ! imaginary part
                    fk(kx,ky,2) = fk(kx,ky,2) - f(j,k)*SIN(twopi*( kx*j*dx/Lx + ky*k*dy/Ly ))
                END DO    
            END DO
            
		END DO
	END DO
	
	
    fk = fk*dx*dy

	
	
END SUBROUTINE slow_DFT



! computes the inverse discrete fourier transform of fk and stores result (real part) in f
SUBROUTINE slow_DFTI(fk, f)

	REAL*8, INTENT(IN) :: fk(-nx/2:nx/2 - 1,-ny/2:ny/2-1 ,2)  
	REAL*8, INTENT(INOUT) :: f(0:nx+1,0:ny+1)

	INTEGER :: kx, ky, j, k
	
	f = 0.d0

    DO j = 1, nx
        DO k = 1, ny
            DO kx = -nx/2, nx/2 -1
                DO ky = -ny/2, ny/2 -1 
                    f(j,k) = f(j,k) + fk(kx,ky,1)*COS(twopi*( kx*j*dx/Lx + ky*k*dy/Ly )) - &
                                      fk(kx,ky,2)*SIN(twopi*( kx*j*dx/Lx + ky*k*dy/Ly ))
                END DO    
            END DO
        END DO    
    END DO
            

    f = f/(Lx*Ly)

END SUBROUTINE slow_DFTI



END MODULE fieldSolver_mod