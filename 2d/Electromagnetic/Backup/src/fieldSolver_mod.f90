MODULE fieldSolver_mod

USE constants_mod
USE data_mod
IMPLICIT NONE

CONTAINS


! This subroutine computes E.M field components by solving Maxwell's equations.
!
! Electric and Magnetic fields are decentered across the grid:
! Bz(j,k) = Bz(x_j+1/2, y_k+1/2)
! Ex(j,k) = Ex(x_j+1/2, y_k)
! Ey(j,k) = Ey(x_j, y_k+1/2)
!
!         Ey    Bz
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
    DO j = 0, nx
        DO k = 0, ny
            Bz(j,k) = Bz(j,k) - 0.5*dt*c*((Ey(j+1,k)-Ey(j,k))/dx - (Ex(j,k+1)-Ex(j,k))/dy  )
        END DO
     END DO
   
    !***************************************************************************
    ! Electric Field field update: Ex_j,k+1/2 (t) -> Ex_j,k+1/2 (t+dt)
    !                              Ey_j+1/2,k (t) -> Ey_j+1/2,k (t+dt)
    !***************************************************************************
    DO j = 0, nx
        DO k = 0, ny
            Ex(j,k) = Ex(j,k) + dt * ( c * (Bz(j,k)-Bz(j,k-1))/dy - Jx(j,k) )
            Ey(j,k) = Ey(j,k) + dt * (-c * (Bz(j,k)-Bz(j-1,k))/dx - Jy(j,k) )             
        END DO
    END DO
   
    !***************************************************************************
    ! Remaining magnetic field half-update: Bz_j+1/2,k+1/2 (t-dt/2) -> Bz_j+1/2,k+1/2 (t)  
    !***************************************************************************
    DO j = 0, nx
        DO k = 0, ny
            Bz(j,k) = Bz(j,k) - 0.5*dt*c*((Ey(j+1,k)-Ey(j,k))/dx - (Ex(j,k+1)-Ex(j,k))/dy  )
        END DO
    END DO
   
    ! apply periodic boundary condition
    Bz(-1,:) = Bz(nx,:)
    Bz(nx+1,:) = Bz(0,:) 
    Ex(-1,:) = Ex(nx,:)
    Ex(nx+1,:) = Ex(0,:)   
    Ey(-1,:) = Ey(nx,:)
    Ey(nx+1,:) = Ey(0,:) 
   
    Bz(:,-1) = Bz(:,ny)
    Bz(:,ny+1) = Bz(:,0) 
    Ex(:,-1) = Ex(:,ny)
    Ex(:,ny+1) = Ex(:,0)   
    Ey(:,-1) = Ey(:,ny)
    Ey(:,ny+1) = Ey(:,0) 
   
   
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
    DO j = 1, nx
        DO k = 1, ny
            Bz_grid(j,k) = 0.25d0*( Bz(j-1,k-1) + Bz(j,k) + Bz(j-1,k) + Bz(j,k-1) )
            Ex_grid(j,k) = 0.5d0*( Ex(j-1,k) + Ex(j,k) )
            Ey_grid(j,k) = 0.5d0*( Ey(j,k-1) + Ey(j,k) )  
        END DO
    END DO    
    ! periodic boundary
    Bz_grid(0,:) = Bz_grid(nx,:) 
    Ex_grid(0,:) = Ex_grid(nx,:)
    Ey_grid(0,:) = Ey_grid(nx,:)
    Bz_grid(nx+1,:) = Bz_grid(1,:) 
    Ex_grid(nx+1,:) = Ex_grid(1,:)
    Ey_grid(nx+1,:) = Ey_grid(1,:)
    Bz_grid(:,0) = Bz_grid(:,ny) 
    Ex_grid(:,0) = Ex_grid(:,ny)
    Ey_grid(:,0) = Ey_grid(:,ny)
    Bz_grid(:,ny+1) = Bz_grid(:,1) 
    Ex_grid(:,ny+1) = Ex_grid(:,1)
    Ey_grid(:,ny+1) = Ey_grid(:,1)    
    

END SUBROUTINE grid_avg_fields



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

	
	! apply periodic boundary conditions
	Ex_grid(nx+1,:) = Ex_grid(1,:)
    Ex_grid(0,:) = Ex_grid(nx,:)
    Ex_grid(:,ny+1) = Ex_grid(:,1)
	Ex_grid(:,0) = Ex_grid(:,ny)
	
    Ey_grid(nx+1,:) = Ey_grid(1,:)
    Ey_grid(0,:) = Ey_grid(nx,:)
    Ey_grid(:,ny+1) = Ey_grid(:,1)
	Ey_grid(:,0) = Ey_grid(:,ny)



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


! enforces field boundary conditions
SUBROUTINE field_boundary()

! periodic boundary conditions



! radiation absorbing boundary conditions


END SUBROUTINE field_boundary


END MODULE fieldSolver_mod