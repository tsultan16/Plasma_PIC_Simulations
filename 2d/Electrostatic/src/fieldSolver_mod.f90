MODULE fieldSolver_mod

USE constants_mod
USE data_mod
IMPLICIT NONE

CONTAINS

SUBROUTINE compute_fields()

	INTEGER :: kx, ky, i, j, k

    ESE2 = ESE

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
    
	
	IF(print_debug) THEN
	PRINT*,''
	PRINT*,'PHI = '
	DO j = ny+1, 0, -1
        DO i = 0, nx+1
            WRITE(*,FMT='(1f5.1)',ADVANCE='NO') phi(i,j)
        END DO
        PRINT*,''
    END DO    
	END IF
	
	! compute electric field
	IF(interpolation_type .EQ. 3) THEN  ! forward difference

		DO j = 1, nx
            DO k = 1, ny
                Ex(j,k) = -(phi(j+1,k)-phi(j,k))/dx
                Ey(j,k) = -(phi(j,k+1)-phi(j,k))/dy
            END DO        
		END DO
		
	ELSE  ! centered difference

		DO j = 1, nx
            DO k = 1, ny
                Ex(j,k) = -0.5d0*(phi(j+1,k)-phi(j-1,k))/dx
                Ey(j,k) = -0.5d0*(phi(j,k+1)-phi(j,k-1))/dy
            END DO        
		END DO
		
	END IF
	
	! apply periodic boundary conditions
	Ex(nx+1,:) = Ex(1,:)
    Ex(0,:) = Ex(nx,:)
    Ex(:,ny+1) = Ex(:,1)
	Ex(:,0) = Ex(:,ny)
	
    Ey(nx+1,:) = Ey(1,:)
    Ey(0,:) = Ey(nx,:)
    Ey(:,ny+1) = Ey(:,1)
	Ey(:,0) = Ey(:,ny)
	
    
	! compute Electrostatic Field energy
	ESE = 0.d0
	DO kx = -nx/2, nx/2-1
        DO ky = -ny/2, ny/2-1
            ESEk(kx,ky) = 0.5*(rhok(kx,ky,1)*phik(kx,ky,1)+rhok(kx,ky,2)*phik(kx,ky,2))
            ESE = ESE + (rhok(kx,ky,1)*phik(kx,ky,1)+rhok(kx,ky,2)*phik(kx,ky,2))
        END DO
    END DO
	ESE = 0.5*ESE/(Lx*Ly)
	
	
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



! computes the discrete fourier transform of f and stores result in fk
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