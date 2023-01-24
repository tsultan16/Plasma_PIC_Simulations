MODULE fieldSolver_mod

USE constants_mod
USE data_mod
IMPLICIT NONE

CONTAINS

SUBROUTINE compute_fields()

	INTEGER :: i,k


    ESE2 = ESE

    ! compute DFT of rho
	CALL slow_DFT(rho, rhok)

	! set k=0 mode to zero (for periodic boundaries, average charge density must be zero)
	rhok(0,:) = 0.d0

    !PRINT*,''
	!PRINT*,'RHOK ='
	!DO k = 1, nx/2 -1
	!    PRINT*,'k,Re[rhok(k)],Im[rhok(k)]=',k,rhok(k,1),rhok(k,2)
	!END DO
	!PRINT*,''

	! compute K squared
	DO k = 1, nx-1
	
	    !attentuation/compensation factor
	    SMk(k)=EXP(a1*SIN(0.5*twopi*k*dx/L)**2-a2*TAN(0.5*twopi*k*dx/L)**2)
	    ! real part
		Ksqr(k,1) = (2.d0*SIN(0.5*twopi*k*dx/L)/dx)**2/(Smk(k)**2)
		! imaginary part
		Ksqr(k,2) = 0.d0
	END DO
   

	!PRINT*,''
	!PRINT*,'Ksqr ='
	!DO k = -nx/2, nx/2 -1
	!    PRINT*,'k,Re[Ksqr(k)]=',k,Ksqr(k,1)
	!END DO
	!PRINT*,''

	
	! compute phik
	DO k = 1, nx-1
		phik(k,1) = rhok(k,1)/(eps0*Ksqr(k,1))
		phik(k,2) = rhok(k,2)/(eps0*Ksqr(k,1))
	END DO
	phik(0,:) = 0 ! get rid of singular value
	
	!PRINT*,''
	!PRINT*,'PHIK ='
	!DO k = 1, nx/2 -1
	!    PRINT*,'k,Re[phik(k)],Im[phik(k)]=',k,phik(k,1),phik(k,2)
	!END DO
	!PRINT*,''

    ! inverse transform phik
    CALL slow_DFTI(phik,phi)
	
	! apply periodic boundary conditions
	phi(nx+1) = phi(1)
	phi(0) = phi(nx)
	
	IF(print_debug) THEN
	PRINT*,''
	DO i = 0, nx+1
		PRINT*,'i,phi(i)=',i,phi(i)
	END DO
	PRINT*,''
	END IF
	
	! compute electric field
	IF(interpolation_type .EQ. 3) THEN  ! forward difference

		DO i = 1,nx
			E(i) = -(phi(i+1)-phi(i))/dx
		END DO
		
	ELSE  ! centered difference

		DO i = 1,nx
			E(i) = -0.5*(phi(i+1)-phi(i-1))/dx
		END DO

	END IF
	
	! apply periodic boundary conditions
	E(nx+1) = E(1)
	E(0) = E(nx)
	
	! compute Electrostatic Field energy
	ESE = 0.d0
	DO k = 1, nx-1
	    ESEk(k) = 0.5*(rhok(k,1)*phik(k,1)+rhok(k,2)*phik(k,2))
		ESE = ESE + (rhok(k,1)*phik(k,1)+rhok(k,2)*phik(k,2))
	END DO
	ESE = 0.5*ESE/L
	
	
	IF(print_debug) THEN
	PRINT*,''
	DO i = 1, nx
		PRINT*,'i,E(i)=',i,E(i)
	END DO
	PRINT*,''
	END IF

END SUBROUTINE compute_fields



! computes the discrete fourier transform of f and stores result in fk
SUBROUTINE slow_DFT(f, fk)

	REAL*8, INTENT(IN) :: f(0:nx+1)
	REAL*8, INTENT(INOUT) :: fk(0:nx,2) 

	INTEGER :: k, j
	
	fk = 0.d0
	
    !DO k = -nx/2, nx/2 -1
	!	DO j = 0, nx-1
	!	    ! real part
	!		fk(k,1) = fk(k,1) + f(j)*COS(twopi*k*j*dx/L) 
	!		! imaginary part
	!		fk(k,2) = fk(k,2) - f(j)*SIN(twopi*k*j*dx/L)
	!	END DO
	!END DO
	
	DO k = 1, nx-1
		DO j = 1, nx
		    ! real part
			fk(k,1) = fk(k,1) + f(j)*COS(twopi*k*j*dx/L) 
			! imaginary part
			fk(k,2) = fk(k,2) - f(j)*SIN(twopi*k*j*dx/L)
		END DO
	END DO
	
	
    fk = fk*dx

	
	
END SUBROUTINE slow_DFT



! computes the inverse discrete fourier transform of fk and stores result (real part) in f
SUBROUTINE slow_DFTI(fk, f)

	REAL*8, INTENT(IN) :: fk(0:nx,2) 
	REAL*8, INTENT(INOUT) :: f(0:nx+1)

	INTEGER :: k, j
	
	f = 0.d0


	DO j = 1, nx
		DO k = 1, nx-1
	        f(j) = f(j) + fk(k,1)*COS(twopi*k*j*dx/L) - fk(k,2)*SIN(twopi*k*j*dx/L)     
		END DO
	END DO

    f = f/L

END SUBROUTINE slow_DFTI






END MODULE fieldSolver_mod