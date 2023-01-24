PROGRAM driver_plasma_1d

USE constants_mod
USE data_mod
USE init_mod
USE fieldSolver_mod
USE particleMover_mod
IMPLICIT NONE


INTEGER :: i_sim, ii
REAL*8 :: w_mode, w_exact

OPEN(UNIT=9, FILE='particles_e.txt')
OPEN(UNIT=10, FILE='particles_i.txt')
OPEN(UNIT=11, FILE='fields.txt')
OPEN(UNIT=12, FILE='fieldsk.txt')
OPEN(UNIT=13, FILE='energy.txt')
OPEN(UNIT=14, FILE='modes.txt')
OPEN(UNIT=15, FILE='frequencies.txt' )

! check that there are more particles than cells
IF(N .LT. nx) THEN
	PRINT*,'Need N > nx.'
	PRINT*,'Terminating program..'
	STOP
END IF

IF(MOD(N,nx) .NE. 0) THEN
	PRINT*,'Need N to be a multiple of nx.'
	PRINT*,'Terminating program..'
	STOP
END If

! compute dx
dx = L/REAL(nx)

! set grid bounds
xmin = 0.5*dx
xmax = xmin+L


! set dt (arbitrary for now)
dt = 0.05*twopi/wp_e

! initialize particle distribution and state
CALL init_particles()

! allocate memory for fields
ALLOCATE(rho(0:nx+1), E(0:nx+1), phi(0:nx+1), ne(0:nx+1), ni(0:nx+1))
ALLOCATE(rhok(0:nx,2), phik(0:nx,2), Ksqr(0:nx,2),ESEk(0:nx),SMk(0:nx))
ALLOCATE(ESEkt(1:maxsteps,0:nx))

! add sinusoidal perturbation to electrons
CALL add_perturbations()

! compute initial charge distribution
CALL setrho()

! compute initial electric field
CALL compute_fields()

! set v(-dt/2)
CALL setv()

CALL output()

PRINT*,'xmin,xmax,dx=',xmin,xmax,dx

! enter simulation loop
DO i_sim = 1, maxsteps
	
	PRINT*,''
	PRINT*,'Time step, dt, % complete = ',i_sim,dt,100.*REAL(i_sim)/REAL(maxsteps)
	PRINT*,''
	
	PRINT*,'Moving particles...'
	! move particles
    CALL move_particles()
	
	PRINT*,'Updating fields...'
	! evolve fields
	CALL compute_fields()
	
	ESEkt(i_sim,:) = ESEk(:)
	
	CALL output()
	
	
	PRINT*,''
	PRINT*,'Total KE, ESE, Px, Py = ',KE, ESE, Px, Py
	PRINT*,''
	
    PRINT*,'Completed time-step.'

END DO

! compute mode frequencies
w_exact = SQRT(wp_e**2+wp_i**2)
PRINT*,''
PRINT*,'SQRT*(wpe^2 + wpi^2)=',w_exact
PRINT*,''
DO ii = 1, 13
	CALL compute_frequency(ii, w_mode)
	WRITE(15,*) ii*twopi/nx, w_mode, w_exact
END DO

! deallocate memory for fields and particles
DEALLOCATE(rho, E, phi, ne, ni, particles)
DEALLOCATE(rhok, phik, Ksqr,ESEk,SMk, ESEkt) 

PRINT*,'Simulation complete!'

CLOSE(UNIT=10)
CLOSE(UNIT=11)
CLOSE(UNIT=12)
CLOSE(UNIT=13)
CLOSE(UNIT=14)
CLOSE(UNIT=15)


CONTAINS


SUBROUTINE compute_frequency(mode, w_mode)

	INTEGER, INTENT(IN)   ::  mode
	REAL*8, INTENT(INOUT) ::  w_mode
	INTEGER :: i,j

    INTEGER :: minpos1, minpos2	
	REAL*8  :: minval1, period
	LOGICAL :: minfound
	
	! find first local min
	minfound = .FALSE.
	minval1= ESEkt(1,mode)		
	DO i = 2, maxSteps
	
	    IF(minfound) EXIT
	 		
		IF(minval1 .GT. ESEkt(i,mode)) THEN
		    minval1 = ESEkt(i,mode)
			minpos1 = i 
			
			IF(ESEkt(i+1,mode) .GT. ESEkt(i,mode)) minfound = .TRUE.

		END IF			
				
	END DO

	
    ! find second local min
	minfound = .FALSE.
	minval1= ESEkt(minpos1+1,mode)		
	DO i = minpos1+1, maxSteps

	    IF(minfound) EXIT
	 		
		IF(minval1 .GT. ESEkt(i,mode)) THEN
		    minval1 = ESEkt(i,mode)
			minpos2 = i 

			IF(ESEkt(i+1,mode) .GT. ESEkt(i,mode)) minfound = .TRUE.

		END IF			
		
		
	END DO

    period = 2.d0 *(minpos2-minpos1)*dt ! multiply by 2 because ESEk(k)~ phik*rhok ~ exp(iwt)*exp(iwt) ~ cos(2wt) 
    w_mode = twopi/period
    PRINT*,'Mode#,',mode,' w = ',w_mode
	
	

END SUBROUTINE compute_frequency


SUBROUTINE output()

	INTEGER :: i, j, k

    IF(particles_out) THEN
	! save electron states
	DO j = 1, N
		WRITE(9,*) particles(1,j)%x,0.5, particles(1,j)%vx,particles(1,j)%vy
		WRITE(10,*) particles(2,j)%x,0.5, particles(2,j)%vx,particles(2,j)%vy
	END DO
	END IF
	
	
	IF(fields_out) THEN
    ! save rho, phi and E
    DO i = 1, nx
		WRITE(11,*) (i+0.5)*dx,rho(i),phi(i),E(i), ne(i),ni(i) 
	END DO
	END IF	
	
	IF(fieldsk_out) THEN
	! save Re[rhok] and Re[phik]
    DO k = 1, nx-1
		WRITE(12,*) k*twopi/L, rhok(k,1),phik(k,1),ESEk(k) 
	END DO
	END IF
	
	IF(energy_out) THEN
		IF(i_sim .GT. 0) WRITE(13,*) i_sim, KE, ESE2, KE+ESE2, Px 
	END IF

	IF(modes_out) THEN
	IF(i_sim .GT. 0) WRITE(14,*) i_sim, ESEk(1),ESEk(2), ESEk(4), ESEK(6),ESEk(8), ESEk(10), ESEk(12), ESEk(14) ,ESEk(16) 
	END IF
	
END SUBROUTINE output

END PROGRAM driver_plasma_1d