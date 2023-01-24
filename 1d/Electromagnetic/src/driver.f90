PROGRAM driver_plasma_1d

USE constants_mod
USE data_mod
USE init_mod
USE fieldSolver_mod
USE particleMover_mod
IMPLICIT NONE

INTEGER :: i_sim, ii
REAL*8 :: w_mode, w_exact


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

IF(ABS(wc_e) .GT. 1.d-15 .AND. interpolation_type .NE. 2) THEN
	PRINT*,'For non-zero uniform magnetic field, need to use interpolation_type = 2 (momentum conserving).'
    PRINT*,'Terminating program...'
	STOP
END IF


! open output files
OPEN(UNIT=9, FILE='particles_e.txt')
OPEN(UNIT=10, FILE='particles_i.txt')
OPEN(UNIT=11, FILE='fields.txt')
OPEN(UNIT=12, FILE='fieldsk.txt')
OPEN(UNIT=13, FILE='energy.txt')
OPEN(UNIT=14, FILE='modes.txt')
OPEN(UNIT=15, FILE='frequencies.txt' )
OPEN(UNIT=16, FILE='single_particle.txt' )
OPEN(UNIT=17, FILE='fv.txt')


! allocate memory for fields
ALLOCATE(rho(0:nx+1), phi(0:nx+1), ne(0:nx+1), ni(0:nx+1))
ALLOCATE( Ex(0:nx+1), Ey(0:nx+1), Bz(0:nx+1), Fr(0:nx+1), Fl(0:nx+1), Jy(0:nx+1), Jold(0:nx+1))
ALLOCATE(rhok(0:nx,2), phik(0:nx,2), Ksqr(0:nx,2),ESEk(0:nx),SMk(0:nx))
ALLOCATE(ESEkt(1:maxsteps,0:nx))
ALLOCATE(fv(1:ns,1:vbins))


! compute dx
dx = L/REAL(nx)

! set grid bounds
xmin = 0.5*dx
xmax = xmin+L


! set dt (needs to satrisfy courant condition dt <= dx/c)
dt = dx/c

! initialize particle distribution and state
CALL init_particles()

! initilaize electromagnetic fields (i.e. Ey and Bz)
CALL init_emfields()

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
	! move particles and deposit charges and currents to grid
    CALL move_particles()
	
	PRINT*,'Updating fields...'
	! evolve fields
	CALL compute_fields()
	
	ESEkt(i_sim,:) = ESEk(:)
	
	! compute fv
	CALL compute_velocity_distribution()
	
	! save to file
	IF(MOD(i_sim-1,tskip) .EQ. 0) CALL output()
	
	
	PRINT*,''
	PRINT*,'Total KE, ESE, Px, Py = ',KE, ESE, Px, Py
	PRINT*,''
	
    PRINT*,'Completed time-step.'

END DO

! compute mode frequencies
!w_exact = SQRT(wp_e**2+wc_e**2)
!PRINT*,''
!PRINT*,'SQRT*(wpe^2 + wce^2)=',w_exact
!PRINT*,''
!DO ii = 1, 13
!	CALL compute_frequency(ii, w_mode)
!	WRITE(15,*) ii*twopi/nx, w_mode, w_exact
!END DO

! deallocate memory for fields and particles
DEALLOCATE(rho, phi, ne, ni, particles, fv)
DEALLOCATE(Ex, Ey, Bz, Fr, Fl, Jy, Jold)
DEALLOCATE(rhok, phik, Ksqr,ESEk,SMk, ESEkt) 

PRINT*,'Simulation complete!'

CLOSE(UNIT=10)
CLOSE(UNIT=11)
CLOSE(UNIT=12)
CLOSE(UNIT=13)
CLOSE(UNIT=14)
CLOSE(UNIT=15)
CLOSE(UNIT=16)
CLOSE(UNIT=17)

CONTAINS


SUBROUTINE compute_frequency(mode, w_mode)

	INTEGER, INTENT(IN)   ::  mode
	REAL*8, INTENT(INOUT) ::  w_mode
	INTEGER :: i,j

    INTEGER :: maxpos1, maxpos2	
	REAL*8  :: maxval1, period
	LOGICAL :: maxfound
	
	! find first local max
	maxfound = .FALSE.
	maxval1= ESEkt(1,mode)	
    maxpos1 = 1	
	DO i = 2, maxSteps
	
	    IF(maxfound) EXIT
	 		
		IF(ESEkt(i,mode) .GT. ESEkt(i-1,mode) .AND. ESEkt(i,mode) .GT. ESEkt(i+1,mode) ) THEN
		    maxval1 = ESEkt(i,mode)
			maxpos1 = i 	
			maxfound = .TRUE.
		END IF			
				
	END DO

	PRINT*,'Mode#',mode,', MAXPOS1 =',maxpos1
	
    ! find second local max
	maxfound = .FALSE.
	maxval1= ESEkt(maxpos1+1,mode)		
	maxpos2 = maxpos1+1
	DO i = maxpos1+1, maxSteps

	    IF(maxfound) EXIT
	 		
		IF(ESEkt(i,mode) .GT. ESEkt(i-1,mode) .AND. ESEkt(i,mode) .GT. ESEkt(i+1,mode) ) THEN
		    maxval1 = ESEkt(i,mode)
			maxpos2 = i 	
			maxfound = .TRUE.
		END IF			
		
	END DO

	PRINT*,'Mode#',mode,', MAXPOS2 =',maxpos2

    period = 2.d0 *(maxpos2-maxpos1)*dt ! multiply by 2 because ESEk(k)~ phik*rhok ~ exp(iwt)*exp(iwt) ~ cos(2wt) 
    w_mode = twopi/period
    PRINT*,'Mode#,',mode,' w = ',w_mode
	
	

END SUBROUTINE compute_frequency


! computes space-averaged velocity distribution function
! and drift and velocity spreads for each species
SUBROUTINE compute_velocity_distribution()

	INTEGER :: i, j, i_v
    REAL*8 :: dv, vv
	
	! clear fv, vd vt
	fv = 0.d0
	vd1 = 0.d0
	vd2 = 0.d0
	vspread1 = 0.d0
	vspread2 = 0.d0
	
	dv = (vmax-vmin)/vbins 
    ! loop over all particles
    DO i = 1, ns
		DO j = 1, N
            ! compute velocity bin index
			vv = (particles(i,j)%vx-vmin)/dv
		    i_v = vv 
			
			! place particle in bin (NOTE: this is a first order interpolation)
			IF(i_v .GE. 1 .AND. i_v .LE. vbins) THEN
				fv(i,i_v)   = fv(i,i_v) + (1.d0 - (vv-i_v))
				IF(i_v .LE. vbins-1) fv(i,i_v+1) = fv(i,i_v+1) + (vv-i_v)
			END IF

            IF(i .EQ. 1) THEN
				vd1 = vd1 + particles(i,j)%vx
			ELSE IF(i .EQ. 2) THEN
				vd2 = vd2 + particles(i,j)%vx
			END IF
			
		END DO
	END DO
	
	vd1 = vd1/N
	vd2 = vd2/N
	
	! compute velocity spreads 
	DO j = 1,N
		vspread1 = vspread1 +  particles(1,j)%vx**2
		vspread2 = vspread2 +  particles(2,j)%vx**2		
	END DO
	vspread1 = vspread1/N
	vspread2 = vspread2/N

	vspread1 = vspread1-vd1**2
	vspread2 = vspread2-vd1**2
	
	
	
END SUBROUTINE compute_velocity_distribution



SUBROUTINE output()

	INTEGER :: i, j, k
    REAL*8 :: dv
	
    IF(particles_out) THEN
	! save electron states
	DO j = 1, N
		WRITE(9,*) particles(1,j)%x,particles(1,j)%y, particles(1,j)%vx,particles(1,j)%vy
		WRITE(10,*) particles(2,j)%x,particles(2,j)%y, particles(2,j)%vx,particles(2,j)%vy
	END DO
	
	WRITE(16,*) particles(1,N/2)%x,particles(1,N/2)%y, particles(2,N/2)%x,particles(2,N/2)%y

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
		IF(i_sim .GT. 0) WRITE(13,*) i_sim, KE, ESE2, KE+ESE2, Px, vd1, vspread1, vd2, vspread2 
	END IF

	IF(modes_out) THEN
	IF(i_sim .GT. 0) WRITE(14,*) i_sim, ESEk(1),ESEk(2), ESEk(4), ESEK(6),ESEk(8), ESEk(10), ESEk(12), ESEk(14) ,ESEk(16) 
	END IF
	
	IF(fv_out)THEN
		dv = (vmax-vmin)/vbins 
		DO i = 1, vbins
			WRITE(17,*) vmin+(i-0.5)*dv, fv(1,i), fv(2,i)
		END DO
	END IF
	
END SUBROUTINE output

END PROGRAM driver_plasma_1d