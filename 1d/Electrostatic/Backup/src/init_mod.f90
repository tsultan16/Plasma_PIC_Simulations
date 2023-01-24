MODULE init_mod

USE constants_mod
USE data_mod
USE particleMover_mod
IMPLICIT NONE

CONTAINS

! intial particle spatial distribution function (same for all species)
FUNCTION fn0(x) RESULT(fx)

	REAL*8, INTENT(IN) :: x
	INTEGER :: fx

    ! uniform spatial distribution
    fx = N/nx

END FUNCTION fn0


! initial particle velocity distribution function (same for all species)
FUNCTION fv0(v) RESULT(fx)

	REAL*8, INTENT(IN) :: v
	INTEGER :: fx

    ! uniform and zero (cold plasma)
    fx = 0.0

END FUNCTION fv0



! sets x(0) and v(0) for all particles 
SUBROUTINE init_particles()

	INTEGER ::  i, j, counter, N_cell, i_v
	REAL*8 :: wp, qm, q, m 


	! allocate memory for particles
	ALLOCATE(particles(ns,N))

    ! loop through particles and set charge and mass
    DO i = 1, ns
	    
		! calculate mass and charge of particles
		! wp^2 = n*q^2/(eps0*m) = (N/L)*q*(q/m)*(1/eps0) => q = (wp^2)*eps0*(L/N)*(m/q)
		IF(i .EQ. 1) THEN
			wp = wp_e
			qm = qm_e
		ELSE IF(i .EQ. 2) THEN
			wp = wp_i
			qm = qm_i
		END IF
		
		q = (wp**2)*eps0*(L/N)/qm
		m = q/qm
		
		PRINT*,'Species# ',i,' q,m = ',q,m
		
	    DO j = 1, N
		    particles(i,j)%q = q
		    particles(i,j)%m = m 
			
		END DO
		
 	END DO	

    ! distribute particles across grid
	counter = 1
	
    DO i = 1,nx

        ! find out how many particles to put in this cell
		N_cell = fn0(i*dx)
		
		!PRINT*,'Cell#, N_cell',i+1,N_cell
		
		IF(counter .GT. N) EXIT	
		
	    !PRINT*,'counter = ',counter
	    DO j = 1, N_cell	

            ! set initial position : x(0)		
			particles(1,counter)%x = i*dx+0.25*dx   ! grid-point centered, offset by 0.25 cell lengths to the right
			particles(2,counter)%x = i*dx+0.25*dx 
            
			! set initial velocity : v(0) sampled uniformly from (vmin,vmax)
			particles(1,counter)%vx = v0 + fv0(vmin+((vmax-vmin)/N_cell)*(j-0.5))
			particles(2,counter)%vx = v0 + fv0(vmin+((vmax-vmin)/N_cell)*(j-0.5)) 
				
			counter = counter + 1
		END DO
	
	END DO

    !PRINT*,''
	!DO i = 1, N
	!    PRINT*,'i,x(i)=',i,particles(1,i)%x 
	!    particles(i,j)%m = m
	!END DO
    !PRINT*,''

	IF(print_debug) THEN
	PRINT*,''
    DO i = 1, ns
		DO j = 1, N
			PRINT*,'PARTICLE TYPE, ID, x, vx, m, q = ',i,j,particles(i,j)%x, particles(i,j)%vx, &
			                                           particles(i,j)%m, particles(i,j)%q 
		END DO
	END DO
	PRINT*,''
	END IF

END SUBROUTINE init_particles


! accumulates charges at grid-points via interpolation
SUBROUTINE setrho()

    INTEGER :: i, j, ix, ileft, iright
	REAL*8 :: qdx, x

    ! clear density
    rho(:) = rho0
	ne(:) = 0.d0
	ni(:) = 0.d0
	
	PRINT*,''
	 
	! accumulate particle chrages at grid points
	
	! zeroth order interpolation
    IF(interpolation_type .EQ. 1) THEN
	    DO i = 1, ns
		    qdx = particles(i,1)%q/dx
			
			DO j = 1, N
			
				! apply periodic boundary condition
				IF(particles(i,j)%x .LT. xmin) THEN
					particles(i,j)%x = particles(i,j)%x + L
				END IF
				IF(particles(i,j)%x .GT. xmax) THEN
					particles(i,j)%x = particles(i,j)%x - L
				END IF

				x = particles(i,j)%x/dx 
				
				!compute index of occupied cell
				ix = x + 0.5  
                
				!PRINT*,'i,j,x,ix=',i,j,x,ix

				! accumulate charge in the cell occupied by the particle
				rho(ix) = rho(ix) + qdx
								
                IF(i .EQ. 1) ne(ix) = ne(ix) + 1.0				
				IF(i .EQ. 2) ni(ix) = ni(ix) + 1.0			
				 
		    END DO
		END DO 
	
	
	! linear interpolation
	ELSE IF(interpolation_type .EQ. 2 .OR. interpolation_type .EQ. 3) THEN

	    DO i = 1, ns
		    qdx = particles(i,1)%q/dx
			
			DO j = 1, N
				
				! apply periodic boundary condition
				IF(particles(i,j)%x .LT. xmin) THEN
					particles(i,j)%x = particles(i,j)%x + L
				END IF
				IF(particles(i,j)%x .GT. xmax) THEN
					particles(i,j)%x = particles(i,j)%x - L
				END IF
				
				x = particles(i,j)%x/dx

				!compute index of nearest left hand grid point
				ix = x 
				
				!PRINT*,'i,j,xpos,x,ix, x-ix=',i,j,particles(i,j)%x,x,ix,x-ix
	
				! accumulate charge in the cells occupied by the cloud (weighted according to the fraction of the cell occupied by the cloud)
				rho(ix)   = rho(ix) + qdx*(1.d0 - (x-ix))
				rho(ix+1) = rho(ix+1) +  qdx*(x - ix)
	           			
			    IF(i .EQ. 1) ne(ix) = ne(ix) + (1.d0 - (x-ix))				
				IF(i .EQ. 2) ni(ix) = ni(ix) + (x - ix)
				
		    END DO
		END DO 
	END IF
	
	rho(nx) = rho(nx) + rho(0)
	rho(1) = rho(1) + rho(nx+1)
	ne(nx) = ne(nx) + ne(0)
	ne(1) = ne(1) + ne(nx+1)
	ni(nx) = ni(nx) + ni(0)
	ni(1) = ni(1) + ni(nx+1)
	
    ! apply periodic boundary conditions
	rho(0) = rho(nx)
    rho(nx+1) = rho(1)	
	ne(0) = ne(nx)
    ne(nx+1) = ne(1)	
	ni(0) = ni(nx)
    ni(nx+1) = ni(1)	
	
	
	IF(print_debug) THEN
	PRINT*,''
	PRINT*,'RHO = '
	DO i = 0, nx+1
	    PRINT*,'i,rho(i)=',i,rho(i)
	END DO
	PRINT*,''
	END IF
	
END SUBROUTINE setrho

	
! sets v(0) to v(-dt/2) as required by the leap-frog integrator
SUBROUTINE setv()

	INTEGER :: i, j
	REAL*8 :: qm, wc, x, a

    ! Integrate v(0) backwards by half a time-step to obtain v(-dt/2).
    ! For a non-zero background magnetic field present, will use the Boris 1970 "half-acceleration" scheme

    DO i = 1, ns
	    IF(i .EQ. 1) THEN
			qm = qm_e
			wc = wc_e
		ELSE IF(i .EQ. 2) THEN 	
			qm = qm_i
			wc = wc_i
		END IF
		
		DO  j = 1, N
			x = particles(i,j)%x/dx 

			IF(wc .GT. 1.d-15) THEN	! w/ magnetic field 	
				! half acceleration
				CALL accel(x, qm, a)			
				particles(i,j)%vx = particles(i,j)%vx + a*(-dt/4.d0)
			
				! rotation
				particles(i,j)%vx = COS(0.5*wc*dt)*particles(i,j)%vx - SIN(0.5*wc*dt)*particles(i,j)%vy
				particles(i,j)%vy = SIN(0.5*wc*dt)*particles(i,j)%vx + COS(0.5*wc*dt)*particles(i,j)%vy
			
				! half acceleration			
				particles(i,j)%vx = particles(i,j)%vx + a*(-dt/4.d0)			

			ELSE  ! no magnetic field
				particles(i,j)%vx = particles(i,j)%vx - 0.5*a*dt
			END IF	
			
		END DO	
	END DO


END SUBROUTINE setv


! adds sinusoidal perturbations to x(0) and v(0) to electrons
SUBROUTINE add_perturbations()

	INTEGER, PARAMETER :: modes = nx
    
	REAL*8 :: x1(modes) 
	REAL*8 :: v1(modes)
	REAL*8 :: thetax(modes) 
	REAL*8 :: thetav(modes) 
	INTEGER :: i, j
    REAL*8 :: theta
	
	x1 = 0.d0
	x1(1) = 0.45d0
	!x1(17) = 0.01d0
    !x1(32) = 0.01d0
	v1 = 0.d0
	thetax = 0.d0
	thetav = 0.d0

	
	DO i = 1, N
	    DO j = 1, modes
			theta = particles(1,i)%x*j*twopi/L
			particles(1,i)%x = particles(1,i)%x + x1(j)*COS(theta+thetax(j))
			particles(1,i)%vx = particles(1,i)%vx + v1(j)*COS(theta+thetav(j))
		END DO
		
		!PRINT*,'Electron ID, x, vx = ',i,particles(1,i)%x, particles(1,i)%vx 
	END DO



END SUBROUTINE add_perturbations


END MODULE init_mod