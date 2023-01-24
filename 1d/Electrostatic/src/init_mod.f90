MODULE init_mod

USE constants_mod
USE data_mod
USE particleMover_mod
IMPLICIT NONE

CONTAINS

! intial particle spatial distribution function (same for all species)
FUNCTION fn0(x) RESULT(fx)

	REAL*8, INTENT(IN) :: x
	REAL*8 :: fx

    ! uniform spatial distribution
    fx = 1/L!N/nx

END FUNCTION fn0


! initial particle velocity distribution function (same for all species)
FUNCTION fv0(v) RESULT(fv)

	REAL*8, INTENT(IN) :: v
	REAL*8 :: fv

    ! uniform and zero (cold plasma)
    !fv = 0.0

    ! Maxwellian
    fv = EXP(-(v/vt)**2)/SQRT(0.5*twopi)/vt

END FUNCTION fv0



! sets x(0) and v(0) for all particles 
SUBROUTINE init_particles()

	INTEGER ::  i, j, counter, j1, N_cell, i_v, ngr, nv, ii
	REAL*8 :: wp, qm, q, m, xx, dxx, vv, dvv, Dv_new, Dv_old, Rs, dRs, temp, p
    REAL*8, ALLOCATABLE :: temp_pos(:) 

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
		
		rho0 =-q*(N/nx)/dx  ! neutralizing background charge
		
	    DO j = 1, N
		    particles(i,j)%q = q
		    particles(i,j)%m = m 
			!IF(i .EQ. 1) particles(i,j)%vx = v0
            !IF(i .EQ. 2) particles(i,j)%vx = -v0
            
		END DO
		
 	END DO	


    !**************************************************************************
    !  Load particles in velocity space
    !**************************************************************************
    
    ! distriubute particles across velocity space (divide into equal sized 
    ! "loading groups", load for one group, then copy into the remaining groups)
    
    ngr = N/nlg   ! number of particles in loading group

    !first compute cumulative velocity distribution (mid-point rule integration)  
    nv = 32*N
    ALLOCATE(CD(-nv:nv))
    CD = 0.d0
    
    dvv = 2*vmax/(2*nv)
    CD(-nv) = 0.d0
    DO i = -nv+1, nv
        vv = -vmax + (i+nv-0.5)*dvv
        CD(i) = CD(i-1) + MAX(fv0(vv),0.d0)*dvv                       
    END DO

    dRs = CD(nv)/(ngr+1)  ! spacing between Rs (ngr equally spaced numbers between 0 and 1)
    Rs = dRs
    counter = 1
    DO i = -nv+1, nv
 
        vv = -vmax + (i+nv-0.5)*dvv 

        IF(Rs .GT. CD(i-1) .AND. Rs .LE. CD(i)) THEN
            !PRINT*,'Rs, vv, counter =', Rs, vv, counter
            
            ! assign this velocity to particle (Note: this is a zeroth order assignment. For a smoother distribution, can use linear weighting)         
            particles(1,counter)%vx = particles(1,counter)%vx + vv
            !particles(2,counter)%vx = particles(2,counter)%vx + vv
            
            counter = counter + 1            
            Rs = Rs + dRs            
        END IF
        
        IF(counter .GT. ngr) EXIT
        
    END DO
    
    DEALLOCATE(CD)
    
    ! copy this distribution into remaining loading groups (if any)
    IF(nlg .GT. 1) THEN
        DO i = 2,nlg
            j1 = ngr*(i-1)  
            DO j = 1, ngr
                particles(1,j1+j)%vx = particles(1,j1+j)%vx  + particles(1,j)%vx
                !particles(2,j1+j)%vx = particles(2,j1+j)%vx  + particles(2,j)%vx
            END DO
        END DO
    END IF
       

    !**************************************************************************
    !  Load particles in position space
    !**************************************************************************


    ! distriubute particles across the grid (divide into equal sized 
    ! "loading groups", load for one group, then copy into the remaining groups)
    
    !first compute cumulative spatial distribution (mid-point rule integration)  
    !nv = 16*N
    ALLOCATE(CD(-nv:nv))
    CD = 0.d0

    dxx = L/(2*nv)
    CD(-nv) = 0.d0
    DO i = -nv+1, nv
        xx = (i+nv-1.5d0)*dxx
        CD(i) = CD(i-1) + MAX(fn0(xx),0.d0)*dxx                     
    END DO

    dRs = CD(nv)/(ngr)  ! spacing between Rs (ngr equally spaced numbers between 0 and 1)
    Rs = dRs
    counter = 1
    DO i = -nv+1, nv
 
        xx = (i+nv-1.5d0)*dxx 
        IF(Rs .GT. CD(i-1) .AND. Rs .LE. CD(i)) THEN
            !PRINT*,'Rs, xx, counter =', Rs, xx, counter
            
            ! assign this position to particle (Note: this is a zeroth order assignment. For a smoother distribution, can use linear weighting)         
            particles(1,counter)%x = xx
            !particles(2,counter)%x = xx
            
            counter = counter + 1            
            Rs = Rs + dRs            
        END IF
        
        IF(counter .GT. ngr) EXIT
        
    END DO
    
    DEALLOCATE(CD)
    
    ! copy this distribution into remaining loading groups (if any)
    IF(nlg .GT. 1) THEN
        DO i = 2,nlg
            j1 = ngr*(i-1)  
            DO j = 1, ngr
                particles(1,j1+j)%x = particles(1,j)%x
                !particles(2,j1+j)%x = particles(2,j)%x
            END DO
        END DO
    END IF
       
    ! scrable the position to fill out phase space more uniformly (Fisher-Yates shuffle)     
    DO i = 1, ns 
        
        DO j = N, 1, -1
            !generate a random integer between 1 and j-1
            CALL RANDOM_NUMBER(p)
            ii = 1 + p*(j-1)
        
            ! swap positions at indices j and ii
            temp = particles(i,j)%x
            particles(i,j)%x = particles(i,ii)%x        
            particles(i,ii)%x = temp 
             
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
			PRINT*,'PARTICLE TYPE, ID, x, vx = ',i,j,particles(i,j)%x, particles(i,j)%vx
		END DO
	END DO
	PRINT*,''
	END IF

END SUBROUTINE init_particles



! accumulates charges at grid-points via interpolation
SUBROUTINE setrho()

    INTEGER :: i, j, ix, ileft, iright
	REAL*8 :: qdx, x, qadd

    ! clear density
    rho(:) = rho0
	ne(:) = 0.d0
	ni(:) = 0.d0
	rho(0) = 0.d0
	rho(nx+1) = 0.d0
	
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
                
				! accumulate charge in the cell occupied by the particle
				rho(ix) = rho(ix) + qdx
				
				!PRINT*,'i,j,x,ix,qadd=',i,j,x,ix,qdx

				
                IF(i .EQ. 1) ne(ix) = ne(ix) + 1.0				
				!IF(i .EQ. 2) ni(ix) = ni(ix) + 1.0			
				 
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
					
				! accumulate charge in the cells occupied by the cloud (weighted according to the fraction of the cell occupied by the cloud)
				rho(ix)   = rho(ix) + qdx*(1.d0 - (x-ix))
				rho(ix+1) = rho(ix+1) +  qdx*(x - ix)
	           	
				!PRINT*,'i,j,ix,ix+1,qadd(ix),qadd(ix+1)=',i,j,ix,ix+1,(1.d0 - (x-ix)),(x - ix)

                
			    IF(i .EQ. 1) ne(ix) = ne(ix) + (1.d0 - (x-ix))				
				!IF(i .EQ. 2) ni(ix) = ni(ix) + (x - ix)
				
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
	REAL*8 :: qm, wc, x, a, vxx

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

			IF(ABS(wc) .GT. 1.d-15) THEN	! w/ magnetic field 	
				! half acceleration
				CALL accel(x, qm, a)			
				particles(i,j)%vx = particles(i,j)%vx + a*(-dt/4.d0)
			
			    vxx = particles(i,j)%vx
				! rotation
				particles(i,j)%vx = COS(-0.5*wc*dt)*particles(i,j)%vx - SIN(-0.5*wc*dt)*particles(i,j)%vy
				particles(i,j)%vy = SIN(-0.5*wc*dt)*vxx + COS(-0.5*wc*dt)*particles(i,j)%vy
			
				! half acceleration			
				particles(i,j)%vx = particles(i,j)%vx + a*(-dt/4.d0)			

			ELSE  ! no magnetic field
				particles(i,j)%vx = particles(i,j)%vx - 0.5*a*dt
			END IF	
			
		    !PRINT*,'i,j,vx,vy = ',i,j,particles(i,j)%vx,particles(i,j)%vy  

		END DO	
	END DO


END SUBROUTINE setv


! adds sinusoidal perturbations to x(0) and v(0) to electrons
SUBROUTINE add_perturbations()

	INTEGER, PARAMETER :: modes = nx-1
    
	REAL*8 :: x1(modes) 
	REAL*8 :: v1(modes)
	REAL*8 :: thetax(modes) 
	REAL*8 :: thetav(modes) 
	INTEGER :: i, j
    REAL*8 :: theta, p
	
	x1 = 0.d0
	x1(1) = 0.1d0
	!x1(17) = 0.01d0
    !x1(32) = 0.01d0
	v1 = 0.d0
	!v1(1) = 0.5d0
	thetax = 0.d0
	thetav = 0.d0

    !DO i = 1,modes-1
	!    CALL RANDOM_NUMBER(p)
	!	v1(i) = p*0.1d0
	!END DO
	
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