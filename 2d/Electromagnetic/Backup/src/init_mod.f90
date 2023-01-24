MODULE init_mod

USE constants_mod
USE data_mod
USE particleMover_mod
IMPLICIT NONE

CONTAINS

! intial particle spatial distribution function (same for all species)
FUNCTION fn0(x,y) RESULT(fx)

	REAL*8, INTENT(IN) :: x,y
	INTEGER :: fx

    ! uniform spatial distribution
    fx = N/(nx)!*ny)

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

	INTEGER ::  i, j, k, counter, N_cell, i_v, jmin, jmax
	REAL*8 :: wp, qm, q, m, p


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
		
		q = (wp**2)*eps0*(Lx*Ly/N)/qm
		m = q/qm
		
		PRINT*,'Species# ',i,' q,m = ',q,m
		
		!rho0 = -2.d0*q*(N/(nx*ny))/(dx*dy)  ! neutralizing background charge
		
	    DO j = 1, N
		    particles(i,j)%q = q
		    particles(i,j)%m = m 
			
		END DO
		
 	END DO	

    ! distribute particles across grid
	counter = 1
	
    jmin =ny/2  ! infinite line current
    jmax =ny/2
    
    DO i = 1, nx
        DO j = jmin, jmax
            ! find out how many particles to put in this cell
            N_cell = fn0(i*dx,j*dy)
		
            !PRINT*,'Cell#, N_cell',i,j,N_cell
		
            IF(counter .GT. N) EXIT	
		
            ! disrtribute particles uniformly across the cell
            DO k = 1, N_cell	

                ! set initial position : x(0)		
                particles(1,counter)%x = (i-0.5d0)*dx + k*dx/(N_cell+1) !- 0.1*dx/N_cell   
                particles(1,counter)%y = (j-0.5d0)*dy + 0.25*dy !+ k*dy/N_cell - 0.1*dy/N_cell   
               
               
                ! set initial velocity : v(0) sampled uniformly from (vmin,vmax)
                particles(1,counter)%ux = u0x !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5))
                particles(1,counter)%uy = u0y !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5))
             
                IF(ns .EQ. 2) THEN
                    particles(2,counter)%x = (i-0.5d0)*dx + k*dx/(N_cell+1) !- 0.1*dx/N_cell
                    particles(2,counter)%y = (j-0.5d0)*dy + 0.25*dy !+ k*dy/N_cell - 0.1*dy/N_cell   
                
                    particles(2,counter)%ux = u0x !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5)) 
                    particles(2,counter)%uy = u0y !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5))
                END IF
                
                counter = counter + 1
            END DO
        
        END DO
	END DO

        !particles(1,1)%x = (nx/2)*dx    
        !particles(1,1)%y = (ny/2)*dy    
        
        !particles(1,1)%ux = u0x !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5))
        !particles(1,1)%uy = u0y !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5))
             
       
        !IF(ns .EQ. 2) THEN
        !    particles(2,1)%x = (nx/2)*dx
        !    particles(2,1)%y = (ny/2)*dy   
       ! 
        !    particles(2,1)%ux = 0.d0! u0x !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5)) 
        !    particles(2,1)%uy = u0y !+ fv0(vmin+((vmax-vmin)/N_cell)*(k-0.5))
        !END IF    
            
            
    IF(randPos) THEN
    DO i = 1, ns
		
	    DO j = 1, N
            CALL RANDOM_NUMBER(p)
		    particles(i,j)%x = xmin+p*Lx
			 CALL RANDOM_NUMBER(p)
		    particles(i,j)%y = ymin+p*Ly
		END DO
		
 	END DO	 
    END IF           
               
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
			PRINT*,'PARTICLE TYPE, ID, x, y, m, q = ',i,j,particles(i,j)%x, particles(i,j)%y, &
			                                           particles(i,j)%m, particles(i,j)%q 
		END DO
	END DO
	PRINT*,''
	END IF

END SUBROUTINE init_particles


! accumulates charges at grid-points via interpolation
SUBROUTINE setrho()

    INTEGER :: i, j, ix, iy
	REAL*8 :: qdxy, x, y

    ! clear density
    rho(:,:) = rho0
	ne(:,:) = 0.d0
	ni(:,:) = 0.d0
	rho(0,:) = 0.d0
	rho(nx+1,:) = 0.d0
	rho(:,0) = 0.d0
	rho(:,ny+1) = 0.d0
    
	PRINT*,''
	 
	! accumulate particle chrages at grid points
	
	! zeroth order interpolation ("all or nothing..")
    IF(interpolation_type .EQ. 1) THEN
	    DO i = 1, ns
		    qdxy = particles(i,1)%q/(dx*dy)
			
			DO j = 1, N
			
				! apply periodic boundary condition
				IF(particles(i,j)%x .LT. xmin) THEN
					particles(i,j)%x = particles(i,j)%x + Lx
				END IF
				IF(particles(i,j)%x .GT. xmax) THEN
					particles(i,j)%x = particles(i,j)%x - Lx
				END IF
                IF(particles(i,j)%y .LT. ymin) THEN
					particles(i,j)%y = particles(i,j)%y + Ly
				END IF
				IF(particles(i,j)%y .GT. ymax) THEN
					particles(i,j)%y = particles(i,j)%y - Ly
				END IF
                
				x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy                
				
				!compute index of occupied cell
				ix = x + 0.5  
                iy = y + 0.5
                
				!PRINT*,'i,j,x,y,ix,iy,rho(ix,iy)=',i,j,x,y,ix,iy,rho(ix,iy)

				! accumulate charge in the cell occupied by the particle (i.e. nearest grid point)
				rho(ix,iy) = rho(ix,iy) + qdxy
								
                IF(i .EQ. 1) ne(ix,iy) = ne(ix,iy) + 1.0				
				IF(i .EQ. 2) ni(ix,iy) = ni(ix,iy) + 1.0			
				 
		    END DO
		END DO 
	
	
	! bi-linear interpolation
	ELSE IF(interpolation_type .EQ. 2 .OR. interpolation_type .EQ. 3) THEN

	    DO i = 1, ns
		    qdxy = particles(i,1)%q/(dx*dy)
			
			DO j = 1, N
				
				! apply periodic boundary condition
				IF(particles(i,j)%x .LT. xmin) THEN
					particles(i,j)%x = particles(i,j)%x + Lx
				END IF
				IF(particles(i,j)%x .GT. xmax) THEN
					particles(i,j)%x = particles(i,j)%x - Lx
				END IF
                IF(particles(i,j)%y .LT. ymin) THEN
					particles(i,j)%y = particles(i,j)%y + Ly
				END IF
				IF(particles(i,j)%y .GT. ymax) THEN
					particles(i,j)%y = particles(i,j)%y - Ly
				END IF
				
				x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy                

				!compute index of nearest bottom-left grid point
				ix = x 
                iy = y
				
				!PRINT*,'i,j,x,y,ix,iy,rho(ix,iy)=',i,j,x,y,ix,iy,rho(ix,iy)
	
				! accumulate charge in the cells occupied by the cloud (weighted according to the fraction of the cell occupied by the cloud)
				rho(ix,iy)     = rho(ix,iy) + qdxy*(1.d0 - (x-ix))*(1.d0 - (y-iy))				
				rho(ix+1,iy)   = rho(ix+1,iy) +  qdxy*(x - ix)*(1.d0 - (y-iy)) 
                rho(ix,iy+1)   = rho(ix,iy+1) + qdxy*(1.d0 - (x-ix))*(y - iy)	
                rho(ix+1,iy+1) = rho(ix+1,iy+1) + qdxy*(x - ix)*(y - iy)	
                
			    IF(i .EQ. 1)THEN
                    ne(ix,iy)     = ne(ix,iy) + (1.d0 - (x-ix))*(1.d0 - (y-iy))				
                    ne(ix+1,iy)   = ne(ix+1,iy) +  (x - ix)*(1.d0 - (y-iy)) 
                    ne(ix,iy+1)   = ne(ix,iy+1) + (1.d0 - (x-ix))*(y - iy)	
                    ne(ix+1,iy+1) = ne(ix+1,iy+1) + (x - ix)*(y - iy)	
                END IF
                
                IF(i .EQ. 2)THEN
                    ni(ix,iy)     = ni(ix,iy) + (1.d0 - (x-ix))*(1.d0 - (y-iy))				
                    ni(ix+1,iy)   = ni(ix+1,iy) +  (x - ix)*(1.d0 - (y-iy)) 
                    ni(ix,iy+1)   = ni(ix,iy+1) + (1.d0 - (x-ix))*(y - iy)	
                    ni(ix+1,iy+1) = ni(ix+1,iy+1) + (x - ix)*(y - iy)	
                END IF
                
		    END DO
		END DO 
	END IF
	
	rho(nx,:) = rho(nx,:) + rho(0,:) 
	rho(1,:) = rho(1,:) + rho(nx+1,:)
	ne(nx,:) = ne(nx,:) + ne(0,:)
	ne(1,:) = ne(1,:) + ne(nx+1,:)
	ni(nx,:) = ni(nx,:) + ni(0,:)
	ni(1,:) = ni(1,:) + ni(nx+1,:)
	
    rho(:,ny) = rho(:,ny) + rho(:,0) 
	rho(:,1) = rho(:,1) + rho(:,ny+1)
	ne(:,ny) = ne(:,ny) + ne(:,0)
	ne(:,1) = ne(:,1) + ne(:,ny+1)
	ni(:,ny) = ni(:,ny) + ni(:,0)
	ni(:,1) = ni(:,1) + ni(:,ny+1)
	
    
    ! apply periodic boundary conditions
	rho(0,:) = rho(nx,:)
    rho(nx+1,:) = rho(1,:)	
	ne(0,:) = ne(nx,:)
    ne(nx+1,:) = ne(1,:)	
	ni(0,:) = ni(nx,:)
    ni(nx+1,:) = ni(1,:)	
	
	rho(:,0) = rho(:,ny)
    rho(:,ny+1) = rho(:,1)	
	ne(:,0) = ne(:,ny)
    ne(:,ny+1) = ne(:,1)	
	ni(:,0) = ni(:,ny)
    ni(:,ny+1) = ni(:,1)	
	
    
	IF(print_debug) THEN
	PRINT*,''
	PRINT*,'RHO(t=0) = '
	DO j = ny, 1, -1
        DO i = 1, nx
            WRITE(*,FMT='(1f8.2)',ADVANCE='NO') rho(i,j)
        END DO
        PRINT*,''
    END DO    
	END IF
	
END SUBROUTINE setrho

	
! sets v(0) to v(-dt/2) as required by the leap-frog integrator
SUBROUTINE setv()

	INTEGER :: i, j
	REAL*8 :: qm, wc, x, y, a(2), uxx
    REAL*8 :: b, t, s, gam

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
			y = particles(i,j)%y/dy 

		    CALL accel(x, y, qm, a, b)			
			
			! half electric push (u^n-1/2 ->u^-) 
			particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*(-0.5*dt)
			particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*(-0.5*dt)
						
            gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2)/c**2 )
            t = 0.5d0*b*(-0.5*dt)/(gam*c)
            s = 2.d0*t/(1+t**2)
			uxx = particles(i,j)%ux
                                
			! magnetic rotation (u^- -> u^+)
			particles(i,j)%ux = (1.d0-2.d0*s*t)*particles(i,j)%ux + s*particles(i,j)%uy
			particles(i,j)%uy = -s*uxx + (1.d0-2.d0*s*t)*particles(i,j)%uy
			
            
			! remaining half electric push (u^+ -> u^n+1/2)		
			particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*(-0.5*dt)
			particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*(-0.5*dt)
			
		    !PRINT*,'i,j,ux,uy = ',i,j,particles(i,j)%ux,particles(i,j)%uy  

		END DO	
	END DO


END SUBROUTINE setv


! adds sinusoidal perturbations to x(0) and v(0) to electrons
SUBROUTINE add_perturbations()

    
	REAL*8 :: x1(-nx/2:nx/2-1) 
	REAL*8 :: v1(-nx/2:nx/2-1)
	REAL*8 :: thetax(-nx/2:nx/2-1) 
	REAL*8 :: thetav(-nx/2:nx/2-1) 
	INTEGER :: i, j
    REAL*8 :: theta, p
	
	x1 = 0.d0
	!x1(1) = 0.05d0
	!x1(17) = 0.01d0
    !x1(32) = 0.01d0
	v1 = 0.d0
	!v1(1) = 0.001*c
	thetax = 0.d0
	thetav = 0.d0

    !DO i = 1,modes-1
	!    CALL RANDOM_NUMBER(p)
	!	v1(i) = p*0.1d0
	!END DO
	
	DO i = 1, N
	    DO j = -nx/2, nx/2-1
			theta = particles(1,i)%x*j*twopi/Lx
			particles(1,i)%x = particles(1,i)%x + x1(j)*COS(theta+thetax(j))
			particles(1,i)%ux = particles(1,i)%ux + v1(j)*COS(theta+thetav(j))
		END DO
		
		!PRINT*,'Electron ID, x, ux = ',i,particles(1,i)%x, particles(1,i)%ux 
	END DO



END SUBROUTINE add_perturbations


END MODULE init_mod