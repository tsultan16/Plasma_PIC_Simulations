MODULE particleMover_mod

USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


SUBROUTINE accel(x, qm, a, bz)

	REAL*8, INTENT(IN) :: x, qm
	REAL*8, INTENT(INOUT) :: a(2), bz
	INTEGER :: ix, ileft, iright

	IF(interpolation_type .EQ. 1 .OR. interpolation_type .EQ. 3) THEN
		ix = x + 0.5
		a(1) = qm*Ex(ix)
		a(2) = qm*Ey(ix) 
		bz = Bz(ix)
		
	ELSE IF(interpolation_type .EQ. 2) THEN
	
		ix = x 
		a(1) = qm*( (1.d0-(x-ix))*Ex(ix) + (x-ix)*Ex(ix+1) )
		a(2) = qm*( (1.d0-(x-ix))*Ey(ix) + (x-ix)*Ey(ix+1) )
		bz =  (1.d0-(x-ix))*Bz(ix) + (x-ix)*Bz(ix+1)
		
	END IF

END SUBROUTINE accel


! updates velocities and positions of particles (using leap-frog integration) and acculmulates charge at the grid-points
SUBROUTINE move_particles()

	INTEGER :: i, j, ix, ileft, iright
	REAL*8 :: qm, wc, x, a(2), qdx, temp, uxx, bz, w, t, s, gam
		
	! clear kinitic energy and momentum accumulator	
	KE = 0.d0
	Px = 0.d0
	Py = 0.d0
	
    ! clear density and current
    rho = rho0	
	ne(:) = 0.d0
	ni(:) = 0.d0
	rho(0) = 0.d0
	rho(nx+1) = 0.d0
	Jy = 0.d0
	Jy(0) = 0.d0
	Jy(nx+1) = 0.d0
	
	DO i= 1,ns   
	
		IF(i .EQ. 1) THEN
			qm = qm_e
			wc = wc_e
		ELSE IF(i .EQ. 2) THEN 	
			qm = qm_i
			wc = wc_i
		END IF
		qdx = particles(i,1)%q/dx
        
		DO j = 1, N

			x = particles(i,j)%x/dx 
		
			!******************************
			! update v(t-dt/2) -> v(t+dt/2)
		    !******************************
			CALL accel(x, qm, a, bz)			
			
			! half acceleration
			particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*dt
            particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*dt				
			
			gam = SQRT(1.d0+(particles(i,j)%ux**2+particles(i,j)%uy**2)/(c**2))
			w = qm*bz + wc  
			t = w*dt/(2.d0*c*gam)
			s = 2*t/(1+t**2)
			
			uxx = particles(i,j)%ux + t*particles(i,j)%uy
			
			! rotation
			particles(i,j)%uy = particles(i,j)%uy - s*uxx
			particles(i,j)%ux = uxx + t*particles(i,j)%uy
			

			! half acceleration			
			particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*dt
            particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*dt				
			
						
		    !KE = KE + 0.5*particles(i,j)%m*temp 
			!Px = Px + particles(i,j)%m*particles(i,j)%vx
			!Py = Py + particles(i,j)%m*particles(i,j)%vy
			
			!***********************
			! update x(t) -> x(t+dt)
		    !***********************
			gam = SQRT(1.d0+(particles(i,j)%ux**2+particles(i,j)%uy**2)/(c**2))
			
		    particles(i,j)%x = particles(i,j)%x + particles(i,j)%ux*dt/gam
			particles(i,j)%y = particles(i,j)%y + particles(i,j)%uy*dt/gam
			
			
			!******************************************
			! deposit charge and current to grid points
		    !******************************************
			
		    IF(interpolation_type .EQ. 1) THEN			
							
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
                
				! accumulate in the cell occupied by the particle
				rho(ix) = rho(ix) + qdx
				
                Jy(ix) = Jy(ix) + qdx*dx*particles(i,j)%vy							
							
                IF(i .EQ. 1) ne(ix) = ne(ix) + 1.0				
				IF(i .EQ. 2) ni(ix) = ni(ix) + 1.0			
				
			ELSE IF(interpolation_type .EQ. 2 .OR. interpolation_type .EQ. 3) THEN
			    
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
			
				! accumulate in the cells occupied by the cloud (weighted according to the fraction of the cell occupied by the cloud)
				rho(ix)  = rho(ix) + qdx*(1.d0 - (x-ix))
				rho(ix+1) = rho(ix+1) +  qdx*(x - ix)
				
				Jy(ix) = Jy(ix) + qdx*dx*particles(i,j)%vy*(1.d0 - (x-ix))	
				Jy(ix+1) = Jy(ix+1) + qdx*dx*particles(i,j)%vy*(x - ix)	
				
				IF(i .EQ. 1) ne(ix) = ne(ix) + (1.d0 - (x-ix))				
				IF(i .EQ. 2) ni(ix) = ni(ix) + (x - ix)
				
			END IF
		

     		!PRINT*,'PARTICLE TYPE, ID, x, vx, vy = ',i,j,particles(i,j)%x, particles(i,j)%vx,particles(i,j)%vy 

		END DO
	END DO
	
	rho(nx) = rho(nx) + rho(0) 
	rho(1) = rho(1) + rho(nx+1)
	Jy(nx) = Jy(nx) + Jy(0) 
	Jy(1) = Jy(1) + Jy(nx+1)
	ne(nx) = ne(nx) + ne(0)
	ne(1) = ne(1) + ne(nx+1)
	ni(nx) = ni(nx) + ni(0)
	ni(1) = ni(1) + ni(nx+1)
	
    ! apply periodic boundary conditions
	rho(0) = rho(nx)
    rho(nx+1) = rho(1)
    Jy(0) = Jy(nx)
    Jy(nx+1) = Jy(1)    	
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
	
END SUBROUTINE move_particles




END MODULE particleMover_mod