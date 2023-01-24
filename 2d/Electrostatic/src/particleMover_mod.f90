MODULE particleMover_mod

USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


SUBROUTINE accel(x, y, qm, a)

	REAL*8, INTENT(IN) :: x, y, qm
	REAL*8, INTENT(INOUT) :: a(2)
	INTEGER :: ix, iy

	IF(interpolation_type .EQ. 1 .OR. interpolation_type .EQ. 3) THEN
		ix = x + 0.5
        iy = y + 0.5
        
		a(1) = qm*Ex(ix,iy)
		a(2) = qm*Ey(ix,iy)
        
	ELSE IF(interpolation_type .EQ. 2) THEN
	
		ix = x 
        iy = y
        
		a(1) = qm*( (1.d0-(x-ix))*(1.d0-(y-iy))*Ex(ix,iy) + (x-ix)*(1.d0-(y-iy))*Ex(ix+1,iy) + &
                    (1.d0-(x-ix))*(y-iy)*Ex(ix,iy+1) + (x-ix)*(y-iy)*Ex(ix+1,iy+1)  )
                    
		a(2) = qm*( (1.d0-(x-ix))*(1.d0-(y-iy))*Ey(ix,iy) + (x-ix)*(1.d0-(y-iy))*Ey(ix+1,iy) + &
                    (1.d0-(x-ix))*(y-iy)*Ey(ix,iy+1) + (x-ix)*(y-iy)*Ey(ix+1,iy+1)  )
        
	END IF

END SUBROUTINE accel


! updates velocities and positions of particles (using leap-frog integration) and acculmulates charge at the grid-points
SUBROUTINE move_particles()

	INTEGER :: i, j, ix, iy
	REAL*8 :: qm, wc, x, y, a(2), qdxy, temp, vxx
		
	! clear kinitic energy and momentum accumulator	
	KE = 0.d0
	Px = 0.d0
	Py = 0.d0
	
    ! clear density
    rho(:,:) = rho0
	ne(:,:) = 0.d0
	ni(:,:) = 0.d0
	rho(0,:) = 0.d0
	rho(nx+1,:) = 0.d0
	rho(:,0) = 0.d0
	rho(:,ny+1) = 0.d0
	
	DO i= 1,ns   
	
		IF(i .EQ. 1) THEN
			qm = qm_e
			wc = wc_e
		ELSE IF(i .EQ. 2) THEN 	
			qm = qm_i
			wc = wc_i
		END IF
		qdxy = particles(i,1)%q/(dx*dy)
        
		DO j = 1, N

			x = particles(i,j)%x/dx
			y = particles(i,j)%y/dy
            
			temp =  SQRT(particles(i,j)%vx**2+particles(i,j)%vy**2)
		
			!******************************
			! update v(t-dt/2) -> v(t+dt/2)
		    !******************************
			CALL accel(x, y, qm, a)			
						
			IF(ABS(wc) .GT. 1.d-15) THEN	! w/ magnetic field 	
		
				! half acceleration
				particles(i,j)%vx = particles(i,j)%vx + 0.5*a(1)*dt
				particles(i,j)%vy = particles(i,j)%vy + 0.5*a(2)*dt
						
			    vxx = particles(i,j)%vx
				! rotation
				particles(i,j)%vx = COS(wc*dt)*particles(i,j)%vx - SIN(wc*dt)*particles(i,j)%vy
				particles(i,j)%vy = SIN(wc*dt)*vxx + COS(wc*dt)*particles(i,j)%vy
			
     			temp =  particles(i,j)%vx**2+particles(i,j)%vy**2

				! half acceleration			
				particles(i,j)%vx = particles(i,j)%vx + 0.5*a(1)*dt
				particles(i,j)%vy = particles(i,j)%vy + 0.5*a(2)*dt
                
			ELSE ! no magnetic field
		
				particles(i,j)%vx = particles(i,j)%vx + a(1)*dt
				particles(i,j)%vy = particles(i,j)%vy + a(2)*dt
                
				temp =  temp*SQRT(particles(i,j)%vx**2+particles(i,j)%vy**2)
		
			END IF
		    			
		    KE = KE + 0.5*particles(i,j)%m*temp 
			Px = Px + particles(i,j)%m*particles(i,j)%vx
			Py = Py + particles(i,j)%m*particles(i,j)%vy
			
			!***********************
			! update x(t) -> x(t+dt)
		    !***********************
		    particles(i,j)%x = particles(i,j)%x + particles(i,j)%vx*dt
			particles(i,j)%y = particles(i,j)%y + particles(i,j)%vy*dt
			
			
			!******************
			! accumulate charge
		    !******************
			!x = particles(i,j)%x/dx 
			
		    IF(interpolation_type .EQ. 1) THEN			
							
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
				
			ELSE IF(interpolation_type .EQ. 2 .OR. interpolation_type .EQ. 3) THEN
			    
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
                
			END IF
		

     		!PRINT*,'PARTICLE TYPE, ID, x, y, vx, vy = ',i,j,particles(i,j)%x, particles(i,j)%y, particles(i,j)%vx,particles(i,j)%vy 

		END DO
	END DO

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
	PRINT*,'RHO = '
	DO j = ny+1, 0, -1
        DO i = 0, nx+1
            WRITE(*,FMT='(1f5.1)',ADVANCE='NO') rho(i,j)
        END DO
        PRINT*,''
    END DO    
	END IF
	
END SUBROUTINE move_particles




END MODULE particleMover_mod