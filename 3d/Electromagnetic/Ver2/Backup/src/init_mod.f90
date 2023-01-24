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
    fx = N/(nx)!*ny)  for line current

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
	REAL*8 :: wp, qm, q, m, p, gam

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
		
		q = (wp**2)*eps0*(Lx*Ly*Lz/N)/qm
		m = q/qm
		
		PRINT*,'Species# ',i,' q,m = ',q,m
		
		!rho0 = -2.d0*q*(N/(nx*ny))/(dx*dy)  ! neutralizing background charge
		
	    DO j = 1, N
		    particles(i,j)%q = q
		    particles(i,j)%m = m 
			
		END DO
		
 	END DO	



    GO TO 101

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

    101 CONTINUE
     
    gam = 1.d0 / SQRT(1.d0 - (u0x**2+u0y**2+u0z**2)/c**2 ) 
     
    particles(1,1)%x = (nx/4)*dx    
    particles(1,1)%y = (ny/4)*dy   
    particles(1,1)%z = (nz/2)*dz
        
    particles(1,1)%ux = gam*u0x 
    particles(1,1)%uy = gam*u0y 
    particles(1,1)%uz = gam*u0z
    
       
    IF(ns .EQ. 2) THEN
        particles(2,1)%x = (nx/4)*dx    
        particles(2,1)%y = (ny/4)*dy    
        particles(2,1)%z = (nz/2)*dz
        particles(2,1)%ux = 0.d0 
        particles(2,1)%uy = 0.d0 
        particles(2,1)%uz = 0.d0
    END IF    
            
            
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

    INTEGER :: i, j, ix, iy, iz
	REAL*8 :: qdxyz, x, y, z
	REAL*8 :: delx, dely, delz, dxx, dyy, dzz

     ! clear density
    rho = 0.d0
    !rho(1:nx,1:ny,1:nz) = rho0
	ne = 0.d0
	ni = 0.d0
    
	PRINT*,''
	 
	! accumulate particle chrages at grid points
	
	! zeroth order interpolation ("all or nothing..")
    IF(interpolation_type .EQ. 1) THEN
	    DO i = 1, ns
		    qdxyz = particles(i,1)%q/(dx*dy*dz)
			
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
                IF(particles(i,j)%z .LT. zmin) THEN
                    particles(i,j)%z = particles(i,j)%z + Lz
                END IF
                IF(particles(i,j)%z .GT. zmax) THEN
                    particles(i,j)%z = particles(i,j)%z - Lz
                END IF
        
				x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy                
                z = particles(i,j)%y/dz                
                
                !compute index of occupied cell
                ix = x + 0.5  
                iy = y + 0.5
                iz = z + 0.5
                
                    !PRINT*,'i,j,x,y,ix,iy,rho(ix,iy)=',i,j,x,y,ix,iy,rho(ix,iy)

                ! accumulate charge in the cell occupied by the particle (i.e. nearest grid point)
                rho(ix,iy,iz) = rho(ix,iy,iz) + qdxyz
								
                IF(i .EQ. 1) ne(ix,iy,iz) = ne(ix,iy,iz) + 1.0				
                IF(i .EQ. 2) ni(ix,iy,iz) = ni(ix,iy,iz) + 1.0			
				 
		    END DO
		END DO 
	
	
	! bi-linear interpolation
	ELSE IF(interpolation_type .EQ. 2 .OR. interpolation_type .EQ. 3) THEN

	    DO i = 1, ns
		    qdxyz = particles(i,1)%q/(dx*dy*dz)
			
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
				IF(particles(i,j)%z .LT. zmin) THEN
                    particles(i,j)%z = particles(i,j)%z + Lz
                END IF
                IF(particles(i,j)%z .GT. zmax) THEN
                    particles(i,j)%z = particles(i,j)%z - Lz
                END IF
                
				x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy                
                z = particles(i,j)%z/dz                

                !compute index of nearest bottom-left grid point
                ix = x 
                iy = y
                iz = z
					         
                delx = x-ix
                dely = y-iy
                delz = z-iz
                dxx = 1.d0 - delx
                dyy = 1.d0 - dely
                dzz = 1.d0 - delz
        
                ! accumulate charge in cells occupied by the cloud (weighted according to the fraction of the cell occupied by the cloud)
                rho(ix,iy,iz)       = rho(ix,iy,iz)       + qdxyz*dxx*dyy*dzz				     
                rho(ix+1,iy,iz)     = rho(ix+1,iy,iz)     + qdxyz*delx*dyy*dzz
                rho(ix,iy+1,iz)     = rho(ix,iy+1,iz)     + qdxyz*dxx*dely*dzz
                rho(ix,iy,iz+1)     = rho(ix,iy,iz+1)     + qdxyz*dxx*dyy*delz
                rho(ix+1,iy+1,iz)   = rho(ix+1,iy+1,iz)   + qdxyz*delx*dely*dzz
                rho(ix+1,iy,iz+1)   = rho(ix+1,iy,iz+1)   + qdxyz*delx*dyy*delz
                rho(ix,iy+1,iz+1)   = rho(ix,iy+1,iz+1)   + qdxyz*dxx*dely*delz
                rho(ix+1,iy+1,iz+1) = rho(ix+1,iy+1,iz+1) + qdxyz*delx*dely*delz				
	
                
                IF(i .EQ. 1)THEN
                    ne(ix,iy,iz)       = ne(ix,iy,iz)       + qdxyz*dxx*dyy*dzz				     
                    ne(ix+1,iy,iz)     = ne(ix+1,iy,iz)     + qdxyz*delx*dyy*dzz
                    ne(ix,iy+1,iz)     = ne(ix,iy+1,iz)     + qdxyz*dxx*dely*dzz
                    ne(ix,iy,iz+1)     = ne(ix,iy,iz+1)     + qdxyz*dxx*dyy*delz
                    ne(ix+1,iy+1,iz)   = ne(ix+1,iy+1,iz)   + qdxyz*delx*dely*dzz
                    ne(ix+1,iy,iz+1)   = ne(ix+1,iy,iz+1)   + qdxyz*delx*dyy*delz
                    ne(ix,iy+1,iz+1)   = ne(ix,iy+1,iz+1)   + qdxyz*dxx*dely*delz
                    ne(ix+1,iy+1,iz+1) = ne(ix+1,iy+1,iz+1) + qdxyz*delx*dely*delz				
                END IF
                
                IF(i .EQ. 2)THEN
                    ni(ix,iy,iz)       = ni(ix,iy,iz)       + qdxyz*dxx*dyy*dzz				     
                    ni(ix+1,iy,iz)     = ni(ix+1,iy,iz)     + qdxyz*delx*dyy*dzz
                    ni(ix,iy+1,iz)     = ni(ix,iy+1,iz)     + qdxyz*dxx*dely*dzz
                    ni(ix,iy,iz+1)     = ni(ix,iy,iz+1)     + qdxyz*dxx*dyy*delz
                    ni(ix+1,iy+1,iz)   = ni(ix+1,iy+1,iz)   + qdxyz*delx*dely*dzz
                    ni(ix+1,iy,iz+1)   = ni(ix+1,iy,iz+1)   + qdxyz*delx*dyy*delz
                    ni(ix,iy+1,iz+1)   = ni(ix,iy+1,iz+1)   + qdxyz*dxx*dely*delz
                    ni(ix+1,iy+1,iz+1) = ni(ix+1,iy+1,iz+1) + qdxyz*delx*dely*delz		
                END IF
                            
		    END DO
		END DO 
	END IF
	
    rho(nx,:,:) = rho(nx,:,:) + rho(0,:,:) 
    rho(1,:,:) = rho(1,:,:) + rho(nx+1,:,:)
    rho(:,ny,:) = rho(:,ny,:) + rho(:,0,:) 
    rho(:,1,:) = rho(:,1,:) + rho(:,ny+1,:)
    rho(:,:,nz) = rho(:,:,nz) + rho(:,:,0) 
    rho(:,:,1) = rho(:,:,1) + rho(:,:,nz+1)
                     
    ne(nx,:,:) = ne(nx,:,:) + ne(0,:,:) 
    ne(1,:,:) = ne(1,:,:) + ne(nx+1,:,:)
    ne(:,ny,:) = ne(:,ny,:) + ne(:,0,:) 
    ne(:,1,:) = ne(:,1,:) + ne(:,ny+1,:)
    ne(:,:,nz) = ne(:,:,nz) + ne(:,:,0) 
    ne(:,:,1) = ne(:,:,1) + ne(:,:,nz+1)
        
    ni(nx,:,:) = ni(nx,:,:) + ni(0,:,:) 
    ni(1,:,:) = ni(1,:,:) + ni(nx+1,:,:)
    ni(:,ny,:) = ni(:,ny,:) + ni(:,0,:) 
    ni(:,1,:) = ni(:,1,:) + ni(:,ny+1,:)
    ni(:,:,nz) = ni(:,:,nz) + ni(:,:,0) 
    ni(:,:,1) = ni(:,:,1) + ni(:,:,nz+1)
	
    
	IF(print_debug) THEN
	PRINT*,''
	PRINT*,'RHO(t=0) = '
	DO j = ny, 1, -1
        DO i = 1, nx
            WRITE(*,FMT='(1f8.2)',ADVANCE='NO') rho(i,j,1)
        END DO
        PRINT*,''
    END DO    
	END IF
	
END SUBROUTINE setrho

	
! sets v(0) to v(-dt/2) as required by the leap-frog integrator
SUBROUTINE setv()

	INTEGER :: i, j    
    REAL*8 :: qm, wc, x, y, z, a(3), b(3), qdxyz, temp, uxx, uyy, uzz
    REAL*8 :: t(3), s(3), gam, tsqr1, udots, tdots
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
            z = particles(i,j)%z/dz
            
		    CALL accel(x, y, z, qm, a, b)						
                   
            ! half electric push (u^n-1/2 ->u^-) 
            particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*(-0.5*dt)
            particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*(-0.5*dt)
            particles(i,j)%uz = particles(i,j)%uz + 0.5*a(3)*(-0.5*dt)
						
            gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2 + particles(i,j)%uz**2 )/c**2 )
            t     = 0.5d0 * b * (-0.5*dt) / (gam*c)
            tsqr1 = 1.d0 + SUM(t**2)
            s     = 2.d0  * t  / tsqr1
            uxx   = particles(i,j)%ux
            uyy   = particles(i,j)%uy
            uzz   = particles(i,j)%uz
            udots = uxx*s(1) + uyy*s(2) + uzz*s(3)
            tdots = 2.d0*SUM(t**2) / tsqr1
                
            ! magnetic rotation (u^- -> u^+)
            particles(i,j)%ux = (1.d0-tdots)*uxx + udots*t(1) + (s(3)*uyy - s(2)*uzz)  
            particles(i,j)%uy = (1.d0-tdots)*uyy + udots*t(2) + (s(1)*uzz - s(3)*uxx)
            particles(i,j)%uz = (1.d0-tdots)*uzz + udots*t(3) + (s(2)*uxx - s(1)*uyy)
			
            
            ! remaining half electric push (u^+ -> u^n+1/2)		
            particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*(-0.5*dt)
            particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*(-0.5*dt)
            particles(i,j)%uz = particles(i,j)%uz + 0.5*a(3)*(-0.5*dt)
                             
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