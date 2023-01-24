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
    !fx = Np_in/(nx*ny*nz)  !for line current

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

	INTEGER ::  i, j, k, l
    INTEGER, PARAMETER :: npass = 2
	REAL*8 :: p(3), gam, x, y, z, uxe, uye, uze, uxi, uyi, uzi
    
    
    ! set mass and charge
    q = (/ q_e, q_i /)
    m = (/ m_e, m_i /)
    
    PRINT*,'q/m : electrons, ions =', q(1)/m(1), q(2)/m(2)
    
    GO TO 111
    
    DO l = 1, npass
    
    ! distribute particles across grid (1 pair per cell)    
    DO k = 2, nz-2, 1 
        DO j = 2, ny-2, 1
            DO i = 2, nx-2, 1
	
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 1*p(1)
            y = j + 1*p(2) 
            z = k + 1*p(3) 
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            uxe = v0x + vth_e * (vth_e + SUM(p))
                                        
            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5
            uye = vth_e * SUM(p)

            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5                
            uze = vth_e * SUM(p)
               
            gam = 1.d0 / SQRT(1.d0 - (uxe**2+uye**2+uze**2)/c**2 ) 
            uxe = gam * uxe
            uye = gam * uye
            uze = gam * uze
               
            IF(gam .LT. 1.D0) THEN
                PRINT*,'ERROR: Gamma < 1...'
                STOP   
            END IF    

                     
            CALL create_particle(1,x,y,z,uxe,uye,uze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
               
            IF(ns .EQ. 2) THEN


                CALL RANDOM_NUMBER(p)                                      

                p = p - 0.5
                


                uxi = v0x + vth_i*(vth_i + SUM(p))
                
                CALL RANDOM_NUMBER(p)                                      
                p = p - 0.5
                uyi = vth_i * SUM(p)

 

                CALL RANDOM_NUMBER(p)                    
                p = p - 0.5                  
                uzi = vth_i * SUM(p)


                gam = 1.d0 / SQRT(1.d0 - (uxi**2+uyi**2+uzi**2)/c**2 ) 
                uxi = gam * uxi
                uyi = gam * uyi
                uzi = gam * uzi
            
                IF(gam .LT. 1.D0) THEN
                    PRINT*,'ERROR: Gamma < 1...'
                    STOP   
                END IF    
            

                CALL create_particle(2,x,y,z,uxi,uyi,uzi)
    
                Ntot_inj(2) = Ntot_inj(2) + 1

            END IF
            
            END DO
        END DO
	END DO


    END DO
    
    111 CONTINUE     
   
    IF(myrank .EQ. 0) THEN
        CALL create_particle(1,DBLE(nx-2),DBLE(ny-2),DBLE(nz/2),0.1*c,0.d0,0.d0)
        CALL create_particle(2,DBLE(nx-2),DBLE(ny-2),DBLE(nz/2),-0.1*c,0.d0,0.d0)
    END IF
   
    WRITE(*,FMT = '(" Inserted ",i8," particle pairs.")') Np_in(1)

         
          
               
    PRINT*,''
	!DO i = 1, Np_in(1)
	!    WRITE(*,FMT='("Electron # ", I4, " x, y, z = ", F5.2, 2X, F5.2, 2X, F5.2 )') &
    !          i,particles(1,i)%x,particles(1,i)%y,particles(1,i)%z 
	!    WRITE(*,FMT='("Ion      # ", I4, " x, y, z = ", F5.2, 2X, F5.2, 2X, F5.2 )') &
    !          i,particles(2,i)%x,particles(2,i)%y,particles(2,i)%z 
    !END DO
    !PRINT*,''

    PRINT*,'Initialization Completed...'

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
		    qdxyz = q(i)/(dx*dy*dz)
			
			DO j = 1, Np_in(1)
			
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
		    qdxyz = q(i)/(dx*dy*dz)
			
			DO j = 1, Np_in(1)
				
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
	   
        qm = q(i)/m(i)
		
		DO  j = 1, Np_in(1)
        
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
	
	DO i = 1, Np_in(1)
	    DO j = -nx/2, nx/2-1
			theta = particles(1,i)%x*j*twopi/Lx
			particles(1,i)%x = particles(1,i)%x + x1(j)*COS(theta+thetax(j))
			particles(1,i)%ux = particles(1,i)%ux + v1(j)*COS(theta+thetav(j))
		END DO
		
		!PRINT*,'Electron ID, x, ux = ',i,particles(1,i)%x, particles(1,i)%ux 
	END DO



END SUBROUTINE add_perturbations


END MODULE init_mod