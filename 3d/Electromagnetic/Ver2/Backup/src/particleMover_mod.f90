MODULE particleMover_mod

USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


! updates velocities and positions of particles (using leap-frog integration) and acculmulates charge at the grid-points
SUBROUTINE move_particles(ts)

    INTEGER, INTENT(IN) :: ts ! timestep index
	INTEGER :: i, j
	REAL*8 :: qm, x, y, z, a(3), b(3), temp, uxx, uyy, uzz, uxp, uyp, uzp
    REAL*8 :: t(3), s(3), gam, tsqr1

		
	! clear kinitic energy and momentum accumulator	
	KE = 0.d0
	Px = 0.d0
	Py = 0.d0
	Pz = 0.d0
    
    ! copy initial positions in buffer array
    pos_buffer(:,:,1) = particles(:,:)%x
    pos_buffer(:,:,2) = particles(:,:)%y
    pos_buffer(:,:,3) = particles(:,:)%z
   
   
    ! loop through all particles and update their positions
    DO i= 1,ns   
	
        IF(i .EQ. 1) THEN
            qm = qm_e
        ELSE IF(i .EQ. 2) THEN 	
            qm = qm_i
        END IF
        
        DO j = 1, N

            ! only update particles that are inside the domain
            IF(.NOT. particles(i,j)%oob) THEN

                x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy
                z = particles(i,j)%z/dz
            
                !******************************
                ! update u(t-dt/2) -> u(t+dt/2)
                !******************************
               ! IF(ts .LT. maxsteps) THEN 
               !     b = 0.d0
               !     a(1) = 0.05*c
               !     a(2) = 0.d0		
               ! ELSE
                    !CALL accel(x, y, z, qm, a, b)
                !END IF
             
             
                ! compute acceleration
                CALL accel(x, y, z, qm, a, b)
             
                !PRINT*,'a = ',a
                !PRINT*,'b = ',b
               
                !a = 0.d0
                !b = 0.d0

                ! half electric push (u^n-1/2 ->u^-) 
                ! half electric push (u^n-1/2 ->u^-) 
                particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*dt
                particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*dt
                particles(i,j)%uz = particles(i,j)%uz + 0.5*a(3)*dt
						
                gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2 + particles(i,j)%uz**2 )/c**2 )
                t     = 0.5d0 * b * dt / (gam*c)
                tsqr1 = 1.d0 + SUM(t**2)
                s     = (2.d0 / tsqr1) * t  
                uxx   = particles(i,j)%ux
                uyy   = particles(i,j)%uy
                uzz   = particles(i,j)%uz
                
                ! magnetic rotation (u^- -> u^+)			
                uxp = uxx +  (uyy*t(3) - uzz*t(2))
                uyp = uyy +  (uzz*t(1) - uxx*t(3))
                uzp = uzz +  (uxx*t(2) - uyy*t(1))
            
                particles(i,j)%ux = uxx +  (uyp*s(3) - uzp*s(2))
                particles(i,j)%uy = uyy +  (uzp*s(1) - uxp*s(3))
                particles(i,j)%uz = uzz +  (uxp*s(2) - uyp*s(1))
			
            
                ! remaining half electric push (u^+ -> u^n+1/2)		
                particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*dt
                particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*dt
                particles(i,j)%uz = particles(i,j)%uz + 0.5*a(3)*dt
                             
                KE = KE + (gam-1.d0)* particles(i,j)%m * (c**2)
                Px = Px + particles(i,j)%m*particles(i,j)%ux
                Py = Py + particles(i,j)%m*particles(i,j)%uy
                Pz = Pz + particles(i,j)%m*particles(i,j)%uz
			
                !***********************
                ! update x(t) -> x(t+dt)
                !***********************
                !gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2 + particles(i,j)%uz**2 )/c**2 )
                particles(i,j)%x = particles(i,j)%x + particles(i,j)%ux*dt/gam
                particles(i,j)%y = particles(i,j)%y + particles(i,j)%uy*dt/gam       
                particles(i,j)%z = particles(i,j)%z + particles(i,j)%uz*dt/gam       
            

            END IF
            
		END DO
	END DO

	
END SUBROUTINE move_particles


! interpolates force from grid to particle location
SUBROUTINE accel(x, y, z, qm, a, b)

	REAL*8, INTENT(IN) :: x, y, z, qm
	REAL*8, INTENT(INOUT) :: a(3), b(3)
	INTEGER :: ix, iy, iz
    REAL*8 :: delx, dely, delz, dxx, dyy, dzz


	IF(interpolation_type .EQ. 1 .OR. interpolation_type .EQ. 3) THEN
		ix = x + 0.5
        iy = y + 0.5
        iz = z + 0.5
        
		a(1) = qm*Ex_grid(ix,iy,iz)
		a(2) = qm*Ey_grid(ix,iy,iz)
		a(3) = qm*Ez_grid(ix,iy,iz)
        
        b(1) = qm*Bx_grid(ix,iy,iz)
        b(2) = qm*By_grid(ix,iy,iz)
        b(3) = qm*Bz_grid(ix,iy,iz)
        
	ELSE IF(interpolation_type .EQ. 2) THEN
	
		ix = x 
        iy = y
        iz = z
        
        delx = x-ix
        dely = y-iy
        delz = z-iz
        dxx = 1.d0 - delx
        dyy = 1.d0 - dely
        dzz = 1.d0 - delz
       
        
		a(1) = qm*( dxx*dyy*dzz  *Ex_grid(ix,iy,iz)     + delx*dyy*dzz  *Ex_grid(ix+1,iy,iz)      + &
                    dxx*dely*dzz *Ex_grid(ix,iy+1,iz)   + dxx*dyy*delz  *Ex_grid(ix,iy,iz+1)      + &
                    dxx*dely*delz*Ex_grid(ix,iy+1,iz+1) + delx*dyy*delz *Ex_grid(ix+1,iy,iz+1)    + &
                    delx*dely*dzz*Ex_grid(ix+1,iy+1,iz) + delx*dely*delz*Ex_grid(ix+1,iy+1,iz+1)  )
                    
		a(2) = qm*( dxx*dyy*dzz  *Ey_grid(ix,iy,iz)     + delx*dyy*dzz  *Ey_grid(ix+1,iy,iz)      + &
                    dxx*dely*dzz *Ey_grid(ix,iy+1,iz)   + dxx*dyy*delz  *Ey_grid(ix,iy,iz+1)      + &
                    dxx*dely*delz*Ey_grid(ix,iy+1,iz+1) + delx*dyy*delz *Ey_grid(ix+1,iy,iz+1)    + &
                    delx*dely*dzz*Ey_grid(ix+1,iy+1,iz) + delx*dely*delz*Ey_grid(ix+1,iy+1,iz+1)  )
                    
        a(3) = qm*( dxx*dyy*dzz  *Ez_grid(ix,iy,iz)     + delx*dyy*dzz  *Ez_grid(ix+1,iy,iz)      + &
                    dxx*dely*dzz *Ez_grid(ix,iy+1,iz)   + dxx*dyy*delz  *Ez_grid(ix,iy,iz+1)      + &
                    dxx*dely*delz*Ez_grid(ix,iy+1,iz+1) + delx*dyy*delz *Ez_grid(ix+1,iy,iz+1)    + &
                    delx*dely*dzz*Ez_grid(ix+1,iy+1,iz) + delx*dely*delz*Ez_grid(ix+1,iy+1,iz+1)  )
                   
        b(1) = qm*( dxx*dyy*dzz  *Bx_grid(ix,iy,iz)     + delx*dyy*dzz  *Bx_grid(ix+1,iy,iz)      + &
                    dxx*dely*dzz *Bx_grid(ix,iy+1,iz)   + dxx*dyy*delz  *Bx_grid(ix,iy,iz+1)      + &
                    dxx*dely*delz*Bx_grid(ix,iy+1,iz+1) + delx*dyy*delz *Bx_grid(ix+1,iy,iz+1)    + &
                    delx*dely*dzz*Bx_grid(ix+1,iy+1,iz) + delx*dely*delz*Bx_grid(ix+1,iy+1,iz+1)  )
                    
		b(2) = qm*( dxx*dyy*dzz  *By_grid(ix,iy,iz)     + delx*dyy*dzz  *By_grid(ix+1,iy,iz)      + &
                    dxx*dely*dzz *By_grid(ix,iy+1,iz)   + dxx*dyy*delz  *By_grid(ix,iy,iz+1)      + &
                    dxx*dely*delz*By_grid(ix,iy+1,iz+1) + delx*dyy*delz *By_grid(ix+1,iy,iz+1)    + &
                    delx*dely*dzz*By_grid(ix+1,iy+1,iz) + delx*dely*delz*By_grid(ix+1,iy+1,iz+1)  )
                    
        b(3) = qm*( dxx*dyy*dzz  *Bz_grid(ix,iy,iz)     + delx*dyy*dzz  *Bz_grid(ix+1,iy,iz)      + &
                    dxx*dely*dzz *Bz_grid(ix,iy+1,iz)   + dxx*dyy*delz  *Bz_grid(ix,iy,iz+1)      + &
                    dxx*dely*delz*Bz_grid(ix,iy+1,iz+1) + delx*dyy*delz *Bz_grid(ix+1,iy,iz+1)    + &
                    delx*dely*dzz*Bz_grid(ix+1,iy+1,iz) + delx*dely*delz*Bz_grid(ix+1,iy+1,iz+1)  )
	END IF


END SUBROUTINE accel



! conservative current deposition using Umeda (2003) Zig-zag algorithm
SUBROUTINE deposit_currents()

    REAL*8 :: x1, y1, z1, x2, y2, z2
    REAL*8 :: qdt
    REAL*8 :: xr, yr, zr    
    REAL*8 :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2
    REAL*8 :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2    
    INTEGER :: sp, np
    INTEGER :: i1, i2, j1, j2, k1, k2
    REAL*8  :: idx, idy, idz, idxyz         


    idx = 1.d0 / dx 
    idy = 1.d0 / dy 
    idz = 1.d0 / dz 
    idxyz = 1.d0 / (dx*dy*dz) 

    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0

    ! loop over all particles
    DO sp = 1, ns   
    
        qdt = particles(sp,1)%q / dt
        
        DO np = 1, N

            ! only contributins from particles that are inside the computational domain
            IF(.NOT. particles(sp,np)%oob) THEN

                ! old position
                x1 = pos_buffer(sp,np,1)
                y1 = pos_buffer(sp,np,2)
                z1 = pos_buffer(sp,np,3)

                i1 = x1 * idx
                j1 = y1 * idy
                k1 = z1 * idz
    
                ! new position
                x2 = particles(sp,np)%x
                y2 = particles(sp,np)%y
                z2 = particles(sp,np)%z
    
                i2 = x2 * idx
                j2 = y2 * idy
                k2 = z2 * idz
    
    
                ! relay point
                
                ! with IF statement
                IF(i1 .EQ. i2) THEN
                    xr = 0.5*(x1+x2)
                ELSE
                    xr = MAX(DBLE(i1*dx), DBLE(i2*dx))
                END IF    
                IF(j1 .EQ. j2) THEN
                    yr = 0.5*(y1+y2)
                ELSE
                    yr = MAX(DBLE(j1*dy), DBLE(j2*dy))
                END IF    
                IF(k1 .EQ. k2) THEN
                    zr = 0.5*(z1+z2)
                ELSE
                    zr = MAX(DBLE(k1*dz), DBLE(k2*dz))
                END IF    
               
                
                ! without IF statement 
                !xr = MIN( MIN(DBLE(i1*dx), DBLE(i2*dx)) + dx, MAX( MAX(DBLE(i1*dx), DBLE(i2*dx)), 0.5*(x1+x2) ) )
                !yr = MIN( MIN(DBLE(j1*dy), DBLE(j2*dy)) + dy, MAX( MAX(DBLE(j1*dy), DBLE(j2*dy)), 0.5*(y1+y2) ) )
                !zr = MIN( MIN(DBLE(k1*dz), DBLE(k2*dz)) + dz, MAX( MAX(DBLE(k1*dz), DBLE(k2*dz)), 0.5*(z1+z2) ) )
    
                ! compute segment flux and shape-factors                
                Fx1 = qdt * (xr - x1)
                Fy1 = qdt * (yr - y1)
                Fz1 = qdt * (zr - z1)
    
                Fx2 = qdt * (x2 - xr)
                Fy2 = qdt * (y2 - yr)
                Fz2 = qdt * (z2 - zr)
    
                Wx1 = 0.5 * idx * (x1 + xr) - i1
                Wy1 = 0.5 * idy * (y1 + yr) - j1
                Wz1 = 0.5 * idz * (z1 + zr) - k1

                Wx2 = 0.5 * idx * (xr + x2) - i2
                Wy2 = 0.5 * idy * (yr + y2) - j2
                Wz2 = 0.5 * idz * (zr + z2) - k2


               ! compute fluxes at all 24 adjacent grid point 
               Jx(i1,j1,k1)     = Jx(i1,j1,k1)  + idxyz * Fx1 * (1.d0 - Wy1) * (1.d0 - Wz1)
               Jx(i1,j1+1,k1)   = Jx(i1,j1+1,k1) + idxyz * Fx1 * Wy1 * (1.d0 - Wz1)
               Jx(i1,j1,k1+1)   = Jx(i1,j1,k1+1)  + idxyz * Fx1 * (1.d0 - Wy1) * Wz1
               Jx(i1,j1+1,k1+1) = Jx(i1,j1+1,k1+1) + idxyz * Fx1 * Wy1 * Wz1

               Jx(i2,j2,k2)     = Jx(i2,j2,k2)  + idxyz * Fx2 * (1.d0 - Wy2) * (1.d0 - Wz2)
               Jx(i2,j2+1,k2)   = Jx(i2,j2+1,k2) + idxyz * Fx2 * Wy2 * (1.d0 - Wz2)
               Jx(i2,j2,k2+1)   = Jx(i2,j2,k2+1)  + idxyz * Fx2 * (1.d0 - Wy2) * Wz2
               Jx(i2,j2+1,k2+1) = Jx(i2,j2+1,k2+1) + idxyz * Fx2 * Wy2 * Wz2

               Jy(i1,j1,k1)     = Jy(i1,j1,k1) + idxyz * Fy1 * (1.d0 - Wx1) * (1.d0 - Wz1)
               Jy(i1+1,j1,k1)   = Jy(i1+1,j1,k1) + idxyz * Fy1 * Wx1 * (1.d0 - Wz1)
               Jy(i1,j1,k1+1)   = Jy(i1,j1,k1+1)  + idxyz * Fy1 * (1.d0 - Wx1) * Wz1
               Jy(i1+1,j1,k1+1) = Jy(i1+1,j1,k1+1) + idxyz * Fy1 * Wx1 * Wz1

               Jy(i2,j2,k2)     = Jy(i2,j2,k2) + idxyz * Fy2 * (1.d0 - Wx2) * (1.d0 - Wz2)
               Jy(i2+1,j2,k2)   = Jy(i2+1,j2,k2) + idxyz * Fy2 * Wx2 * (1.d0 - Wz2)
               Jy(i2,j2,k2+1)   = Jy(i2,j2,k2+1)  + idxyz * Fy2 * (1.d0 - Wx2) * Wz2
               Jy(i2+1,j2,k2+1) = Jy(i2,j2+1,k2+1) + idxyz * Fy2 * Wx2 * Wz2

               Jz(i1,j1,k1)     = Jz(i1,j1,k1) + idxyz * Fz1 * (1.d0 - Wx1) * (1.d0 - Wy1)
               Jz(i1+1,j1,k1)   = Jz(i1+1,j1,k1) + idxyz * Fz1 * Wx1 * (1.d0 - Wy1)
               Jz(i1,j1+1,k1)   = Jz(i1,j1+1,k1)  + idxyz * Fz1 * (1.d0 - Wx1) * Wy1
               Jz(i1+1,j1+1,k1) = Jz(i1+1,j1+1,k1) + idxyz * Fz1 * Wx1 * Wy1

               Jz(i2,j2,k2)     = Jz(i2,j2,k2) + idxyz * Fz2 * (1.d0 - Wx2) * (1.d0 - Wy2)
               Jz(i2+1,j2,k2)   = Jz(i2+1,j2,k2) + idxyz * Fz2 * Wx2 * (1.d0 - Wy2)
               Jz(i2,j2+1,k2)   = Jz(i2,j2+1,k2)  + idxyz * Fz2 * (1.d0 - Wx2) * Wy2
               Jz(i2+1,j2+1,k2) = Jz(i2+1,j2+1,k2) + idxyz * Fz2 * Wx2 * Wy2
               
            END IF
            
        END DO        
    END DO
    
    
    ! If periodic boundaries, need to touch up the grid deposited charges and currents at boundaries
    IF (bndry .EQ. 1) THEN

        buffer = Jx
        Jx(0,:,:) = Jx(0,:,:) + Jx(nx,:,:)	
        Jx(1,:,:) = Jx(1,:,:) + Jx(nx+1,:,:)	
        Jx(nx,:,:) = Jx(nx,:,:) + buffer(0,:,:)
        Jx(nx-1,:,:) = Jx(nx-1,:,:) + Jx(-1,:,:)

        buffer = Jy
        Jy(:,0,:) = Jy(:,0,:) + Jy(:,ny,:)	
        Jy(:,1,:) = Jy(:,1,:) + Jy(:,ny+1,:)	
        Jy(:,ny,:) = Jy(:,ny,:) + buffer(:,0,:)
        Jy(:,ny-1,:) = Jy(:,ny-1,:) + Jy(:,-1,:)
        
        buffer = Jz
        Jz(:,:,0) = Jz(:,:,0) + Jz(:,:,nz)	
        Jz(:,:,1) = Jz(:,:,1) + Jz(:,:,nz+1)	
        Jz(:,:,nz) = Jz(:,:,nz) + buffer(:,:,0)
        Jz(:,:,nz-1) = Jz(:,:,nz-1) + Jz(:,:,-1)
        
    END IF
    
    
    ! Filtering the current makes the spatial distriution smoother and attenuates short wavelength (high k) components of the fourier transform.
    ! Note: Filtering is absolutely crucial. The electromagnetic waves generated by our finite difference equations are dispersive, the wave speed
    ! decreases monotonically for higher k. This results in numerical Cerenkov radiation from particles moving at relativistic speeds. 
    ! (e.g. try moving a single point charge at relativistic speed. It will generate what looks like a bow shock leaving a triangular region
    ! in it's wake containing the outward going e.m. waves. Of course, the particle is moving faster than those waves hence the triangular wake..)
    ! Filtering will attentuate those (aliased) waves at higher frequencies and so mitagate this unphysical effect to a significant extent.
    IF(current_filter_on)THEN
        CALL filter(Jx)
        CALL filter(Jy)
        CALL filter(Jz)
    END IF
    

END SUBROUTINE deposit_currents


! 1D particle charge shape function, first order (i.e. piecewise linear)
FUNCTION S_1D(x) RESULT (sx)

    REAL*8, INTENT(IN) :: x
    REAL*8 :: sx

    IF(ABS(x) .LT. 1.d0) THEN
        sx = 1.d0 - ABS(x)     
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S_1D


! deposits charge to grid
SUBROUTINE deposit_charges()

	INTEGER :: i, j, ix, iy, iz
    REAL*8 :: delx, dely, delz, dxx, dyy, dzz
	REAL*8 :: x, y, z, qdxyz

    ! clear density
    rho = 0.d0
    !rho(1:nx,1:ny,1:nz) = rho0   -> uniform background 
	ne = 0.d0
	ni = 0.d0
	

    ! loop over all particles
    DO i= 1,ns   
    
        qdxyz = particles(i,1)%q/(dx*dy*dz)
        
        DO j = 1, N

            ! only contributins from particles that are inside the computational domain
            IF(.NOT. particles(i,j)%oob) THEN

                x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy
                z = particles(i,j)%z/dz

                !**********************************
                ! deposit charge at grid points
                ! rho_j,k
                !**********************************
			
                IF(interpolation_type .EQ. 1) THEN			             
				
                    !compute index of occupied cell
                    ix = x + 0.5  
                    iy = y + 0.5
                    iz = z + 0.5
                

                    ! accumulate charge in the cell occupied by the particle (i.e. nearest grid point)
                    rho(ix,iy,iz) = rho(ix,iy,iz) + qdxyz
								
                    IF(i .EQ. 1) ne(ix,iy,iz) = ne(ix,iy,iz) + 1.0				
                    IF(i .EQ. 2) ni(ix,iy,iz) = ni(ix,iy,iz) + 1.0				
				
                ELSE IF(interpolation_type .EQ. 2 .OR. interpolation_type .EQ. 3) THEN

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
        
                    ! accumulate charge in cells occupied by the cloud weighted according to the fraction of the cell occupied by the cloud, i.e. linear interpolation
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
            
                END IF
            END IF
        END DO
    END DO


    ! If periodic boundaries, need to touch up the grid deposited charges at boundaries
    IF (bndry .EQ. 1) THEN
        
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
    END IF         


END SUBROUTINE deposit_charges



! 3D biomial filter to smooth out the current so as to elimoinate high-frequency spatial fourier components.
! Removing high frequency components helps mitigate spurious numerical Cherenkov radiaiton (which are generated
! due to the fact that electromagnetic wave solutions in the finite differenced approximation have some numerical
! dispersion, with lower phase velocities at higher frequencies.)
SUBROUTINE filter(f)

    REAL*8, INTENT(INOUT) ::  f(-1:nx+3,-1:ny+3,-1:nz+3)
    INTEGER :: i, j, k

    
    ! Successively apply 1D binomial filters separately in x, y, and z directions 
    ! e.g. 1D_x_Filter[ f_i,j,k ] = 0.25*f_i-1,j,k + 0.5*f_i,j,k + 0.25*f_i+1,j,k 
    !      1D_y_Filter[ f_i,j,k ] = 0.25*f_i,j-1,k + 0.5*f_i,j,k + 0.25*f_i,j+1,k 
    !      1D_z_Filter[ f_i,j,k ] = 0.25*f_i,j,k-1 + 0.5*f_i,j,k + 0.25*f_i,j,k+1 
    
       
    ! x-filter pass
    buffer(-1:nx+1,-1:ny+1, -1:nz+1) = f(-1:nx+1,-1:ny+1, -1:nz+1)

    DO k = 0, nz
        DO j = 0, ny
            DO i = 0, nx
                f(i,j,k) = 0.25*buffer(i-1,j,k) + 0.5*buffer(i,j,k) + 0.25*buffer(i+1,j,k)
            END DO
        END DO        
    END DO    
	
    ! y-filter pass
    buffer(-1:nx+1,-1:ny+1, -1:nz+1) = f(-1:nx+1,-1:ny+1, -1:nz+1)

    DO k = 0, nz
        DO j = 0, ny
            DO i = 0, nx
                f(i,j,k) = 0.25*buffer(i,j-1,k) + 0.5*buffer(i,j,k) + 0.25*buffer(i,j+1,k)
            END DO
        END DO        
    END DO    
	
	! z-filter pass
    buffer(-1:nx+1,-1:ny+1, -1:nz+1) = f(-1:nx+1,-1:ny+1, -1:nz+1)

    DO k = 0, nz
        DO j = 0, ny
            DO i = 0, nx
                f(i,j,k) = 0.25*buffer(i,j,k-1) + 0.5*buffer(i,j,k) + 0.25*buffer(i,j,k+1)
            END DO
        END DO        
    END DO   
    
    
END SUBROUTINE filter



! enforces particle boundary conditions
SUBROUTINE particle_boundary()

    INTEGER :: i, j
    
    
    ! loop over all particles
    DO i= 1,ns           
        DO j = 1, N

            ! periodic boundary conditions
            IF(bndry .eq. 1) THEN
    
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
                
            END IF

            ! outflow boundary conditions (in this case, need to distribute ions and electrons uniformly across boundary
            ! to avoid net charge buildup)
            IF(bndry .EQ. 2) THEN    
                ! set out of bounds flag
                IF(particles(i,j)%x .LT. xmin .OR. particles(i,j)%x .GT. xmax .OR. &
                   particles(i,j)%y .LT. ymin .OR. particles(i,j)%y .GT. ymax .OR. &
                   particles(i,j)%z .LT. zmin .OR. particles(i,j)%z .GT. zmax) THEN

                    particles(i,j)%oob = .TRUE.   
                END IF
     
            END IF
    
        END DO
    END DO

END SUBROUTINE particle_boundary



END MODULE particleMover_mod