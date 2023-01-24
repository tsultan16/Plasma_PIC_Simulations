MODULE particleMover_mod

USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


! updates velocities and positions of particles (using leap-frog integration) and acculmulates charge at the grid-points
SUBROUTINE move_particles(ts)

    INTEGER, INTENT(IN) :: ts ! timestep index
	INTEGER :: i, j
	REAL*8 :: qm, x, y, z, a(3), b(3), uxx, uyy, uzz, uxp, uyp, uzp
    REAL*8 :: gam, igam, f

		
	! clear kinitic energy and momentum accumulator	
	K_E = 0.d0
	Px = 0.d0
	Py = 0.d0
	Pz = 0.d0
    
    !PRINT*,''
    !PRINT*,'Np_in = ',Np_in
    !PRINT*,''
       
    ! loop through all particles and update their positions
    DO i= 1,ns   

        qm = q(i)/m(i)
        
        !$OMP PARALLEL DO REDUCTION(+:KE,Px,Py,Pz) &
        !$OMP PRIVATE(x, y, z, a, b, gam, igam, f, uxx, uyy, uzz, uxp, uyp, uzp)
        DO j = 1, Np_in(i)

            !PRINT*,'i,j = ',i,j
 
            ! only update particles that are inside the domain
            !IF(.NOT. particles(i,j)%oob) THEN

            x = particles(i,j)%x
            y = particles(i,j)%y
            z = particles(i,j)%z
                       

                !******************************
                ! update u(t-dt/2) -> u(t+dt/2)
                !******************************
                               
                ! compute acceleration
                CALL accel(x, y, z, qm, a, b)
             
                !PRINT*,'a = ',a
                !PRINT*,'b = ',b
                !a = 0.d0
                !b = 0.d0
     
                gam = c / SQRT( c**2 - (particles(i,j)%vx**2 + particles(i,j)%vy**2 + particles(i,j)%vz**2 ) )
     
                ! half electric push (u^n-1/2 ->u^-) 
                uxx = gam * particles(i,j)%vx + 0.5d0 * a(1)
                uyy = gam * particles(i,j)%vy + 0.5d0 * a(2)
                !IF(nz .GT. 1) THEN 
                    uzz = gam * particles(i,j)%vz + 0.5d0 * a(3)
                !ELSE
                !    uzz = gam * particles(i,j)%vz
                !END IF
                
                
                !IF(i .EQ. 1) particles(i,j)%vx = 0.3

                igam = c / SQRT( c**2 + (uxx**2 + uyy**2 + uzz**2 ) )			
                b = 0.5d0 * igam*b/c
                f=2.d0 / (1.d0+SUM(b**2))

                ! magnetic rotation (u^- -> u^+)  +  remaining half electric push (u^+ -> u^n+1/2)			
                uxp = (uxx + uyy*b(3) - uzz*b(2))*f
                uyp = (uyy + uzz*b(1) - uxx*b(3))*f
                uzp = (uzz + uxx*b(2) - uyy*b(1))*f
            
                uxx = uxx + uyp*b(3) - uzp*b(2) + 0.5d0 * a(1)
                uyy = uyy + uzp*b(1) - uxp*b(3) + 0.5d0 * a(2)
                !IF(nz .GT. 1) THEN  
                    uzz = uzz + uxp*b(2) - uyp*b(1) + 0.5d0 * a(3)
                !ELSE
                !    uzz = uzz + uxp*b(2) - uyp*b(1) 
                !END IF
            
                igam = c / SQRT( c**2 + (uxx**2 + uyy**2 + uzz**2 ) )			                
                particles(i,j)%vx = igam * uxx
                particles(i,j)%vy = igam * uyy
                particles(i,j)%vz = igam * uzz
                             

                ! ckeck for unsphysical velocity
                IF(particles(i,j)%vx**2+ particles(i,j)%vy**2 + particles(i,j)%vz**2 .GE. c**2)THEN
                    PRINT*,'ERROR. Speed of light exceeded in particle mover...vx,vy,vz =', &
                             particles(i,j)%vx, particles(i,j)%vy, particles(i,j)%vz
                    PRINT*,'x,y,z, =',x,y,z         
                    STOP        
                END IF
                
                
                K_E = K_E + (gam-1.d0)* m(i) * (c**2)
                !Px = Px + m(i)*particles(i,j)%vx
                !Py = Py + m(i)*particles(i,j)%vy
                !Pz = Pz + m(i)*particles(i,j)%vz
			
            
                !***********************
                ! update x(t) -> x(t+dt)
                !***********************
                
                particles(i,j)%x = particles(i,j)%x + particles(i,j)%vx 
                particles(i,j)%y = particles(i,j)%y + particles(i,j)%vy       
                IF(nz .GT. 1) particles(i,j)%z = particles(i,j)%z + particles(i,j)%vz    ! For 2.5d (i.e. nz=1),  don't update z-position  
 
                
                !IF(i .EQ. 1) THEN
                !    particles(i,j)%vz =  (5.d0*twopi/100.d0)*COS(twopi*ts/100.d0)
                !    particles(i,j)%z =    (nz-8) + 5.d0*SIN(twopi*ts/100.d0)
                !ELSE IF(i .EQ. 2) THEN
                !    particles(i,j)%vz = -(5.d0*twopi/200.d0)*COS(twopi*ts/100.d0)
                !    particles(i,j)%z = (nz-8) - 5.d0*SIN(twopi*ts/100.d0)
                !END IF
            
		END DO
        !$OMP END PARALLEL DO

	END DO

	
END SUBROUTINE move_particles


! interpolates force from grid to particle location
SUBROUTINE accel(x, y, z, qm, a, b)

	REAL*8, INTENT(IN) :: x, y, z, qm
	REAL*8, INTENT(INOUT) :: a(3), b(3)
	INTEGER :: ix, iy, iz, i, j, k
    REAL*8 :: Sint(-1:2,3)

    
    a = 0.d0
    b = 0.d0
    Sint = 0.d0


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
	
        ! first order force interpolation will involve two grid-points, one on either side of the particle, 
        ! on each dimension of space
       
		ix = x 
        iy = y
        iz = z       
        
        ! compute the 1st order interpolation co-efficients
        
        DO i = 0,1
            Sint(i,1) = S1_1D((ix+i)-x)
            Sint(i,2) = S1_1D((iy+i)-y)
            Sint(i,3) = S1_1D((iz+i)-z)
        END DO       
       
        ! now interpolate from grid to particle location
        DO k = 0,1
            DO j = 0,1
                DO i = 0,1
                    a(1) = a(1) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Ex_grid(ix+i,iy+j,iz+k)) 
                    a(2) = a(2) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Ey_grid(ix+i,iy+j,iz+k)) 
                    a(3) = a(3) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Ez_grid(ix+i,iy+j,iz+k)) 
                    
                    b(1) = b(1) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Bx_grid(ix+i,iy+j,iz+k)) 
                    b(2) = b(2) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*By_grid(ix+i,iy+j,iz+k)) 
                    b(3) = b(3) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Bz_grid(ix+i,iy+j,iz+k)) 
                    
                END DO
            END DO
        END DO        
                    
    ELSE IF(interpolation_type .EQ. 3) THEN
	  
        ! second order force interpolation will involve four grid-points , two on either side of the particle, 
        ! on each dimension of space)
       
		ix = x 
        iy = y
        iz = z       
        
        
        !###### WARNING: If ix ( or iy or iz) = 1 or nx (ny, nz) , then the outermost grid values won't be available.
        ! Should probably switch to first order interpolation in that case..
        
        
        ! compute the 2nd order interpolation co-efficients
        DO i = -1,2
            Sint(i,1) = S2_1D((ix+i)-x)
            Sint(i,2) = S2_1D((iy+i)-y)
            Sint(i,3) = S2_1D((iz+i)-z)
        END DO       
        
       
        ! now interpolate from grid to particle location
        DO k = -1,2
            DO j = -1,2
                DO i = -1,2
                    a(1) = a(1) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Ex_grid(ix+i,iy+j,iz+k)) 
                    a(2) = a(2) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Ey_grid(ix+i,iy+j,iz+k)) 
                    a(3) = a(3) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Ez_grid(ix+i,iy+j,iz+k)) 
                    
                    b(1) = b(1) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Bx_grid(ix+i,iy+j,iz+k)) 
                    b(2) = b(2) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*By_grid(ix+i,iy+j,iz+k)) 
                    b(3) = b(3) + qm * ( Sint(i,1)*Sint(j,2)*Sint(k,3)*Bz_grid(ix+i,iy+j,iz+k)) 
                    
                END DO
            END DO
        END DO      
       
       
	END IF


END SUBROUTINE accel



! 1D particle 1st order shape function (piecewise linear)
FUNCTION S1_1D(x) RESULT (sx)

    REAL*8, INTENT(IN) :: x
    REAL*8 :: sx

    IF(ABS(x) .LT. 1.d0) THEN
        sx = 1.d0 - ABS(x)     
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S1_1D


! 1D particle 2nd order shape function (quadratic spline)
FUNCTION S2_1D(x) RESULT (sx)

    REAL*8, INTENT(IN) :: x
    REAL*8 :: sx

    IF(ABS(x) .LT. 0.5d0) THEN
        sx = (3.d0/4.d0) - x**2     
    ELSE IF(ABS(x) .GE. 0.5d0 .AND. ABS(x) .LE. 1.5d0)THEN
        sx = 0.5d0 * ( 1.5d0 - ABS(x))**2
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S2_1D



! conservative current deposition using Umeda (2003) Zig-zag algorithm
SUBROUTINE deposit_currents_zigzag()

    REAL*8 :: x1, y1, z1, x2, y2, z2
    REAL*8 :: qdt
    REAL*8 :: xr, yr, zr    
    REAL*8 :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2
    REAL*8 :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2    
    INTEGER :: sp, np
    INTEGER :: i1, i2, j1, j2, k1, k2
    REAL*8  :: idx, idy, idz, idxyz, J_max         


    J_max = 0.d0

    idx = 1.d0 
    idy = 1.d0 
    idz = 1.d0  
    idxyz = 1.d0 

    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0

    
    ! loop over all particles
    DO sp = 1, ns   
    
        qdt = q(sp) 
        
        !PRINT*,'sp, qdt = ', sp, qdt
        
        DO np = 1, Np_in(sp)

            !PRINT*,'np, Np_in', np, Np_in

            ! only contributins from particles that are inside the computational domain

                ! old position
                x1 = particles(sp,np)%x - particles(sp,np)%vx
                y1 = particles(sp,np)%y - particles(sp,np)%vy
                z1 = particles(sp,np)%z - particles(sp,np)%vz

                i1 = x1 
                j1 = y1
                k1 = z1 
    
                ! new position
                x2 = particles(sp,np)%x
                y2 = particles(sp,np)%y
                z2 = particles(sp,np)%z
    
                i2 = x2 
                j2 = y2 
                k2 = z2 
 
    
                ! relay point
                
                ! with IF statement
                !IF(i1 .EQ. i2) THEN
                !    xr = 0.5*(x1+x2)
                !ELSE
                !    xr = MAX(DBLE(i1*dx), DBLE(i2*dx))
                !END IF    
                !IF(j1 .EQ. j2) THEN
                !    yr = 0.5*(y1+y2)
                !ELSE
                !    yr = MAX(DBLE(j1*dy), DBLE(j2*dy))
                !END IF    
                !IF(k1 .EQ. k2) THEN
                !    zr = 0.5*(z1+z2)
                !ELSE
                !    zr = MAX(DBLE(k1*dz), DBLE(k2*dz))
                !END IF    
               
                
                ! without IF statement 
                xr = MIN( MIN(DBLE(i1), DBLE(i2)) + 1.d0, MAX( MAX(DBLE(i1), DBLE(i2)), 0.5d0*(x1+x2) ) )
                yr = MIN( MIN(DBLE(j1), DBLE(j2)) + 1.d0, MAX( MAX(DBLE(j1), DBLE(j2)), 0.5d0*(y1+y2) ) )
                zr = MIN( MIN(DBLE(k1), DBLE(k2)) + 1.d0, MAX( MAX(DBLE(k1), DBLE(k2)), 0.5d0*(z1+z2) ) )
    
    
                ! compute segment flux and shape-factors                
                Fx1 = qdt * (xr - x1)
                Fy1 = qdt * (yr - y1)
                Fz1 = qdt * (zr - z1)
    
                Fx2 = qdt * (x2 - xr)
                Fy2 = qdt * (y2 - yr)
                Fz2 = qdt * (z2 - zr)
    
                Wx1 = 0.5d0 * (x1 + xr) - i1
                Wy1 = 0.5d0 * (y1 + yr) - j1
                Wz1 = 0.5d0 * (z1 + zr) - k1

                Wx2 = 0.5d0 * (xr + x2) - i2
                Wy2 = 0.5d0 * (yr + y2) - j2
                Wz2 = 0.5d0 * (zr + z2) - k2


               ! compute fluxes at all 24 adjacent grid point 
               Jx(i1,j1,k1)     = Jx(i1,j1,k1)  + Fx1 * (1.d0 - Wy1) * (1.d0 - Wz1)
               Jx(i1,j1+1,k1)   = Jx(i1,j1+1,k1) + Fx1 * Wy1 * (1.d0 - Wz1)
               Jx(i1,j1,k1+1)   = Jx(i1,j1,k1+1)  + Fx1 * (1.d0 - Wy1) * Wz1
               Jx(i1,j1+1,k1+1) = Jx(i1,j1+1,k1+1) + Fx1 * Wy1 * Wz1

               Jx(i2,j2,k2)     = Jx(i2,j2,k2)  + Fx2 * (1.d0 - Wy2) * (1.d0 - Wz2)
               Jx(i2,j2+1,k2)   = Jx(i2,j2+1,k2) + Fx2 * Wy2 * (1.d0 - Wz2)
               Jx(i2,j2,k2+1)   = Jx(i2,j2,k2+1)  + Fx2 * (1.d0 - Wy2) * Wz2
               Jx(i2,j2+1,k2+1) = Jx(i2,j2+1,k2+1) +  Fx2 * Wy2 * Wz2

               Jy(i1,j1,k1)     = Jy(i1,j1,k1) + Fy1 * (1.d0 - Wx1) * (1.d0 - Wz1)
               Jy(i1+1,j1,k1)   = Jy(i1+1,j1,k1) + Fy1 * Wx1 * (1.d0 - Wz1)
               Jy(i1,j1,k1+1)   = Jy(i1,j1,k1+1)  + Fy1 * (1.d0 - Wx1) * Wz1
               Jy(i1+1,j1,k1+1) = Jy(i1+1,j1,k1+1) + Fy1 * Wx1 * Wz1

               Jy(i2,j2,k2)     = Jy(i2,j2,k2) + Fy2 * (1.d0 - Wx2) * (1.d0 - Wz2)
               Jy(i2+1,j2,k2)   = Jy(i2+1,j2,k2) + Fy2 * Wx2 * (1.d0 - Wz2)
               Jy(i2,j2,k2+1)   = Jy(i2,j2,k2+1)  + Fy2 * (1.d0 - Wx2) * Wz2
               Jy(i2+1,j2,k2+1) = Jy(i2,j2+1,k2+1) + Fy2 * Wx2 * Wz2

               Jz(i1,j1,k1)     = Jz(i1,j1,k1) + Fz1 * (1.d0 - Wx1) * (1.d0 - Wy1)
               Jz(i1+1,j1,k1)   = Jz(i1+1,j1,k1) + Fz1 * Wx1 * (1.d0 - Wy1)
               Jz(i1,j1+1,k1)   = Jz(i1,j1+1,k1)  + Fz1 * (1.d0 - Wx1) * Wy1
               Jz(i1+1,j1+1,k1) = Jz(i1+1,j1+1,k1) + Fz1 * Wx1 * Wy1

               Jz(i2,j2,k2)     = Jz(i2,j2,k2) + Fz2 * (1.d0 - Wx2) * (1.d0 - Wy2)
               Jz(i2+1,j2,k2)   = Jz(i2+1,j2,k2) + Fz2 * Wx2 * (1.d0 - Wy2)
               Jz(i2,j2+1,k2)   = Jz(i2,j2+1,k2)  + Fz2 * (1.d0 - Wx2) * Wy2
               Jz(i2+1,j2+1,k2) = Jz(i2+1,j2+1,k2) + Fz2 * Wx2 * Wy2
               
            
        END DO        
    END DO
  
  
    IF (bndry .EQ. 1) THEN
     
        IF(nranks_x .EQ. 1) THEN  
       
            buffer = Jx  
            Jx(0,:,:) = Jx(0,:,:) + Jx(nx,:,:)	
            Jx(1,:,:) = Jx(1,:,:) + Jx(nx+1,:,:)	
            Jx(nx,:,:) = Jx(nx,:,:) + buffer(0,:,:)
            Jx(nx-1,:,:) = Jx(nx-1,:,:) + Jx(-1,:,:)
        END IF
        
        IF(nranks_y .EQ. 1) THEN  
       
            buffer = Jy
            Jy(:,0,:) = Jy(:,0,:) + Jy(:,ny,:)	
            Jy(:,1,:) = Jy(:,1,:) + Jy(:,ny+1,:)	
            Jy(:,ny,:) = Jy(:,ny,:) + buffer(:,0,:)
            Jy(:,ny-1,:) = Jy(:,ny-1,:) + Jy(:,-1,:)
        
        END IF
        
        
        buffer = Jz

        Jz(:,:,0) = Jz(:,:,0) + Jz(:,:,nz)	
        Jz(:,:,1) = Jz(:,:,1) + Jz(:,:,nz+1)	

        Jz(:,:,nz) = Jz(:,:,nz) + buffer(:,:,0)
        Jz(:,:,nz-1) = Jz(:,:,nz-1) + Jz(:,:,-1)
        
    END IF
    
    

END SUBROUTINE deposit_currents_zigzag



! conservative current deposition using Esirkepov (2001) algorithm
SUBROUTINE deposit_currents_esirkepov()

    REAL*8 :: x0, y0, z0, x1, y1, z1
    INTEGER :: ix0, iy0, iz0, i, j, k
    REAL*8  :: qdt
    REAL*8  :: S0(-2:2,3), S1(-2:2,3), DS(-2:2,3)
    REAL*8  :: W(-2:2, -2:2, -2:2, 3), temp(-2:2,-2:2,-2:2,3)
    INTEGER :: sp, np



    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0
    
    
    ! clear auxilliary arrays
    S0 = 0.d0
    S1 = 0.d0
    DS = 0.d0
    W = 0.d0
    temp = 0.d0


    ! loop over all particles
    DO sp = 1, ns   
    
        qdt = q(sp) 
        
        DO np = 1, Np_in(sp)

                ! old position
                x0 = particles(sp,np)%x - particles(sp,np)%vx
                y0 = particles(sp,np)%y - particles(sp,np)%vy
                z0 = particles(sp,np)%z - particles(sp,np)%vz


                ! new position
                x1 = particles(sp,np)%x
                y1 = particles(sp,np)%y
                z1 = particles(sp,np)%z
    
         
                ! check for bad step size
                IF(ABS(x1-x0) .GE. c .OR. ABS(y1-y0) .GE. c .OR. ABS(y1-y0) .GE. c) THEN
                    PRINT*,'ERROR. Courant condition violated by particle...'
                    STOP
                END IF
                
    
                ! compute index of occupied cell at initial position
                ix0 = x0 + 0.5  
                iy0 = y0 + 0.5
                iz0 = z0 + 0.5
        
                ! apply protection (otherwise risk of out of bounds J-array access)
                ix0 = MIN( MAX(ix0, 1), nx+1 )
                iy0 = MIN( MAX(iy0, 1), ny+1 )
                iz0 = MIN( MAX(iz0, 1), nz+1 )
                

                ! prepare array of 1D form factors corresponding to initial and final positions and their difference 
                DO i = -2, 2
    
                 	IF(interpolation_type .EQ. 2) THEN

                        S0(i,1) = S1_1D(DBLE(ix0+i)-x0)
                        S0(i,2) = S1_1D(DBLE(iy0+i)-y0)
                        S0(i,3) = S1_1D(DBLE(iz0+i)-z0)
        
                        S1(i,1) = S1_1D(DBLE(ix0+i)-x1)
                        S1(i,2) = S1_1D(DBLE(iy0+i)-y1)
                        S1(i,3) = S1_1D(DBLE(iz0+i)-z1)
                        
                     ELSE IF(interpolation_type .EQ. 3) THEN  
                        
                        S0(i,1) = S2_1D(DBLE(ix0+i)-x0)
                        S0(i,2) = S2_1D(DBLE(iy0+i)-y0)
                        S0(i,3) = S2_1D(DBLE(iz0+i)-z0)
        
                        S1(i,1) = S2_1D(DBLE(ix0+i)-x1)
                        S1(i,2) = S2_1D(DBLE(iy0+i)-y1)
                        S1(i,3) = S2_1D(DBLE(iz0+i)-z1)
                        
                    END IF    
        
                    DS(i,1) = S1(i,1) - S0(i,1)
                    DS(i,2) = S1(i,2) - S0(i,2)
                    DS(i,3) = S1(i,3) - S0(i,3)
        
                END DO
        


                IF(nz .EQ. 1) THEN
                
                ! compute the "density decomposition" weights and current density contributions
                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = 0, 0
                                   
                            temp(i,j,k,1) = Jx(ix0+i,iy0+j,iz0+k)  
                            temp(i,j,k,2) = Jy(ix0+i,iy0+j,iz0+k)
                            temp(i,j,k,3) = Jz(ix0+i,iy0+j,iz0+k)  
                            
                        END DO
                    END DO    
                END DO 
                
                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = 0, 0
            
                            ! W1_i,j,k
                            W(i, j, k, 1) = DS(i,1) * (S0(j,2) +  0.5*DS(j,2))
                            ! W2_i,j,k
                            W(i, j, k, 2) = DS(j,2) * (S0(i,1) +  0.5*DS(i,1))
                            ! W3_i,j,k
                            W(i, j, k, 3) = S0(i,1) * S0(j,2) + 0.5*S0(j,2)*DS(i,1) + &
                                            0.5*S0(i,1)*DS(j,2) + third*DS(i,1)*DS(j,2)
                
                            ! add current contributions resulting from motion of this particle                
                            Jx(ix0+i,iy0+j,iz0+k) = Jx(ix0+i-1,iy0+j,iz0+k) - qdt * W(i, j, k, 1) 
                            Jy(ix0+i,iy0+j,iz0+k) = Jy(ix0+i,iy0+j-1,iz0+k) - qdt * W(i, j, k, 2)
                            Jz(ix0+i,iy0+j,iz0+k) = particles(sp,np)%vz * qdt * W(i, j, k, 3) 
                            
                        END DO
                    END DO    
                END DO 

                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = 0, 0
            
                            Jx(ix0+i,iy0+j,iz0+k) = Jx(ix0+i,iy0+j,iz0+k) + temp(i,j,k,1)  
                            Jy(ix0+i,iy0+j,iz0+k) = Jy(ix0+i,iy0+j,iz0+k) + temp(i,j,k,2)  
                            Jz(ix0+i,iy0+j,iz0+k) = Jz(ix0+i,iy0+j,iz0+k) + temp(i,j,k,3)  
                            
                        END DO
                    END DO    
                END DO

                ELSE
                
                
                ! compute the "density decomposition" weights and current density contributions
                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = -2, 2
                        
                            temp(i,j,k,1) = Jx(ix0+i,iy0+j,iz0+k)  
                            temp(i,j,k,2) = Jy(ix0+i,iy0+j,iz0+k)
                            temp(i,j,k,3) = Jz(ix0+i,iy0+j,iz0+k)  
                            
                        END DO
                    END DO    
                END DO 

                
                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = -2, 2
            
                            ! W1_i,j,k
                            W(i, j, k, 1) = DS(i,1) * (S0(j,2) +  0.5*DS(j,2)*S0(k,3) + &
                                            0.5*S0(j,2)*DS(k,3) + third*DS(j,2)*DS(k,3))
                            ! W2_i,j,k
                            W(i, j, k, 2) = DS(j,2) * (S0(i,1)*S0(k,3) +  0.5*DS(i,1)*S0(k,3) + &
                                            0.5*S0(i,1)*DS(k,3) + third*DS(i,1)*DS(k,3))
                            ! W3_i,j,k
                            W(i, j, k, 3) = DS(k,3) * (S0(i,1)*S0(j,2) +  0.5*DS(i,1)*S0(j,2) + &
                                            0.5*S0(i,1)*DS(j,2) + third*DS(i,1)*DS(j,2))
                
                
                            ! add current contributions resulting from motion of this particle                
                            Jx(ix0+i,iy0+j,iz0+k) = Jx(ix0+i-1,iy0+j,iz0+k) - qdt * W(i, j, k, 1) 
                            Jy(ix0+i,iy0+j,iz0+k) = Jy(ix0+i,iy0+j-1,iz0+k) - qdt * W(i, j, k, 2)
                            Jz(ix0+i,iy0+j,iz0+k) = Jz(ix0+i,iy0+j,iz0+k-1) - qdt * W(i, j, k, 3) 
                            
                        END DO
                    END DO    
                END DO 
   

                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = -2, 2
            
                            Jx(ix0+i,iy0+j,iz0+k) = Jx(ix0+i,iy0+j,iz0+k) + temp(i,j,k,1)  
                            Jy(ix0+i,iy0+j,iz0+k) = Jy(ix0+i,iy0+j,iz0+k) + temp(i,j,k,2)  
                            Jz(ix0+i,iy0+j,iz0+k) = Jz(ix0+i,iy0+j,iz0+k) + temp(i,j,k,3)  
                            
                        END DO
                    END DO    
                END DO 
                
                END IF
                         
        END DO
        
    END DO
    
    
    ! If periodic boundaries, need to touch up the grid deposited charges and currents at boundaries
    IF (bndry .EQ. 1) THEN
     
        IF(nranks_x .EQ. 1) THEN  
       
            buffer = Jx  
            Jx(0,:,:) = Jx(0,:,:) + Jx(nx,:,:)	
            Jx(1,:,:) = Jx(1,:,:) + Jx(nx+1,:,:)	
            Jx(nx,:,:) = Jx(nx,:,:) + buffer(0,:,:)
            Jx(nx-1,:,:) = Jx(nx-1,:,:) + Jx(-1,:,:)
        END IF
        
        IF(nranks_y .EQ. 1) THEN  
       
            buffer = Jy
            Jy(:,0,:) = Jy(:,0,:) + Jy(:,ny,:)	
            Jy(:,1,:) = Jy(:,1,:) + Jy(:,ny+1,:)	
            Jy(:,ny,:) = Jy(:,ny,:) + buffer(:,0,:)
            Jy(:,ny-1,:) = Jy(:,ny-1,:) + Jy(:,-1,:)
        
        END IF
        
        
        buffer = Jz

        Jz(:,:,0) = Jz(:,:,0) + Jz(:,:,nz)	
        Jz(:,:,1) = Jz(:,:,1) + Jz(:,:,nz+1)	

        Jz(:,:,nz) = Jz(:,:,nz) + buffer(:,:,0)
        Jz(:,:,nz-1) = Jz(:,:,nz-1) + Jz(:,:,-1)
        
    END IF
    
    
END SUBROUTINE deposit_currents_esirkepov



! conservative current deposition using Villasenor-Buneman algorithm 
SUBROUTINE deposit_currents_buneman()

    INTEGER :: sp, np
    REAL*8 :: x0, y0, z0, x1, y1, z1, igam, qq
    LOGICAL :: in
    
    
    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0
    
    ! Loop over all particles 
    DO sp = 1, ns   
    
        qq = q(sp) 
        
        DO np = 1, Np_in(sp)

            ! old position
            x0 = particles(sp,np)%x -  particles(sp,np)%vx
            y0 = particles(sp,np)%y -  particles(sp,np)%vy
            z0 = particles(sp,np)%z -  particles(sp,np)%vz

            ! new position
            x1 = particles(sp,np)%x
            y1 = particles(sp,np)%y
            z1 = particles(sp,np)%z
    
            ! deposit currents (done using 1D passes: x followed by y, then z)
            CALL xsplit(x1+2.d0,y1+2.d0,z1+2.d0,x0+2.d0,y0+2.d0,z0+2.d0,qq,in)

        END DO
    END DO


END SUBROUTINE deposit_currents_buneman


SUBROUTINE xsplit(x,y,z,x0,y0,z0,qq, in)

    REAL*8, INTENT(IN) :: x, y, z, x0, y0, z0, qq
    LOGICAL, INTENT(INOUT) :: in
    REAL*8 :: x1, y1, z1


    in = .TRUE.

    IF((INT(x) .NE. INT(x0)) .AND. ((x-x0) .NE. 0.d0)) THEN

        x1 = 0.5d0 * (1 + INT(x) + INT(x0))
        y1 = y0 + (y-y0) * ((x1-x0) / (x-x0))
        z1 = z0 + (z-z0) * ((x1-x0) / (x-x0))  
        CALL ysplit(x1,y1,z1,x0,y0,z0,qq,in) 
        
        IF(.NOT. in) RETURN
        
        in = (x1 .GT. 3.d0) .AND. (x1 .LT. mx-2.d0)
        IF(in) CALL ysplit(x,y,z,x1,y1,z1,qq,in)
        
 
    ELSE
        CALL ysplit(x,y,z,x0,y0,z0,qq,in)         
    END IF

END SUBROUTINE xsplit


SUBROUTINE ysplit(x,y,z,x0,y0,z0,qq,in)

    REAL*8, INTENT(IN) :: x, y, z, x0, y0, z0, qq
    LOGICAL, INTENT(INOUT) :: in
    REAL*8 :: x1, y1, z1

     IF((INT(y) .NE. INT(y0)) .AND. ((y-y0) .NE. 0.d0)) THEN

        y1 = 0.5d0 * (1 + INT(y) + INT(x0))
        z1 = z0 + (z-z0) * ((y1-y0) / (y-y0))  
        x1 = x0 + (x-x0) * ((y1-y0) / (y-y0))

        CALL zsplit(x1,y1,z1,x0,y0,z0,qq,in) 
        
        IF(.NOT. in) RETURN
        
        in = (y1 .GT. 3.d0) .AND. (y1 .LT. my-2.d0)
        IF(in) CALL zsplit(x,y,z,x1,y1,z1,qq,in)
        
 
    ELSE
        CALL zsplit(x,y,z,x0,y0,z0,qq,in)         
    END IF


END SUBROUTINE ysplit


SUBROUTINE zsplit(x,y,z,x0,y0,z0,qq,in)

    REAL*8, INTENT(IN) :: x, y, z, x0, y0, z0, qq
    LOGICAL, INTENT(INOUT) :: in
    REAL*8 :: x1, y1, z1

     IF((INT(z) .NE. INT(z0)) .AND. ((z-z0) .NE. 0.d0)) THEN

        z1 = 0.5d0 * (1 + INT(z) + INT(z0))
        x1 = x0 + (x-x0) * ((z1-z0) / (z-z0))   
        y1 = y0 + (y-y0) * ((z1-z0) / (z-z0))  
        
        CALL depsit(x1,y1,z1,x0,y0,z0,qq,Jx,Jy,Jz) 
        
        IF(.NOT. in) RETURN
        
        in = (z1 .GT. 3.d0) .AND. (z1 .LT. mz-2.d0)
        IF(in) CALL depsit(x,y,z,x1,y1,z1,qq,Jx,Jy,Jz)
        
    ELSE
        CALL depsit(x,y,z,x0,y0,z0,qq,Jx,Jy,Jz)         
    END IF


END SUBROUTINE zsplit


SUBROUTINE depsit(x,y,z,x0,y0,z0,qq,jx,jy,jz)

    REAL*8, INTENT(IN) :: x, y, z, x0, y0, z0, qq
    REAL*8, INTENT(INOUT) :: jx(1), jy(1), jz(1)  ! arrays treated as 1D
    INTEGER :: i, j, k, l, m, n
    REAL*8 :: delx, dely, delz, cx, cy, cz, &
              qu, qv, qw, delt, s, &
              su, sv, sw
    
    ! cell indices of half-way point
    i = 0.5d0 * (x+x0) 
    j = 0.5d0 * (y+y0) 
    k = 0.5d0 * (z+z0) 

    ! displacement in cell of half-way point
    delx = 0.5d0 * (x+x0) - i
    cx = 1.d0 - delx
    dely = 0.5d0 * (y+y0) - j
    cy = 1.d0 - dely
    delz = 0.5d0 * (z+z0) - k
    cz = 1.d0 - delz
    
    
    ! compute 1D equivalent cell index
    l = i + iy * (j-1) + iz * (k-1) 

    ! current elements
    qu = qq * (x-x0)
    qv = qq * (y-y0)
    qw = qq * (z-z0)
    delt = 0.08333333d0 * qu * (y-y0) * (z-z0)

    ! Add current contributions from this charge
    ! (OR can directly decrement smoothed currents from the electric field)
    !IF(current_filter_on) THEN
    
    GO TO 123
    DO n = 1, 27
        m  = ms(n) + l
        s  = sm(n) * delt
        su = sm(n) * qu
        sv = sm(n) * qv
        sw = sm(n) * qw
        
        !ex(m+iy+iz) = ex(m+iy+iz) - su*dely*delz-s
        !ex(m+iz) = ex(m+iz) - su*cy*delz+s
        !ex(m+iy) = ex(m+iy) - su*dely*cz+s
        !ex(m) = ex(m) - su*cy*cz-s
        
        !ey(m+iz+ix) = ey(m+iz+ix)-sv*delz*delx-s
        !ey(m+ix) = ey(m+ix)-sv*cz*delx+s
        !ey(m+iz) =  ey(m+iz)-sv*delz*cx+s
        !ey(m) = ey(m)-sv*cz*cx-s
      
        !ez(m+ix+iy) = ez(m+ix+iy)-sw*delx*dely-s
        !ez(m+iy) = ez(m+iy)-sw*cx*dely+s
        !ez(m+ix) = ez(m+ix)-sw*delx*cy+s
        !ez(m) = ez(m)-sw*cx*cy-s
        
        
        jx(m+iy+iz) = jx(m+iy+iz) + sm(n)*(qu*dely*delz+delt)
        jx(m+iz) = jx(m+iz) + sm(n)*(qu*(1.0-dely)*delz-delt)
        jx(m+iy) = jx(m+iy) + sm(n)*(qu*dely*(1.0-delz)-delt)
        jx(m) = jx(m) + sm(n)*(qu*(1.0-dely)*(1.0-delz)+delt)

        
        jy(m+iz+ix) = jy(m+iz+ix) + sm(n)*(qv*delz*delx+delt)
        jy(m+ix) = jy(m+ix) + sm(n)*(qv*(1.0-delz)*delx-delt)
        jy(m+iz) =  jy(m+iz) + sm(n)*(qv*delz*(1.0-delx)-delt)
        jy(m) = jy(m) + sm(n)*(qv*(1.0-delz)*(1.0-delx)+delt)
      
       
        jz(m+ix+iy) = jz(m+ix+iy) + sm(n)*(qw*delx*dely+delt)
        jz(m+iy) = jz(m+iy) + sm(n)*(qw*(1.0-delx)*dely-delt)
        jz(m+ix) = jz(m+ix)+ sm(n)*(qw*delx*(1.0-dely)-delt)
        jz(m) = jz(m) + sm(n)*(qw*(1.0-delx)*(1.0-dely)+delt)
    END DO
    123 CONTINUE
    
    !ELSE
        n = 14
        m  = ms(n) + l
        qu = 8.d0 * qu
        qv = 8.d0 * qv
        qw = 8.d0 * qw
        !ex(m+iy+iz) = ex(m+iy+iz) - su*dely*delz-s
        !ex(m+iz) = ex(m+iz) - su*cy*delz+s
        !ex(m+iy) = ex(m+iy) - su*dely*cz+s
        !ex(m) = ex(m) - su*cy*cz-s
        
        !ey(m+iz+ix) = ey(m+iz+ix)-sv*delz*delx-s
        !ey(m+ix) = ey(m+ix)-sv*cz*delx+s
        !ey(m+iz) =  ey(m+iz)-sv*delz*cx+s
        !ey(m) = ey(m)-sv*cz*cx-s
      
        !ez(m+ix+iy) = ez(m+ix+iy)-sw*delx*dely-s
        !ez(m+iy) = ez(m+iy)-sw*cx*dely+s
        !ez(m+ix) = ez(m+ix)-sw*delx*cy+s
        !ez(m) = ez(m)-sw*cx*cy-s
        
        
        jx(m+iy+iz) = jx(m+iy+iz) + sm(n)*(qu*dely*delz+delt)
        jx(m+iz) = jx(m+iz) + sm(n)*(qu*(1.0-dely)*delz-delt)
        jx(m+iy) = jx(m+iy) + sm(n)*(qu*dely*(1.0-delz)-delt)
        jx(m) = jx(m) + sm(n)*(qu*(1.0-dely)*(1.0-delz)+delt)

        
        jy(m+iz+ix) = jy(m+iz+ix) + sm(n)*(qv*delz*delx+delt)
        jy(m+ix) = jy(m+ix) + sm(n)*(qv*(1.0-delz)*delx-delt)
        jy(m+iz) =  jy(m+iz) + sm(n)*(qv*delz*(1.0-delx)-delt)
        jy(m) = jy(m) + sm(n)*(qv*(1.0-delz)*(1.0-delx)+delt)
      
       
        jz(m+ix+iy) = jz(m+ix+iy) + sm(n)*(qw*delx*dely+delt)
        jz(m+iy) = jz(m+iy) + sm(n)*(qw*(1.0-delx)*dely-delt)
        jz(m+ix) = jz(m+ix)+ sm(n)*(qw*delx*(1.0-dely)-delt)
        jz(m) = jz(m) + sm(n)*(qw*(1.0-delx)*(1.0-dely)+delt)
   !END IF    

END SUBROUTINE depsit




! deposits charge to grid
SUBROUTINE deposit_charges()

	INTEGER :: i, j, k, ix, iy, iz, sp, np
	REAL*8 :: x, y, z, qdxyz, Sint(-1:2,3)

    ! clear density
    rho = 0.d0
    !rho(1:nx,1:ny,1:nz) = rho0   -> uniform background 
	ne = 0.d0
	ni = 0.d0
	

    ! loop over all particles
    DO sp = 1,ns   
    
        qdxyz = q(sp)
        
        DO np = 1, Np_in(sp)

                x = particles(sp,np)%x/dx
                y = particles(sp,np)%y/dy
                z = particles(sp,np)%z/dz

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
				
                ELSE IF(interpolation_type .EQ. 2) THEN

                    ix = x 
                    iy = y
                    iz = z       
        
                    ! compute the 1st order interpolation co-efficients
                    DO i = -1,2
                        Sint(i,1) = S1_1D((ix+i)-x)
                        Sint(i,2) = S1_1D((iy+i)-y)
                        Sint(i,3) = S1_1D((iz+i)-z)
                    END DO       
       
                    ! now interpolate from grid to particle location
                    DO k = -1,2
                        DO j = -1,2
                            DO i = -1,2
                                rho(ix+i,iy+j,iz+k) = rho(ix+i,iy+j,iz+k) + qdxyz * Sint(i,1)*Sint(j,2)*Sint(k,3)
                                IF(i .EQ. 1) ne(ix+i,iy+j,iz+k) = ne(ix+i,iy+j,iz+k) +  Sint(i,1)*Sint(j,2)*Sint(k,3)
                                IF(i .EQ. 2) ni(ix+i,iy+j,iz+k) = ni(ix+i,iy+j,iz+k) +  Sint(i,1)*Sint(j,2)*Sint(k,3) 
                            END DO
                        END DO
                    END DO        

                ELSE IF( interpolation_type .EQ. 3) THEN

                    ix = x 
                    iy = y
                    iz = z       
        
                    ! compute the 2nd order interpolation co-efficients
                    DO i = -1,2
                        Sint(i,1) = S2_1D((ix+i)-x)
                        Sint(i,2) = S2_1D((iy+i)-y)
                        Sint(i,3) = S2_1D((iz+i)-z)
                    END DO       
       
                    ! now interpolate from grid to particle location
                    DO k = -1,2
                        DO j = -1,2
                            DO i = -1,2
                                rho(ix+i,iy+j,iz+k) = rho(ix+i,iy+j,iz+k) + qdxyz * Sint(i,1)*Sint(j,2)*Sint(k,3)
                                IF(i .EQ. 1) ne(ix+i,iy+j,iz+k) = ne(ix+i,iy+j,iz+k) +  Sint(i,1)*Sint(j,2)*Sint(k,3)
                                IF(i .EQ. 2) ni(ix+i,iy+j,iz+k) = ni(ix+i,iy+j,iz+k) +  Sint(i,1)*Sint(j,2)*Sint(k,3) 
                            END DO
                        END DO
                    END DO        

     
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
! Implemented as three successive 1D sweeps (x,y and z directions)
SUBROUTINE filter(f)

    REAL*8, INTENT(INOUT) ::  f(-1:nx+3,-1:ny+3,-1:nz+3)
    INTEGER :: i, j, k, ix1, ix2, iy1, iy2, iz1, iz2

    
    ! Equivalent to successively apply 1D binomial filters separately in x, y, and z directions 
    ! e.g. 1D_x_Filter[ f_i,j,k ] = 0.25*f_i-1,j,k + 0.5*f_i,j,k + 0.25*f_i+1,j,k 
    !      1D_y_Filter[ f_i,j,k ] = 0.25*f_i,j-1,k + 0.5*f_i,j,k + 0.25*f_i,j+1,k 
    !      1D_z_Filter[ f_i,j,k ] = 0.25*f_i,j,k-1 + 0.5*f_i,j,k + 0.25*f_i,j,k+1 
    


    IF(mycoord(1)  .EQ. 0) THEN
        ix1 = 1
    ELSE
       ix1 = 0
    END IF

    IF(mycoord(1)  .EQ. nranks_x-1) THEN
        ix2 = 1
    ELSE
       ix2 = 0
    END IF
    
    IF(mycoord(2)  .EQ. 0) THEN
        iy1 = 1
    ELSE
       iy1 = 0
    END IF

    IF(mycoord(2)  .EQ. nranks_y-1) THEN
        iy2 = 1
    ELSE
       iy2 = 0
    END IF

       
        buffer = f

        !$OMP PARALLEL DO
        DO k = 0, nz+1
            DO j = 1-iy1, ny+iy2
                DO i = 1-ix1, nx+ix2
                    f(i,j,k) = 0.015625* ( &
                               buffer(i-1,j-1,k-1) + 2.d0*buffer(i,j-1,k-1) + buffer(i+1,j-1,k-1)    + &
                               2.d0*(buffer(i-1,j,k-1) + 2.d0*buffer(i,j,k-1) + buffer(i+1,j,k-1))   + &
                               buffer(i-1,j+1,k-1) + 2.d0*buffer(i,j+1,k-1) + buffer(i+1,j+1,k-1)    + &
                               
                               2.d0 * ( buffer(i-1,j-1,k) + 2.d0*buffer(i,j-1,k) + buffer(i+1,j-1,k) + &
                               2.d0*( buffer(i-1,j,k) + 2.d0*buffer(i,j,k) + buffer(i+1,j,k))        + &
                               buffer(i-1,j+1,k) + 2.d0*buffer(i,j+1,k) + buffer(i+1,j+1,k) )        + &
                               
                               buffer(i-1,j-1,k+1) + 2.d0*buffer(i,j-1,k+1) + buffer(i+1,j-1,k+1)    + &
                               2.d0*(buffer(i-1,j,k+1) + 2.d0*buffer(i,j,k+1) + buffer(i+1,j,k+1))   + &
                               buffer(i-1,j+1,k+1) + 2.d0*buffer(i,j+1,k+1) + buffer(i+1,j+1,k+1)  )                     
                END DO
            END DO        
        END DO    
        !$OMP END PARALLEL DO



    
    
END SUBROUTINE filter


! enforces particle boundary conditions
SUBROUTINE particle_boundary()

    INTEGER :: i, counter
    
    ! reset escape particles counter
    Np_esc_xm = 0
    Np_esc_xp = 0
    Np_esc_ym = 0
    Np_esc_yp = 0
    Np_esc_bl = 0
    Np_esc_br = 0
    Np_esc_tl = 0
    Np_esc_tr = 0
    
    ! clear boundary particle buffers
    bparticles_xm = 0.d0
    bparticles_xp = 0.d0
    bparticles_ym = 0.d0
    bparticles_yp = 0.d0
    bparticles_bl = 0.d0
    bparticles_br = 0.d0
    bparticles_tl = 0.d0
    bparticles_tr = 0.d0
    
    
    Ntot_esc(1) = Np_in(1)                         
    Ntot_esc(2) = Np_in(2)
    
    !PRINT*,''
    !PRINT*,'Np_in = ',Np_in
    !PRINT*,''
    
    
    
    ! loop over all particles
    DO i= 1,ns 
    
        counter = 1        
        DO WHILE(counter .LE. Np_in(i))

            ! periodic boundary conditions
            IF(bndry .eq. 1) THEN
      
                
                IF(particles(i,counter)%x .LT. xmin) THEN
                    IF(nranks_x .EQ. 1) THEN
                        particles(i,counter)%x = particles(i,counter)%x + Lx
                    ELSE
                        particles(i,counter)%oob = .TRUE.                      
                        ! delete out of bounds particle from array
                        CALL delete_particle(i,counter)                  
                        counter = counter - 1
                    END IF
               END IF
                
                IF(particles(i,counter)%x .GT. xmax) THEN
                    IF(nranks_x .EQ. 1) THEN
                        particles(i,counter)%x = particles(i,counter)%x - Lx                    
                    ELSE
                        particles(i,counter)%oob = .TRUE.                      
                        ! delete out of bounds particle from array
                        CALL delete_particle(i,counter)                  
                        counter = counter - 1
                    END IF
                END IF
                
                
                IF(particles(i,counter)%y .LT. ymin) THEN
                    IF(nranks_y .EQ. 1) THEN
                        particles(i,counter)%y = particles(i,counter)%y + Ly
                    ELSE
                        particles(i,counter)%oob = .TRUE.                      
                        ! delete out of bounds particle from array
                        CALL delete_particle(i,counter)                  
                        counter = counter - 1                   
                    END IF
                END IF
                
                IF(particles(i,counter)%y .GT. ymax) THEN
                    IF(nranks_y .EQ. 1) THEN
                        particles(i,counter)%y = particles(i,counter)%y - Ly
                    ELSE
                        particles(i,counter)%oob = .TRUE.                      
                        ! delete out of bounds particle from array
                        CALL delete_particle(i,counter)                  
                        counter = counter - 1
                    END IF        
                END IF
                
                
                
                IF(particles(i,counter)%z .LT. zmin) THEN
                    particles(i,counter)%z = particles(i,counter)%z + Lz
                END IF
                
                IF(particles(i,counter)%z .GT. zmax) THEN
                    particles(i,counter)%z = particles(i,counter)%z - Lz
                END IF
                 
                counter = counter + 1
         

            ! outflow boundary conditions (in this case, need to distribute ions and electrons uniformly across boundary
            ! to avoid net charge buildup)
            ELSE IF(bndry .EQ. 2) THEN 

            
                ! set out of bounds flag
                IF( particles(i,counter)%x .LT. xmin .OR. particles(i,counter)%x .GT. xmax .OR. &
                    particles(i,counter)%y .LT. ymin .OR. particles(i,counter)%y .GT. ymax .OR. &
                    particles(i,counter)%z .LT. zmin .OR. particles(i,counter)%z .GT. zmax) THEN

     
                    !PRINT*,'######################################################################'
                    !PRINT*,'#######################################################################'                     
                    !PRINT*,'OOB PARTICLE!!! Species, counter, x, y, z = ',i,counter,particles(i,counter)%x,particles(i,counter)%y
                    !PRINT*,'#######################################################################'
                    !PRINT*,'######################################################################'


                    particles(i,counter)%oob = .TRUE.                      
                
                    ! delete out of bounds particle from array
                    CALL delete_particle(i,counter)                    
             

                ELSE
                    counter = counter + 1 
                END IF 
            
            ELSE    
                EXIT    
            END IF
    
        END DO
        
    END DO
    
    
    !PRINT*,''
    
    Ntot_esc(1) = Ntot_esc(1) - Np_in(1)
    Ntot_esc(2) = Ntot_esc(2) - Np_in(2)
    
    !PRINT*,'Myrank, N_esc(electrons), N_esc(ions) = ',Myrank, Ntot_esc(1), Ntot_esc(2)
    
    !PRINT*,''

END SUBROUTINE particle_boundary



! particle deletion involvles re-sorting: the deleted particle is exchanged with the last particle on the list
SUBROUTINE delete_particle(i,j)

    INTEGER, INTENT(IN) :: i, j
    REAL*8 :: temp(6)
    INTEGER :: ix
   
    ! check if theres space in the array
    IF(Np_esc_xm .GE. max_particles .OR. Np_esc_xp .GE. max_particles .OR. &
       Np_esc_ym .GE. max_particles .OR. Np_esc_yp .GE. max_particles .OR. &
       Np_esc_bl .GE. max_particles .OR. Np_esc_br .GE. max_particles .OR. &
       Np_esc_tl .GE. max_particles .OR. Np_esc_tr .GE. max_particles ) THEN
       
        PRINT*,'ERROR! Particle boundary array full... Terminating program.'
        STOP        
    END IF
        
        
    ! replace this particle with the last particle in our array, and decrement Np_in    
    temp(1) = particles(i,j)%x    
    temp(2) = particles(i,j)%y
    temp(3) = particles(i,j)%z
    temp(4) = particles(i,j)%vx
    temp(5) = particles(i,j)%vy
    temp(6) = particles(i,j)%vz
    
    particles(i,j)%x = particles(i,Np_in(i))%x    
    particles(i,j)%y = particles(i,Np_in(i))%y
    particles(i,j)%z = particles(i,Np_in(i))%z
    particles(i,j)%vx = particles(i,Np_in(i))%vx
    particles(i,j)%vy = particles(i,Np_in(i))%vy
    particles(i,j)%vz = particles(i,Np_in(i))%vz
    
    
    IF(temp(3) .GT. zmin .AND. temp(3) .LT. zmax) THEN
    
    
    ! x- boundary    
    IF(temp(1) .LT. xmin .AND. temp(2) .LT. ymax .AND. temp(2) .GT. ymin)THEN
 
        IF(neighbor_rank(1) .NE. -2) THEN
 
        Np_esc_xm = Np_esc_xm + 1 
         
        ix =  1+ (Np_esc_xm-1) * nvars_particles
         
        bparticles_xm(ix)   = i    
        bparticles_xm(ix+1) = temp(1) + Lx    
        bparticles_xm(ix+2) = temp(2) 
        bparticles_xm(ix+3) = temp(3)
        bparticles_xm(ix+4) = temp(4)
        bparticles_xm(ix+5) = temp(5)
        bparticles_xm(ix+6) = temp(6)
   
        END IF
  
    ! x+ boundary
    ELSE IF(temp(1) .GT. xmax .AND. temp(2) .LT. ymax .AND. temp(2) .GT. ymin)THEN

        IF(neighbor_rank(2) .NE. -2) THEN

        Np_esc_xp = Np_esc_xp + 1 

        ix =  1 + (Np_esc_xp-1) * nvars_particles
         
        bparticles_xp(ix)   = i    
        bparticles_xp(ix+1) = temp(1) - Lx   
        bparticles_xp(ix+2) = temp(2) 
        bparticles_xp(ix+3) = temp(3)
        bparticles_xp(ix+4) = temp(4)
        bparticles_xp(ix+5) = temp(5)
        bparticles_xp(ix+6) = temp(6)

        END IF 
         
    ! y- boundary        
    ELSE IF(temp(2) .LT. ymin .AND. temp(1) .LT. xmax .AND. temp(1) .GT. xmin)THEN

        IF(neighbor_rank(3) .NE. -2) THEN

        Np_esc_ym = Np_esc_ym + 1 
         
        ix =  1 + (Np_esc_ym-1) * nvars_particles
         
        bparticles_ym(ix)   = i    
        bparticles_ym(ix+1) = temp(1)    
        bparticles_ym(ix+2) = temp(2) + Ly
        bparticles_ym(ix+3) = temp(3) 
        bparticles_ym(ix+4) = temp(4)
        bparticles_ym(ix+5) = temp(5)
        bparticles_ym(ix+6) = temp(6)
              
        END IF
        
    ! y+ boundary
    ELSE IF(temp(2) .GT. ymax  .AND. temp(1) .LT. xmax .AND. temp(1) .GT. xmin)THEN
        
        IF(neighbor_rank(4) .NE. -2) THEN

        Np_esc_yp = Np_esc_yp + 1 
         
        ix =  1+ (Np_esc_yp-1) * nvars_particles
         
        bparticles_yp(ix)   = i    
        bparticles_yp(ix+1) = temp(1)    
        bparticles_yp(ix+2) = temp(2) - Ly
        bparticles_yp(ix+3) = temp(3) 
        bparticles_yp(ix+4) = temp(4)
        bparticles_yp(ix+5) = temp(5)
        bparticles_yp(ix+6) = temp(6)

        END IF

    ! x-y- edge 
    ELSE IF(temp(2) .LT. ymin  .AND. temp(1) .LT. xmin)THEN
        
        IF(neighbor_rank(5) .NE. -2) THEN

        Np_esc_bl = Np_esc_bl + 1 
         
        ix =  1+ (Np_esc_bl-1) * nvars_particles
         
        bparticles_bl(ix)   = i    
        bparticles_bl(ix+1) = temp(1) + Lx   
        bparticles_bl(ix+2) = temp(2) + Ly
        bparticles_bl(ix+3) = temp(3) 
        bparticles_bl(ix+4) = temp(4)
        bparticles_bl(ix+5) = temp(5)
        bparticles_bl(ix+6) = temp(6)
        
        END IF    

    ! x+y- edge 
    ELSE IF(temp(2) .LT. ymin  .AND. temp(1) .GT. xmax)THEN
        
        IF(neighbor_rank(6) .NE. -2) THEN

        Np_esc_br = Np_esc_br + 1 
         
        ix =  1+ (Np_esc_br-1) * nvars_particles
         
        bparticles_br(ix)   = i    
        bparticles_br(ix+1) = temp(1) - Lx   
        bparticles_br(ix+2) = temp(2) + Ly
        bparticles_br(ix+3) = temp(3) 
        bparticles_br(ix+4) = temp(4)
        bparticles_br(ix+5) = temp(5)
        bparticles_br(ix+6) = temp(6)
 
        END IF

    ! x-y+ edge 
    ELSE IF(temp(2) .GT. ymax  .AND. temp(1) .LT. xmin)THEN
        
        IF(neighbor_rank(7) .NE. -2) THEN

        Np_esc_tl = Np_esc_tl + 1 
         
        ix =  1+ (Np_esc_tl-1) * nvars_particles
         
        bparticles_tl(ix)   = i    
        bparticles_tl(ix+1) = temp(1) + Lx   
        bparticles_tl(ix+2) = temp(2) - Ly
        bparticles_tl(ix+3) = temp(3) 
        bparticles_tl(ix+4) = temp(4)
        bparticles_tl(ix+5) = temp(5)
        bparticles_tl(ix+6) = temp(6)
        
        END IF
            
    ! x+y+ edge 
    ELSE IF(temp(2) .GT. ymax  .AND. temp(1) .GT. xmax)THEN
        
        IF(neighbor_rank(8) .NE. -2) THEN

        Np_esc_tr = Np_esc_tr + 1 
         
        ix =  1+ (Np_esc_tr-1) * nvars_particles
         
        bparticles_tr(ix)   = i    
        bparticles_tr(ix+1) = temp(1) - Lx   
        bparticles_tr(ix+2) = temp(2) - Ly
        bparticles_tr(ix+3) = temp(3) 
        bparticles_tr(ix+4) = temp(4)
        bparticles_tr(ix+5) = temp(5)
        bparticles_tr(ix+6) = temp(6)        

        END IF

    END IF        
    
    
    ELSE
    ! Dont need to care about z boundary
    ! since domain decomposition is 2d in x,y   
    END IF    
       
    Np_in(i) = Np_in(i) - 1

END SUBROUTINE delete_particle



SUBROUTINE create_particle(sp,x,y,z,vx,vy,vz)

    INTEGER, INTENT(IN) :: sp
    REAL*8, INTENT(IN) :: x,y,z,vx,vy,vz

    ! insert new particle at the end of our array, and increment Np_in    
    Np_in(sp) = Np_in(sp) + 1
   
    particles(sp,Np_in(sp))%x = x   
    particles(sp,Np_in(sp))%y = y
    IF(nz .GT. 1) THEN
        particles(sp,Np_in(sp))%z = z
    ELSE
       particles(sp,Np_in(sp))%z = nz
    END IF

    particles(sp,Np_in(sp))%vx = vx   
    particles(sp,Np_in(sp))%vy = vy
    particles(sp,Np_in(sp))%vz = vz
    
    IF(vx**2 + vy**2 + vz**2 .GE. c**2) THEN
        PRINT*,'ERROR! Attempted to create particle with v>c..'
        PRINT*,'vx, vy, vz, |v| = ', vx, vy, vz, SQRT(vx**2 + vy**2 + vz**2)
        STOP
    END IF    
    
    particles(sp,Np_in(sp))%oob = .FALSE.
    particles(sp,Np_in(sp))%species = sp
   
    !IF( x .LT. xmin .OR. x .GT. xmax .OR. y .LT. ymin .OR. y .GT. ymax .OR. &
    !    z .LT. zmin .OR. z .GT. zmax) THEN
      
    !    PRINT*,'ERROR...out of bound particle creation.'
    !    PRINT*,'x, y, z =',x,y,z
    !    STOP      
                    
    !END IF                 
    
END SUBROUTINE create_particle




END MODULE particleMover_mod