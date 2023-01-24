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
    REAL*8 :: t(3), s(3), gam, tsqr1, w_d

		
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
        !$OMP PRIVATE(x, y, z, a, b, gam, t, tsqr1, s, uxx, uyy, uzz, uxp, uyp, uzp)
        DO j = 1, Np_in(i)

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
        
                ! half electric push (u^n-1/2 ->u^-) 
                particles(i,j)%ux = particles(i,j)%ux + 0.5 * a(1)
                particles(i,j)%uy = particles(i,j)%uy + 0.5 * a(2)
                particles(i,j)%uz = particles(i,j)%uz + 0.5 * a(3)
			
                !IF(i .EQ. 1) particles(i,j)%ux = 0.3

			
                gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2 + particles(i,j)%uz**2 )/c**2 )
                t     = 0.5d0 * b / (gam*c)
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
                particles(i,j)%ux = particles(i,j)%ux + 0.5 * a(1)
                particles(i,j)%uy = particles(i,j)%uy + 0.5 * a(2)
                particles(i,j)%uz = particles(i,j)%uz + 0.5 * a(3)
                             
                K_E = K_E + (gam-1.d0)* m(i) * (c**2)
                !Px = Px + m(i)*particles(i,j)%ux
                !Py = Py + m(i)*particles(i,j)%uy
                !Pz = Pz + m(i)*particles(i,j)%uz
			
            
                !***********************
                ! update x(t) -> x(t+dt)
                !***********************
                !gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2 + particles(i,j)%uz**2 )/c**2 )
                particles(i,j)%x = particles(i,j)%x + particles(i,j)%ux / gam
                particles(i,j)%y = particles(i,j)%y + particles(i,j)%uy / gam       
                particles(i,j)%z = particles(i,j)%z + particles(i,j)%uz / gam       
 
                
                !IF(i .EQ. 1) THEN
                !    particles(i,j)%ux =  (5.d0*twopi/100.d0)*COS(twopi*ts/100.d0)
                !    particles(i,j)%x =    nx/2 + 5.d0*SIN(twopi*ts/100.d0)
                !ELSE IF(i .EQ. 2) THEN
                    !particles(i,j)%ux = -(5.d0*twopi/200.d0)*COS(twopi*ts/100.d0)
                    !particles(i,j)%x = nx/2 - 5.d0*SIN(twopi*ts/100.d0)
                !END IF
            
		END DO
        !$OMP END PARALLEL DO

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
            igam = 1.d0 /  &
            SQRT( 1.d0 + (particles(sp,np)%ux**2 + particles(sp,np)%uy**2 + particles(sp,np)%uz**2 )/c**2 )

            x0 = particles(sp,np)%x - igam * particles(sp,np)%ux
            y0 = particles(sp,np)%y - igam * particles(sp,np)%uy
            z0 = particles(sp,np)%z - igam * particles(sp,np)%uz

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
    IF(current_filter_on) THEN
    
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
    
    ELSE
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
    END IF    

END SUBROUTINE depsit




! conservative current deposition using Umeda (2003) Zig-zag algorithm
SUBROUTINE deposit_currents_zigzag()

    REAL*8 :: x1, y1, z1, x2, y2, z2
    REAL*8 :: qdt
    REAL*8 :: xr, yr, zr    
    REAL*8 :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2
    REAL*8 :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2    
    INTEGER :: sp, np
    INTEGER :: i1, i2, j1, j2, k1, k2
    REAL*8  :: idx, idy, idz, idxyz, igam         


    idx = 1.d0 / dx 
    idy = 1.d0 / dy 
    idz = 1.d0 / dz 
    idxyz = 1.d0 / (dx*dy*dz) 

    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0

    !GO TO 1111
    
    ! loop over all particles
    DO sp = 1, ns   
    
        qdt = q(sp) 
        
        !PRINT*,'sp, qdt = ', sp, qdt
        
        DO np = 1, Np_in(sp)

            !PRINT*,'np, Np_in', np, Np_in

            ! only contributins from particles that are inside the computational domain

                igam = 1.d0 /  &
                SQRT( 1.d0 + (particles(sp,np)%ux**2 + particles(sp,np)%uy**2 + particles(sp,np)%uz**2 )/c**2 )

                ! old position
                x1 = particles(sp,np)%x - igam * particles(sp,np)%ux
                y1 = particles(sp,np)%y - igam * particles(sp,np)%uy
                z1 = particles(sp,np)%z - igam * particles(sp,np)%uz

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
                xr = MIN( MIN(DBLE(i1*dx), DBLE(i2*dx)) + dx, MAX( MAX(DBLE(i1*dx), DBLE(i2*dx)), 0.5d0*(x1+x2) ) )
                yr = MIN( MIN(DBLE(j1*dy), DBLE(j2*dy)) + dy, MAX( MAX(DBLE(j1*dy), DBLE(j2*dy)), 0.5d0*(y1+y2) ) )
                zr = MIN( MIN(DBLE(k1*dz), DBLE(k2*dz)) + dz, MAX( MAX(DBLE(k1*dz), DBLE(k2*dz)), 0.5d0*(z1+z2) ) )
    
    
                ! compute segment flux and shape-factors                
                Fx1 = qdt * (xr - x1)
                Fy1 = qdt * (yr - y1)
                Fz1 = qdt * (zr - z1)
    
                Fx2 = qdt * (x2 - xr)
                Fy2 = qdt * (y2 - yr)
                Fz2 = qdt * (z2 - zr)
    
                Wx1 = 0.5d0 * idx * (x1 + xr) - i1
                Wy1 = 0.5d0 * idy * (y1 + yr) - j1
                Wz1 = 0.5d0 * idz * (z1 + zr) - k1

                Wx2 = 0.5d0 * idx * (xr + x2) - i2
                Wy2 = 0.5d0 * idy * (yr + y2) - j2
                Wz2 = 0.5d0 * idz * (zr + z2) - k2


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
               
            
        END DO        
    END DO
    
     
    !1111 CONTINUE
   
    ! If periodic boundaries, need to touch up the grid deposited charges and currents at boundaries
    !IF (bndry .EQ. 1) THEN

        !buffer = Jx
        !Jx(0,:,:) = Jx(0,:,:) + Jx(nx,:,:)	
        !Jx(1,:,:) = Jx(1,:,:) + Jx(nx+1,:,:)	
        !Jx(nx,:,:) = Jx(nx,:,:) + buffer(0,:,:)
        !Jx(nx-1,:,:) = Jx(nx-1,:,:) + Jx(-1,:,:)

        !buffer = Jy
        !Jy(:,0,:) = Jy(:,0,:) + Jy(:,ny,:)	
        !Jy(:,1,:) = Jy(:,1,:) + Jy(:,ny+1,:)	
        !Jy(:,ny,:) = Jy(:,ny,:) + buffer(:,0,:)
        !Jy(:,ny-1,:) = Jy(:,ny-1,:) + Jy(:,-1,:)
        
        !buffer = Jz
        !Jz(:,:,0) = Jz(:,:,0) + Jz(:,:,nz)	
        !Jz(:,:,1) = Jz(:,:,1) + Jz(:,:,nz+1)	
        !Jz(:,:,nz) = Jz(:,:,nz) + buffer(:,:,0)
        !Jz(:,:,nz-1) = Jz(:,:,nz-1) + Jz(:,:,-1)
        
    !END IF
    
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
    

END SUBROUTINE deposit_currents_zigzag



! conservative current deposition using Esirkepov (2001) algorithm
SUBROUTINE deposit_currents_esirkepov()

    REAL*8 :: x0, y0, z0, x1, y1, z1
    INTEGER :: ix0, iy0, iz0, i, j, k
    REAL*8  :: qdt, idxy, idyz, idxz, igam
    REAL*8  :: S0(-2:2,3), S1(-2:2,3), DS(-2:2,3)
    REAL*8  :: W(-2:2, -2:2, -2:2, 3)
    INTEGER :: sp, np
         
    idxy = 1.d0 / (dx * dy)
    idxz = 1.d0 / (dx * dz)
    idyz = 1.d0 / (dy * dz)

    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0


    ! loop over all particles
    DO sp = 1, ns   
    
        qdt = q(sp) 
        
        DO np = 1, Np_in(sp)

            ! only contributins from particles that are inside the computational domain
            IF(.NOT. particles(sp,np)%oob) THEN

                ! old position
                igam = 1.d0 /  &
                SQRT( 1.d0 + (particles(sp,np)%ux**2 + particles(sp,np)%uy**2 + particles(sp,np)%uz**2 )/c**2 )


                x0 = particles(sp,np)%x - igam * particles(sp,np)%ux
                y0 = particles(sp,np)%y - igam * particles(sp,np)%uy
                z0 = particles(sp,np)%z - igam * particles(sp,np)%uz


                ! new position
                x1 = particles(sp,np)%x
                y1 = particles(sp,np)%y
                z1 = particles(sp,np)%z
    
                !compute index of occupied cell at initial position
                ix0 = x0 + 0.5  
                iy0 = y0 + 0.5
                iz0 = z0 + 0.5
        

                ! prepare array of 1D form factors corresponding to initial and final positions and their difference 
                DO i = -2, 2
    
                    S0(i,1) = S_1D(DBLE(ix0+i)-x0)
                    S0(i,2) = S_1D(DBLE(iy0+i)-y0)
                    S0(i,3) = S_1D(DBLE(iz0+i)-z0)
        
                    S1(i,1) = S_1D(DBLE(ix0+i)-x1)
                    S1(i,2) = S_1D(DBLE(iy0+i)-y1)
                    S1(i,3) = S_1D(DBLE(iz0+i)-z1)
        
                    DS(i,1) = S1(i,1) - S0(i,1)
                    DS(i,2) = S1(i,2) - S0(i,2)
                    DS(i,3) = S1(i,3) - S0(i,3)
        
                END DO
        
                ! compute the "density decomposition" weights and current density contributions
                DO i = -2, 2 
                    DO j = -2, 2 
                        DO k = -2, 2
            
                            ! W1_i,j,k
                            W(i, j, k, 1) = DS(i,1) * (S0(j,2)*S0(k,3) +  0.5*DS(j,2)*S0(k,3) + &
                                            0.5*S0(j,2)*DS(k,3) + third*DS(j,2)*DS(k,3))
                            ! W2_i,j,k
                            W(i, j, k, 2) = DS(j,2) * (S0(i,1)*S0(k,3) +  0.5*DS(i,1)*S0(k,3) + &
                                            0.5*S0(i,1)*DS(k,3) + third*DS(i,1)*DS(k,3))
                            ! W3_i,j,k
                            W(i, j, k, 3) = DS(k,3) * (S0(i,1)*S0(j,2) +  0.5*DS(i,1)*S0(j,2) + &
                                            0.5*S0(i,1)*DS(j,2) + third*DS(i,1)*DS(j,2))
                
                            ! add current contributions resulting from motion of this particle                
                            Jx(ix0+i,iy0+j,iz0+k) = Jx(ix0+i,iy0+j,iz0+k) + Jx(ix0+i-1,iy0+j,iz0+k) &
                                                    - qdt * W(i, j, k, 1) 
                            Jy(ix0+i,iy0+j,iz0+k) = Jy(ix0+i,iy0+j,iz0+k) + Jy(ix0+i,iy0+j-1,iz0+k) &
                                                    - qdt * W(i, j, k, 2)
                            Jz(ix0+i,iy0+j,iz0+k) = Jz(ix0+i,iy0+j,iz0+k) + Jz(ix0+i,iy0+j,iz0+k-1) &
                                                    - qdt * W(i, j, k, 3) 
                            
                        END DO
                    END DO    
                END DO 

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
    

END SUBROUTINE deposit_currents_esirkepov



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
    
        qdxyz = q(i)/(dx*dy*dz)
        
        DO j = 1, Np_in(i)

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
                        ne(ix,iy,iz)       = ne(ix,iy,iz)       + dxx*dyy*dzz				     
                        ne(ix+1,iy,iz)     = ne(ix+1,iy,iz)     + delx*dyy*dzz
                        ne(ix,iy+1,iz)     = ne(ix,iy+1,iz)     + dxx*dely*dzz
                        ne(ix,iy,iz+1)     = ne(ix,iy,iz+1)     + dxx*dyy*delz
                        ne(ix+1,iy+1,iz)   = ne(ix+1,iy+1,iz)   + delx*dely*dzz
                        ne(ix+1,iy,iz+1)   = ne(ix+1,iy,iz+1)   + delx*dyy*delz
                        ne(ix,iy+1,iz+1)   = ne(ix,iy+1,iz+1)   + dxx*dely*delz
                        ne(ix+1,iy+1,iz+1) = ne(ix+1,iy+1,iz+1) + qdxyz*delx*dely*delz				
                    END IF
                
                    IF(i .EQ. 2)THEN
                        ni(ix,iy,iz)       = ni(ix,iy,iz)       + dxx*dyy*dzz				     
                        ni(ix+1,iy,iz)     = ni(ix+1,iy,iz)     + delx*dyy*dzz
                        ni(ix,iy+1,iz)     = ni(ix,iy+1,iz)     + dxx*dely*dzz
                        ni(ix,iy,iz+1)     = ni(ix,iy,iz+1)     + dxx*dyy*delz
                        ni(ix+1,iy+1,iz)   = ni(ix+1,iy+1,iz)   + delx*dely*dzz
                        ni(ix+1,iy,iz+1)   = ni(ix+1,iy,iz+1)   + delx*dyy*delz
                        ni(ix,iy+1,iz+1)   = ni(ix,iy+1,iz+1)   + dxx*dely*delz
                        ni(ix+1,iy+1,iz+1) = ni(ix+1,iy+1,iz+1) + delx*dely*delz		
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
! Implemented as three successive 1D sweeps (x,y and z directions)
SUBROUTINE filter(f)

    REAL*8, INTENT(INOUT) ::  f(-1:nx+3,-1:ny+3,-1:nz+3)
    INTEGER :: i, j, k, l

    
    ! Equivalent to successively apply 1D binomial filters separately in x, y, and z directions 
    ! e.g. 1D_x_Filter[ f_i,j,k ] = 0.25*f_i-1,j,k + 0.5*f_i,j,k + 0.25*f_i+1,j,k 
    !      1D_y_Filter[ f_i,j,k ] = 0.25*f_i,j-1,k + 0.5*f_i,j,k + 0.25*f_i,j+1,k 
    !      1D_z_Filter[ f_i,j,k ] = 0.25*f_i,j,k-1 + 0.5*f_i,j,k + 0.25*f_i,j,k+1 
    
    ! repeat 'N_filter' times. (Note: In the limit N_filter -> infinity, the binomial filter becomes a gaussian filter)
    DO l = 1 , N_filter
       
        ! x-filter pass
        buffer = f

        DO k = 0, nz
            DO j = 0, ny
                DO i = 0, nx
                    f(i,j,k) = 0.015625* ( &
                               buffer(i-1,j-1,k-1) + 2.0*buffer(i,j-1,k-1) + buffer(i+1,j-1,k-1)    + &
                               2.d0*(buffer(i-1,j,k-1) + 2.0*buffer(i,j,k-1) + buffer(i+1,j,k-1))   + &
                               buffer(i-1,j+1,k-1) + 2.0*buffer(i,j+1,k-1) + buffer(i+1,j+1,k-1)    + &
                               
                               2.d0 * ( buffer(i-1,j-1,k) + 2.0*buffer(i,j-1,k) + buffer(i+1,j-1,k) + &
                               2.d0*( buffer(i-1,j,k) + 2.0*buffer(i,j,k) + buffer(i+1,j,k))        + &
                               buffer(i-1,j+1,k) + 2.0*buffer(i,j+1,k) + buffer(i+1,j+1,k) )        + &
                               
                               buffer(i-1,j-1,k+1) + 2.0*buffer(i,j-1,k+1) + buffer(i+1,j-1,k+1)    + &
                               2.d0*(buffer(i-1,j,k+1) + 2.0*buffer(i,j,k+1) + buffer(i+1,j,k+1))   + &
                               buffer(i-1,j+1,k+1) + 2.0*buffer(i,j+1,k+1) + buffer(i+1,j+1,k+1)  )                     
                END DO
            END DO        
        END DO    



    END DO
    
    
END SUBROUTINE filter


! enforces particle boundary conditions
SUBROUTINE particle_boundary()

    INTEGER :: i, counter
    
    ! reset escape particles counter
    Np_esc_xm = 0
    Np_esc_xp = 0
    Np_esc_ym = 0
    Np_esc_yp = 0
    Np_esc_zm = 0
    Np_esc_zp = 0
    
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
                    particles(i,counter)%x = particles(i,counter)%x + Lx
                END IF
                IF(particles(i,counter)%x .GT. xmax) THEN
                    particles(i,counter)%x = particles(i,counter)%x - Lx
                END IF
                IF(particles(i,counter)%y .LT. ymin) THEN
                    particles(i,counter)%y = particles(i,counter)%y + Ly
                END IF
                IF(particles(i,counter)%y .GT. ymax) THEN
                    particles(i,counter)%y = particles(i,counter)%y - Ly
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

                    particles(i,counter)%oob = .TRUE.  
                     
                    !PRINT*,'Species, counter, x, y, z = ',i,counter,particles(i,counter)%x,particles(i,counter)%y

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



SUBROUTINE inject_particles()

    REAL*8 :: x, y, z, uxe, uye, uze, uxi, uyi, uzi, p(3), p1, p2, p3, &
              gam
    INTEGER :: i, j
    
    !*****************************************
    ! inject particle pairs through x- face
    !*****************************************            
    DO j = 2, ny-2, 2     
        DO i = 2, nz-2, 2 
           
   
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = 2 + p(1) 
            y = j + 2.0*p(2) 
            z = i + 2.0*p(3) 
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            uxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5
            uye = vth_e*(p(1) + p(2) + p(3))

            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5                
            uze = vth_e*(p(1) + p(2) + p(3)) 
               
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
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
            
            
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uyi = vth_i*(p(1) + p(2) + p(3))

            CALL RANDOM_NUMBER(p)                    
            p = p - 0.5                  
            uzi = vth_i*(p(1) + p(2) + p(3))

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
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            

            
        END DO
    END DO

    !*********************************************************************
    ! inject particles through z- face (roughly 4 electrons for every ion)
    !*********************************************************************
    
    ! first electrons
    DO j = 4, ny-4, 4         
        DO i = 4, nx-4, 4
            
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 4.0*p(1)
            y = j + 4.0*p(2) 
            z = 2 + p(3)
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            uxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5
            uye = vth_e*(p(1) + p(2) + p(3))

            CALL RANDOM_NUMBER(p)                                       
            uze = vth_e*table(1+64*INT(p(1))) 
               
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
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = 4, ny-4, 4 
        DO i = 8, nx-8, 8
                
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 8.0*p(1)
            y = j + 4.0*p(2) 
            z = 2 + p(3) 
        
            
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uyi = vth_i*(p(1) + p(2) + p(3))

            CALL RANDOM_NUMBER(p)                    
            uzi = vth_i*table(1+64*INT(p(1))) 

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
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            
                    
        END DO
    END DO    
    
    !*********************************************************************
    ! inject particles through z+ face (roughly 4 electrons for every ion)
    !*********************************************************************
    
    ! first electrons
    DO j = 4, ny-4, 4 
        DO i = 4, nx-4, 4

            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 4.0*p(1)
            y = j + 4.0*p(2) 
            z = nz -1  + p(3)
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            uxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5
            uye = vth_e*(p(1) + p(2) + p(3))

            CALL RANDOM_NUMBER(p)                                       
            uze = -vth_e*table(1+64*INT(p(1))) 
               
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
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = 4, ny-4, 4
        DO i = 8, nx-8, 8
       
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 8.0*p(1)
            y = j + 4.0*p(2) 
            z = nz -1 + p(3) 
        
            
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uyi = vth_i*(p(1) + p(2) + p(3))

            CALL RANDOM_NUMBER(p)                    
            uzi = -vth_i*table(1+64*INT(p(1))) 

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
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            
                    
        END DO
    END DO    
    
    
    !*********************************************************************
    ! inject particles through y- face (roughly 4 electrons for every ion)
    !*********************************************************************
    
    ! first electrons
    DO j = 4, nz-4, 4
        DO i = 4, nx-4, 4

            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 4.0*p(1)
            y = 2 + p(2)
            z = j + 4.0*p(3) 
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            uxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                       
            uye = vth_e*table(1+64*INT(p(1))) 

            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5                
            uze = vth_e*((p(1) + p(2) + p(3))) 
               
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
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = 4, nz-4, 4 
        DO i = 8, nx-8, 8
   
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 8.0*p(1)
            y = 2 + p(2)
            z = j + 4.0*p(3) 
        
            
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                      
            uyi = vth_i*table(1+64*INT(p(1))) 

            CALL RANDOM_NUMBER(p)                    
            p = p - 0.5                  
            uzi = vth_i*(p(1) + p(2) + p(3))

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
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            
                    
        END DO
    END DO    
    
    !*********************************************************************
    ! inject particles through y+ face (roughly 4 electrons for every ion)
    !*********************************************************************
    
    ! first electrons
    DO j = 4, nz-4, 4 
        DO i = 4, nx-4, 4
  
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 4.0*p(1)
            y = ny - 1 + p(2)
            z = j + 4.0*p(3)  
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            uxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                       
            uye = -vth_e*table(1+64*INT(p(1))) 

            CALL RANDOM_NUMBER(p)                                       
            p = p - 0.5                
            uze = vth_e*((p(1) + p(2) + p(3))) 
               
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
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = 4, nz-4, 4 
        DO i = 8, nx-8, 8
    
            CALL RANDOM_NUMBER(p)                                        
            p = p - 0.5
            x = i + 8.0*p(1)
            y = ny - 1 + p(2)
            z = j + 4.0*p(3) 
        
            
            CALL RANDOM_NUMBER(p)                                      
            p = p - 0.5
            uxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RANDOM_NUMBER(p)                                      
            uyi = -vth_i*table(1+64*INT(p(1))) 

            CALL RANDOM_NUMBER(p)                    
            p = p - 0.5                  
            uzi = vth_i*(p(1) + p(2) + p(3))

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
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            
                    
        END DO
    END DO    
     
    
    PRINT*,''
    PRINT*,'Np_in(electrons) = ',Np_in(1)
    PRINT*,'Np_in(ions) = ',Np_in(2)
    PRINT*,''

END SUBROUTINE inject_particles



! particle deletion involvles re-sorting: the deleted particle is exchanged with the last particle on the list
SUBROUTINE delete_particle(i,j)

    INTEGER, INTENT(IN) :: i, j
    REAL*8 :: temp(6)
    ! replace this particle with the last particle in our array, and decrement Np_in    
    temp(1) = particles(i,j)%x    
    temp(2) = particles(i,j)%y
    temp(3) = particles(i,j)%z
    temp(4) = particles(i,j)%ux
    temp(5) = particles(i,j)%uy
    temp(6) = particles(i,j)%uz
    
    particles(i,j)%x = particles(i,Np_in(i))%x    
    particles(i,j)%y = particles(i,Np_in(i))%y
    particles(i,j)%z = particles(i,Np_in(i))%z
    particles(i,j)%ux = particles(i,Np_in(i))%ux
    particles(i,j)%uy = particles(i,Np_in(i))%uy
    particles(i,j)%uz = particles(i,Np_in(i))%uz

    
    ! copy the deleted particle into the appropriate boundary particles array

    ! x- boundary
    IF(temp(1) .LT. xmin)THEN
 
        Np_esc_xm = Np_esc_xm + 1 
         
        bparticles_xm(Np_esc_xm)%species = i    
        bparticles_xm(Np_esc_xm)%x = temp(1)    
        bparticles_xm(Np_esc_xm)%y = temp(2)
        bparticles_xm(Np_esc_xm)%z = temp(3)
        bparticles_xm(Np_esc_xm)%ux = temp(4)
        bparticles_xm(Np_esc_xm)%uy = temp(5)
        bparticles_xm(Np_esc_xm)%uz = temp(6)

    ! x+ boundary
    ELSE IF(temp(1) .GT. xmax)THEN

        Np_esc_xp = Np_esc_xp + 1 

        bparticles_xp(Np_esc_xp)%species = i    
        bparticles_xp(Np_esc_xp)%x = temp(1)    
        bparticles_xp(Np_esc_xp)%y = temp(2)
        bparticles_xp(Np_esc_xp)%z = temp(3)
        bparticles_xp(Np_esc_xp)%ux = temp(4)
        bparticles_xp(Np_esc_xp)%uy = temp(5)
        bparticles_xp(Np_esc_xp)%uz = temp(6)

         
         
    ! y- boundary        
    ELSE IF(temp(2) .LT. ymin)THEN

        Np_esc_ym = Np_esc_ym + 1 
         
        bparticles_ym(Np_esc_ym)%species = i    
        bparticles_ym(Np_esc_ym)%x = temp(1)    
        bparticles_ym(Np_esc_ym)%y = temp(2)
        bparticles_ym(Np_esc_ym)%z = temp(3)
        bparticles_ym(Np_esc_ym)%ux = temp(4)
        bparticles_ym(Np_esc_ym)%uy = temp(5)
        bparticles_ym(Np_esc_ym)%uz = temp(6)

    ! y+ boundary
    ELSE IF(temp(2) .GT. ymax)THEN
        
        Np_esc_yp = Np_esc_yp + 1 
         
        bparticles_yp(Np_esc_yp)%species = i    
        bparticles_yp(Np_esc_yp)%x = temp(1)    
        bparticles_yp(Np_esc_yp)%y = temp(2)
        bparticles_yp(Np_esc_yp)%z = temp(3)
        bparticles_yp(Np_esc_yp)%ux = temp(4)
        bparticles_yp(Np_esc_yp)%uy = temp(5)
        bparticles_yp(Np_esc_yp)%uz = temp(6)

    ! z- boundary
    ELSE IF(temp(3) .LT. zmin) THEN
    
        Np_esc_zm = Np_esc_zm + 1 
         
        bparticles_zm(Np_esc_zm)%species = i    
        bparticles_zm(Np_esc_zm)%x = temp(1)    
        bparticles_zm(Np_esc_zm)%y = temp(2)
        bparticles_zm(Np_esc_zm)%z = temp(3)
        bparticles_zm(Np_esc_zm)%ux = temp(4)
        bparticles_zm(Np_esc_zm)%uy = temp(5)
        bparticles_zm(Np_esc_zm)%uz = temp(6)
    
    ! z+ boundary    
    ELSE IF(temp(3) .GT. zmax) THEN
     
        Np_esc_zp = Np_esc_zp + 1 
         
        bparticles_zp(Np_esc_zp)%species = i    
        bparticles_zp(Np_esc_zp)%x = temp(1)    
        bparticles_zp(Np_esc_zp)%y = temp(2)
        bparticles_zp(Np_esc_zp)%z = temp(3)
        bparticles_zp(Np_esc_zp)%ux = temp(4)
        bparticles_zp(Np_esc_zp)%uy = temp(5)
        bparticles_zp(Np_esc_zp)%uz = temp(6)

    END IF                
                            
    Np_in(i) = Np_in(i) - 1

END SUBROUTINE delete_particle



SUBROUTINE create_particle(sp,x,y,z,ux,uy,uz)

    INTEGER, INTENT(IN) :: sp
    REAL*8, INTENT(IN) :: x,y,z,ux,uy,uz

    ! insert new particle at the end of our array, and increment Np_in    
    Np_in(sp) = Np_in(sp) + 1
   
    particles(sp,Np_in(sp))%x = x   
    particles(sp,Np_in(sp))%y = y
    particles(sp,Np_in(sp))%z = z

    particles(sp,Np_in(sp))%ux = ux   
    particles(sp,Np_in(sp))%uy = uy
    particles(sp,Np_in(sp))%uz = uz
    
    particles(sp,Np_in(sp))%oob = .FALSE.
    particles(sp,Np_in(sp))%species = sp
   
    IF( x .LT. xmin .OR. x .GT. xmax .OR. y .LT. ymin .OR. y .GT. ymax .OR. &
        z .LT. zmin .OR. z .GT. zmax) THEN
      
        PRINT*,'ERROR...out of bound particle creation.'
        PRINT*,'x, y, z =',x,y,z
        STOP      
                    
    END IF                 
    
END SUBROUTINE create_particle



SUBROUTINE add_ring_current(ts)

    INTEGER, INTENT(IN) :: ts
    REAL*8 :: x, y, r, Rc
    INTEGER :: n

    
    ! radius of current loop
    Rc = 1.0


    !gradually ramp up ring current
    IF(ABS(o1) .GT. 0.d0) THEN
        o2 = o2 + o3
        o1 = o1 + o2
        o = o + o1
    END IF
    
    PRINT*,'Ring Current = ',o
    
    
    
    ! sinusiodal ring current
    o = SIN(twopi*ts/10.d0)
    
    
    ! decrement smoothed ring current from electric field
    IF(current_filter_on) THEN
        DO n = 1, 27
            Ex(ms(n)+ie,je,ke) = Ex(ms(n)+ie,je,ke) - (sm(n)*o)
            Ex(ms(n)+ie,je+1,ke) = Ex(ms(n)+ie,je+1,ke) + (sm(n)*o)
            Ey(ms(n)+ie,je,ke) = Ey(ms(n)+ie,je,ke) + (sm(n)*o)
            Ey(ms(n)+ie+1,je,ke) = Ey(ms(n)+ie+1,je,ke) - (sm(n)*o)
        END DO
   ELSE
        n = 14
        Ex(ms(n)+ie,je,ke) = Ex(ms(n)+ie,je,ke) - 8.d0*(sm(n)*o)
        Ex(ms(n)+ie,je+1,ke) = Ex(ms(n)+ie,je+1,ke) + 8.d0*(sm(n)*o)
        Ey(ms(n)+ie,je,ke) = Ey(ms(n)+ie,je,ke) + 8.d0*(sm(n)*o)
        Ey(ms(n)+ie+1,je,ke) = Ey(ms(n)+ie+1,je,ke) - 8.d0*(sm(n)*o)
    END IF         
    
    
END SUBROUTINE add_ring_current


END MODULE particleMover_mod