MODULE particleMover_mod

USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


SUBROUTINE accel(x, y, qm, a, b)

	REAL*8, INTENT(IN) :: x, y, qm
	REAL*8, INTENT(INOUT) :: a(2), b
	INTEGER :: ix, iy

	IF(interpolation_type .EQ. 1 .OR. interpolation_type .EQ. 3) THEN
		ix = x + 0.5
        iy = y + 0.5
        
		a(1) = qm*Ex_grid(ix,iy)
		a(2) = qm*Ey_grid(ix,iy)
        b = qm*Bz_grid(ix,iy)
        
	ELSE IF(interpolation_type .EQ. 2) THEN
	
		ix = x 
        iy = y
        
		a(1) = qm*( (1.d0-(x-ix))*(1.d0-(y-iy))*Ex_grid(ix,iy) + (x-ix)*(1.d0-(y-iy))*Ex_grid(ix+1,iy) + &
                    (1.d0-(x-ix))*(y-iy)*Ex_grid(ix,iy+1) + (x-ix)*(y-iy)*Ex_grid(ix+1,iy+1)  )
                    
		a(2) = qm*( (1.d0-(x-ix))*(1.d0-(y-iy))*Ey_grid(ix,iy) + (x-ix)*(1.d0-(y-iy))*Ey_grid(ix+1,iy) + &
                    (1.d0-(x-ix))*(y-iy)*Ey_grid(ix,iy+1) + (x-ix)*(y-iy)*Ey_grid(ix+1,iy+1)  )
                    
        b    = qm*( (1.d0-(x-ix))*(1.d0-(y-iy))*Bz_grid(ix,iy) + (x-ix)*(1.d0-(y-iy))*Bz_grid(ix+1,iy) + &
                    (1.d0-(x-ix))*(y-iy)*Bz_grid(ix,iy+1) + (x-ix)*(y-iy)*Bz_grid(ix+1,iy+1)  )        
        
	END IF

END SUBROUTINE accel



! Deposits currents to grid-cell faces (charge conservative) using Villasenor-Buneman (1992) algorithm
! Jx(j,k) = Jx(x_j+1/2, y_k)
! Jy(j,k) = Jy(x_j, y_k+1/2)

SUBROUTINE current_deposit(x1,y1,x2,y2,q)

    REAL*8, INTENT(IN) :: x1, y1 ! initial position of macro particle
    REAL*8, INTENT(IN) :: x2, y2 ! final position of macro particle
    REAL*8, INTENT(IN) :: q
    
    INTEGER :: i, j
    REAL*8 :: xx, yy, v1_i, v1_j, v2_i, v2_j, vint_i, vint_j  
    REAL*8 :: Jx1, Jy1, Jx2, Jy2
    REAL*8 :: Jx1p, Jy1p, Jx2p, Jy2p
    REAL*8 :: Jx1pp, Jy1pp, Jx2pp, Jy2pp
    REAL*8 :: delx, dely, delx1, dely1, delx2, dely2,  delx3, dely3
    REAL*8 :: x1p, y1p
    REAL*8 :: x1pp, y1pp

    !PRINT*,''
    !PRINT*,'Current Deposition:'
    !PRINT*,'x1,y1=',x1,y1
    !PRINT*,'x2,y2=',x2,y2
    
    
    ! Find nearest vertex location at initial position       
    xx = x1/dx
    yy = y1/dy
    i = 0.5 + xx  ! occupying cell indices
    j = 0.5 + yy

    !PRINT*,'Initial cell index = ',i,j
    
    v1_i = i + 0.5 * SIGN(1.d0,xx-i)
    v1_j = j + 0.5 * SIGN(1.d0,yy-j)
   
    !PRINT*,'Nearest Cell vertex index =',v1_i,v1_j
    
    ! Find nearest vertex location at final position
    xx = x2/dx
    yy = y2/dy
    i = 0.5 + xx  ! occupying cell indices
    j = 0.5 + yy

    !PRINT*,'Final cell index = ',i,j
   
   
    v2_i = i + 0.5d0 * SIGN(1.d0,xx-i)
    v2_j = j + 0.5d0 * SIGN(1.d0,yy-j)
    
    !PRINT*,'Nearest Cell vertex index =',v2_i,v2_j
    !PRINT*,''
  
    delx = x2 - x1
    dely = y2 - y1
    
    !****************************************************************************
    ! Case 1: 4 boundary move
    ! Occurs when nearest cell vertex is the same for initial and final positions.
    !****************************************************************************
    IF(v1_i .EQ. v2_i .AND. v1_j .EQ. v2_j) THEN
    
        !PRINT*,'Case 1. 4 boundary move.'
    
        ! Compute flues through the 4 boundaries
   
        ! Jx1, Jx2 are fluxes through the bottom and top edges of the vertex 
        ! Jy1, Jy2 are fluxes through the left and right edges of the vertex
        ! (charge * area of parallelogram swept out by each boundary face as the square macro particle moves
        ! along a straight line connecting the initial and final positions)
        Jx1 = q*dy*( v1_j + 0.5 - (y1 + 0.5*dely)/dy )*delx       
        Jx2 = q*dy*(-v1_j + 0.5 + (y1 + 0.5*dely)/dy )*delx
        
        Jy1 = q*dx*( v1_i + 0.5 - (x1 + 0.5*delx)/dx )*dely       
        Jy2 = q*dx*(-v1_i + 0.5 + (x1 + 0.5*delx)/dx )*dely
        
       !PRINT*,'Jx1, Jx2 = ', Jx1, Jx2
       !PRINT*,'Jy1, Jy2 = ', Jy1, Jy2

       !*******************       
       ! Add to flux arrays
       !*******************
       i = v1_i
       j = v1_j
       Jx(i,j)   = Jx1
       Jy(i,j)   = Jy1
       Jx(i,j+1) = Jx2
       Jy(i+1,j) = Jy2
       
      
    !****************************************************************************
    ! Case 2: 7 boundary move (Particle cuts through two cell vertices)
    ! Occurs when either the nearest vertex gets displaced by either 1 unit in
    ! the x direction or  1 unit in y direction, but not both.    
    !****************************************************************************
    ELSE IF( (ABS(v1_i-v2_i) .GT. 0.01d0 .AND. ABS(v1_j-v2_j) .LT. 0.99d0) &
        .OR. (ABS(v1_j-v2_j) .GT. 0.01d0 .AND. ABS(v1_i-v2_i) .LT. 0.99d0) ) THEN
       
        !PRINT*,'Case2. 7 boundary move'     
   
        ! Can decompose the displacements into two pieces, this becomes two successive 4 boundary moves
        ! In the first move, we compute fluxes Jx1, Jx2, Jy1, Jy2 through the four edges around the initial vertex
        ! In the second move, we compute fluxes Jx1' , Jx2', Jy1' , Jy2' through the four edges around the final vertex
        ! Break up the displacements delx = delx1 + dex2 and dely = dely1 + dely2
        ! Then use (delx1,dely1) for the first move and (delx2,dely2) for the second move. 
   
        !PRINT*,'delx,dely=',delx,dely
        !PRINT*,''
        
        ! determine the displcement direction of the final vertex
        IF( (v2_i-v1_i) .GT. 0.01d0 ) THEN ! right displacement of vertex
            !PRINT*,'Right displacement of vertex.'
            delx1 = dx*v1_i-(x1-0.5*dx)
            dely1 = ((y2-y1)/(x2-x1))*delx1
            
        ELSE IF( (v2_i-v1_i) .LT. -0.01d0) THEN ! left displacement
            !PRINT*,'Left displacement of vertex.'
            delx1 = -(x1+0.5*dx - dx*v1_i) 
            dely1 = ((y2-y1)/(x2-x1))*delx1
            
        ELSE IF( (v2_j-v1_j) .GT. 0.01d0 )THEN ! upward displacement 
            !PRINT*,'Up displacement of vertex.'
            dely1 = dy*v1_j-(y1-0.5*dy)
            delx1 = ((x2-x1)/(y2-y1))*dely1
            
        ELSE IF( (v2_j-v1_j) .LT. -0.01d0) THEN ! downward displacement
            !PRINT*,'Down displacement of vertex.'
            dely1 = -(y1+0.5*dy - dy*v1_j) 
            delx1 = ((x2-x1)/(y2-y1))*dely1
        END IF    
        !PRINT*,''
        
        !*******
        ! Move 1
        !*******
        Jx1 = q*dy*( v1_j + 0.5 - (y1 + 0.5*dely1)/dy )*delx1       
        Jx2 = q*dy*(-v1_j + 0.5 + (y1 + 0.5*dely1)/dy )*delx1
        
        Jy1 = q*dx*( v1_i + 0.5 - (x1 + 0.5*delx1)/dx )*dely1       
        Jy2 = q*dx*(-v1_i + 0.5 + (x1 + 0.5*delx1)/dx )*dely1
   
                
        !*******
        ! Move 2
        !*******     
        x1p = x1 + delx1 ! updated positions after Move 1
        y1p = y1 + dely1
        
        delx2 = delx - delx1
        dely2 = dely - dely1
               
        Jx1p = q*dy*( v2_j + 0.5 - (y1p + 0.5*dely2)/dy )*delx2       
        Jx2p = q*dy*(-v2_j + 0.5 + (y1p + 0.5*dely2)/dy )*delx2
        
        Jy1p = q*dx*( v2_i + 0.5 - (x1p + 0.5*delx2)/dx )*dely2       
        Jy2p = q*dx*(-v2_i + 0.5 + (x1p + 0.5*delx2)/dx )*dely2

   
      	IF(print_debug) THEN

        PRINT*,'delx1,delx2=',delx1,delx2
        PRINT*,'dely1,dely2=',dely1,dely2
        PRINT*,'Jx1, Jx2 = ', Jx1, Jx2
        PRINT*,'Jy1, Jy2 = ', Jy1, Jy2
        PRINT*,'Jx1p, Jx2p = ', Jx1p, Jx2p
        PRINT*,'Jy1p, Jy2p = ', Jy1p, Jy2p       
        PRINT*,'sum(Jx) = ',Jx1+Jx2+Jx1p+Jx2p
        PRINT*,'sum(Jy) = ',Jy1+Jy2+Jy1p+Jy2p
        PRINT*,'delx1+delx2, dely1+dely2=',delx1+delx2, dely1+dely2
        
        END IF
        
       !*******************       
       ! Add to flux arrays
       !*******************
       i = v1_i
       j = v1_j
       Jx(i,j)   = Jx1
       Jy(i,j)   = Jy1
       Jx(i,j+1) = Jx2
       Jy(i+1,j) = Jy2
       
       i = v2_i
       j = v2_j
       Jx(i,j)   = Jx1p
       Jy(i,j)   = Jy1p
       Jx(i,j+1) = Jx2p
       Jy(i+1,j) = Jy2p       
       
    !****************************************************************************
    ! Case 3: 10 boundary move (particle cuts through three cell vertices)
    ! Occurs when the nearest vertex gets displaced by 1 unit in the x direction 
    ! and 1 unit in y direction.
    !**************************************************************************** 
    ELSE IF(ABS(v1_i-v2_i) .GT. 0.01d0 .AND. ABS(v1_j-v2_j) .GT. 0.01d0 ) THEN


        ! Similar to the 7 boundary move case, can decompose the displacements into three pieces, i.e three 
        ! successive 4 boundary moves. In the first move, we compute fluxes Jx1, Jx2, Jy1, Jy2 through the four edges around the initial vertex.
        ! In the second move, we compute fluxes Jx1' , Jx2', Jy1' , Jy2' through the four edges around the intermediate vertex (the final vertex is
        ! on the along the diagonal and the intermediate vertex is one of the two remaining vertices in the square containing the initial and final 
        ! vertices, which intermediate vertex it cuts through depends on the initial position of the particle). In the third move, we compute 
        ! fluxes Jx1'' , Jx2'', Jy1'' , Jy2'' through the four edges around the final vertex.
        ! Break up the displacements delx = delx1 + dex2 + delx3 and dely = dely1 + dely2 + dely3
        ! Then use (delx1,dely1, delz1) for the first move, (delx2,dely2,delz2) for the second and (delx3,dely3,delz3) for the third. 

        !PRINT*,'Case3. 10 boundary move'     

        ! Determine intermediate vertex location. (8 possiblities)
        ! Easily done by checking which side of the diagonal line, connecting the intial and final vertices, the
        ! particle initially lies on.
        
        IF( (v2_i-v1_i) .GT. 0.01d0 .AND. (v2_j-v1_j) .GT. 0.01d0) THEN  ! top-right final vertex      
           
            !PRINT*,'Top-right final vertex.'
            
            IF( (y1-dy*v1_j) .GT. (x1-dx*v1_i) ) THEN ! particle above diagonal
                !PRINT*,'Particle above diagonal.'
                vint_i = v1_i              
                vint_j = v1_j + 1
                
                dely1 = dy*v1_j-(y1-0.5*dy)
                delx1 = ((x2-x1)/(y2-y1))*dely1
                
                x1p = x1 + delx1
                y1p = y1 + dely1
                
                delx2 = dx*vint_i-(x1p-0.5*dx)
                dely2 = ((y2-y1)/(x2-x1))*delx2   ! Note how the displacements reflect the directions 
                                                   ! going from initial to intermediate to final vertex

                
            ELSE ! particle below diagonal
                !PRINT*,'Particle below diagonal.'
                vint_i = v1_i + 1             
                vint_j = v1_j 
                
                delx1 = dx*v1_i-(x1-0.5*dx)
                dely1 = ((y2-y1)/(x2-x1))*delx1
            
                x1p = x1 + delx1
                y1p = y1 + dely1
                           
                dely2 = dy*vint_j-(y1p-0.5*dy)
                delx2 = ((x2-x1)/(y2-y1))*dely2   ! Note how the displacements reflect the directions 
                                                   ! going from initial to intermediate to final vertex
            
            END IF

        ELSE IF( (v2_i-v1_i) .LT. -0.01d0 .AND. (v2_j-v1_j) .LT. -0.01d0) THEN ! bottom-left final vertex      
           
            !PRINT*,'Bottom-left final vertex.'
            
            IF( (y1-dy*v1_j) .GT. (x1-dx*v1_i) ) THEN ! particle above diagonal
                !PRINT*,'Particle above diagonal.'
                vint_i = v1_i - 1             
                vint_j = v1_j 
                
                delx1 = -(x1+0.5*dx - dx*v1_i) 
                dely1 = ((y2-y1)/(x2-x1))*delx1
                
                x1p = x1 + delx1
                y1p = y1 + dely1
                
                dely2 = -(y1p+0.5*dy - dy*vint_j) 
                delx2 = ((x2-x1)/(y2-y1))*dely2
                
            ELSE ! particle below diagonal
                !PRINT*,'Particle below diagonal.'
                vint_i = v1_i              
                vint_j = v1_j - 1
                
                dely1 = -(y1+0.5*dy - dy*v1_j) 
                delx1 = ((x2-x1)/(y2-y1))*dely1
            
                x1p = x1 + delx1
                y1p = y1 + dely1
                
                delx2 = -(x1p+0.5*dx - dx*vint_i) 
                dely2 = ((y2-y1)/(x2-x1))*delx2
            
            END IF
 
        ELSE IF( (v2_i-v1_i) .LT. -0.01d0 .AND. (v2_j-v1_j) .GT. 0.01d0) THEN ! top-left final vertex      
           
            !PRINT*,'Top-left final vertex.'
            
            IF( (y1-dy*v1_j) .GT. -(x1-dx*v1_i) ) THEN ! particle above diagonal
                !PRINT*,'Particle above diagonal.'
                vint_i = v1_i              
                vint_j = v1_j + 1
                
                dely1 = dy*v1_j-(y1-0.5*dy)
                delx1 = ((x2-x1)/(y2-y1))*dely1
                
                x1p = x1 + delx1
                y1p = y1 + dely1     
                
                delx2 = -(x1p+0.5*dx - dx*vint_i) 
                dely2 = ((y2-y1)/(x2-x1))*delx2
                               
                
            ELSE ! particle below diagonal
                !PRINT*,'Particle below diagonal.'
                vint_i = v1_i - 1             
                vint_j = v1_j 
            
                delx1 = -(x1+0.5*dx - dx*v1_i) 
                dely1 = ((y2-y1)/(x2-x1))*delx1
            
                x1p = x1 + delx1
                y1p = y1 + dely1 
                
                dely2 = dy*vint_j-(y1p-0.5*dy)
                delx2 = ((x2-x1)/(y2-y1))*dely2
                
            END IF

        ELSE IF( (v2_i-v1_i) .GT. 0.01d0 .AND. (v2_j-v1_j) .LT. -0.01d0) THEN  ! bottom-right final vertex      
           
            !PRINT*,'Bottom-right final vertex.'
            
            IF( (y1-dy*v1_j) .GT. -(x1-dx*v1_i) ) THEN ! particle above diagonal
                !PRINT*,'Particle above diagonal.'
                vint_i = v1_i + 1             
                vint_j = v1_j 
                
                delx1 = dx*v1_i-(x1-0.5*dx)
                dely1 = ((y2-y1)/(x2-x1))*delx1
                
                x1p = x1 + delx1
                y1p = y1 + dely1
                
                dely2 = -(y1p+0.5*dy - dy*vint_j) 
                delx2 = ((x2-x1)/(y2-y1))*dely2
            
            ELSE ! particle below diagonal
                !PRINT*,'Particle below diagonal.'
                vint_i = v1_i              
                vint_j = v1_j - 1
            
                dely1 = -(y1+0.5*dy - dy*v1_j) 
                delx1 = ((x2-x1)/(y2-y1))*dely1
            
                x1p = x1 + delx1
                y1p = y1 + dely1
            
                delx2 = dx*vint_i-(x1p-0.5*dx)
                dely2 = ((y2-y1)/(x2-x1))*delx2
                        
            END IF

        END IF

        !PRINT*,''
        !PRINT*,'Intermediate Cell vertex index =',vint_i,vint_j
        !PRINT*,''
        
        x1pp = x1p + delx2 ! update positions after Move 2
        y1pp = y1p + dely2
        
        delx3 = delx - delx1 - delx2
        dely3 = dely - dely1 - dely2
    
    
        !*******
        ! Move 1
        !*******
        Jx1 = q*dy*( v1_j + 0.5 - (y1 + 0.5*dely1)/dy )*delx1       
        Jx2 = q*dy*(-v1_j + 0.5 + (y1 + 0.5*dely1)/dy )*delx1
        
        Jy1 = q*dx*( v1_i + 0.5 - (x1 + 0.5*delx1)/dx )*dely1       
        Jy2 = q*dx*(-v1_i + 0.5 + (x1 + 0.5*delx1)/dx )*dely1


        !*******
        ! Move 2
        !*******
        Jx1p = q*dy*( vint_j + 0.5 - (y1p + 0.5*dely2)/dy )*delx2       
        Jx2p = q*dy*(-vint_j + 0.5 + (y1p + 0.5*dely2)/dy )*delx2
        
        Jy1p = q*dx*( vint_i + 0.5 - (x1p + 0.5*delx2)/dx )*dely2       
        Jy2p = q*dx*(-vint_i + 0.5 + (x1p + 0.5*delx2)/dx )*dely2


        !*******
        ! Move 3
        !*******
        Jx1pp = q*dy*( v2_j + 0.5 - (y1pp + 0.5*dely3)/dy )*delx3       
        Jx2pp = q*dy*(-v2_j + 0.5 + (y1pp + 0.5*dely3)/dy )*delx3
        
        Jy1pp = q*dx*( v2_i + 0.5 - (x1pp + 0.5*delx3)/dx )*dely3      
        Jy2pp = q*dx*(-v2_i + 0.5 + (x1pp + 0.5*delx3)/dx )*dely3
        
      
        IF(print_debug) THEN

            PRINT*,'delx1,delx2=',delx1,delx2
            PRINT*,'dely1,dely2=',dely1,dely2
            PRINT*,'dely3,dely3=',dely3,dely3
            PRINT*,''
            PRINT*,'Jx1, Jx2 = ', Jx1, Jx2
            PRINT*,'Jy1, Jy2 = ', Jy1, Jy2
            PRINT*,'Jx1p, Jx2p = ', Jx1p, Jx2p
            PRINT*,'Jy1p, Jy2p = ', Jy1p, Jy2p
            PRINT*,'Jx1pp, Jx2pp = ', Jx1pp, Jx2pp
            PRINT*,'Jy1pp, Jy2pp = ', Jy1pp, Jy2pp       
            PRINT*,''
            PRINT*,'sum(Jx) = ',Jx1+Jx2+Jx1p+Jx2p+Jx1pp+Jx2pp
            PRINT*,'sum(Jy) = ',Jy1+Jy2+Jy1p+Jy2p+Jy1pp+Jy2pp
            PRINT*,'delx1+delx2+delx3=',delx1+delx2+delx3
            PRINT*,'dely1+dely2+dely3=',dely1+dely2+dely3
        
        END IF
        
        !*******************       
        ! Add to flux arrays
        !*******************
        i = v1_i
        j = v1_j
        
        Jx(i,j)   = Jx1
        Jy(i,j)   = Jy1
        Jx(i,j+1) = Jx2
        Jy(i+1,j) = Jy2
       
        i = vint_i
        j = vint_j
        Jx(i,j)   = Jx1p
        Jy(i,j)   = Jy1p
        Jx(i,j+1) = Jx2p
        Jy(i+1,j) = Jy2p

        i = v2_i
        j = v2_j
        Jx(i,j)   = Jx1pp
        Jy(i,j)   = Jy1pp
        Jx(i,j+1) = Jx2pp
        Jy(i+1,j) = Jy2pp
         
    END IF

    !PRINT*,''

END SUBROUTINE


! updates velocities and positions of particles (using leap-frog integration) and acculmulates charge at the grid-points
SUBROUTINE move_particles(ts)

    INTEGER, INTENT(IN) :: ts ! timestep index
	INTEGER :: i, j, k, ix, iy
	REAL*8 :: qm, wc, x, y, a(2), qdxy, temp, uxx
    REAL*8 :: b, t, s, gam
    REAL*8 :: x1, y1, x2, y2
		
	! clear kinitic energy and momentum accumulator	
	KE = 0.d0
	Px = 0.d0
	Py = 0.d0
	
    ! clear density
    rho = 0.d0
    rho(1:nx,1:ny) = rho0
	ne(:,:) = 0.d0
	ni(:,:) = 0.d0
	
    ! clear currents
    Jx = 0.d0
    Jy = 0.d0
 
        ! particle update loop 
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

                ! only update particles that are inside the domain
                IF(.NOT. particles(i,j)%oob .OR. bndry .EQ. 1) THEN

                x = particles(i,j)%x/dx
                y = particles(i,j)%y/dy
            
                x1 = particles(i,j)%x           
                y1 = particles(i,j)%y         
            
                !******************************
                ! update u(t-dt/2) -> u(t+dt/2)
                !******************************
               ! IF(ts .LT. maxsteps) THEN 
               !     b = 0.d0
               !     a(1) = 0.05*c
               !     a(2) = 0.d0		
               ! ELSE
                    CALL accel(x, y, qm, a, b)
                !END IF
             
                ! half electric push (u^n-1/2 ->u^-) 
                particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*dt
                particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*dt
						
                gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2)/c**2 )
                t = 0.5d0*b*dt/(gam*c)
                s = 2.d0*t/(1.d0+t**2)
                uxx = particles(i,j)%ux
                                
                ! magnetic rotation (u^- -> u^+)
                particles(i,j)%ux = (1.d0-s*t)*particles(i,j)%ux + s*particles(i,j)%uy
                particles(i,j)%uy = -s*uxx + (1.d0-s*t)*particles(i,j)%uy
			
            
                ! remaining half electric push (u^+ -> u^n+1/2)		
                particles(i,j)%ux = particles(i,j)%ux + 0.5*a(1)*dt
                particles(i,j)%uy = particles(i,j)%uy + 0.5*a(2)*dt
                						
                KE = KE + (gam-1.d0)* particles(i,j)%m * (c**2)
                Px = Px + particles(i,j)%m*particles(i,j)%ux
                Py = Py + particles(i,j)%m*particles(i,j)%uy
			
                !***********************
                ! update x(t) -> x(t+dt)
                !***********************
                gam = SQRT( 1.d0 + (particles(i,j)%ux**2 + particles(i,j)%uy**2)/c**2 )
                particles(i,j)%x = particles(i,j)%x + particles(i,j)%ux*dt/gam
                particles(i,j)%y = particles(i,j)%y + particles(i,j)%uy*dt/gam       
            
                !******************************************
                ! deposit currents at grid-cell interfaces 
                ! Jx_j+1/2,k and Jy_j,k+1/2
                !******************************************
                x2 = particles(i,j)%x           
                y2 = particles(i,j)%y         

                CALL current_deposit(x1,y1,x2,y2, particles(i,1)%q)
                                       
                !**********************************
                ! deposit charge at grid points
                ! rho_j,k
                !**********************************
			
                IF(interpolation_type .EQ. 1) THEN			
                
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
		
                !************************************
                ! apply particle boundary conditions
                !************************************

                CALL particle_boundary(particles(i,j))
            
                !PRINT*,'PARTICLE TYPE, ID, x, y, ux, uy = ',i,j,particles(i,j)%x, particles(i,j)%y, particles(i,j)%ux,particles(i,j)%uy 

            END IF
            
		END DO
	END DO


    ! If periodic boundaries, need to touch up the grid deposited charges and currents at boundaries
    IF (bndry .EQ. 1) THEN
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

        buffer = Jx
        Jx(0,:) = Jx(0,:) + Jx(nx,:)	
        Jx(1,:) = Jx(1,:) + Jx(nx+1,:)	
        Jx(nx,:) = Jx(nx,:) + buffer(0,:)
        Jx(nx-1,:) = Jx(nx-1,:) + Jx(-1,:)

        buffer = Jy
        Jy(:,0) = Jy(:,0) + Jy(:,ny)	
        Jy(:,1) = Jy(:,1) + Jy(:,ny+1)	
        Jy(:,ny) = Jy(:,ny) + buffer(:,0)
        Jy(:,ny-1) = Jy(:,ny-1) + Jy(:,-1)
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
    END IF
    
    
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



! 2D biomial filter (equivalent to applying a y-binomial filter followed by an x-binomial filter)
SUBROUTINE filter(f)

    REAL*8, INTENT(INOUT) ::  f(-1:nx+3,-1:ny+3)
    REAL*8, PARAMETER :: M = 4.0, S = 2.0, K = 1.0   ! filter weights 
    INTEGER :: i, j

    ! copy input into buffer
    buffer(-1:nx+1,-1:ny+1) = f(-1:nx+1,-1:ny+1)

    ! apply 2d binomial filter
    ! Filter(f_i,j) = [ M*f_i,j + S*(f_i-1,j + f_i+1,j + f_i,j-1 + f_i,j+1) + K*(f_i-1,j-1 + f_i-1,j+1 + f_i+1,j-1 + f_i+1,j+1) ]/( M + 4*(S+K) )
    
    DO i = 0, nx
        DO j = 0, ny
            f(i,j) = M*buffer(i,j) + S*( buffer(i-1,j)+buffer(i+1,j)+buffer(i,j-1)+buffer(i,j+1) ) + &
                     K*( buffer(i-1,j-1)+buffer(i-1,j+1)+buffer(i+1,j-1)+buffer(i+1,j+1) )
            f(i,j) = f(i,j)/(M+4.d0*(S+K))
        END DO            
    END DO    
	
	
END SUBROUTINE filter



! enforces particle boundary conditions
SUBROUTINE particle_boundary(p)

    TYPE(particle), INTENT(INOUT) :: p


    ! periodic boundary conditions
    IF(bndry .eq. 1) THEN
    
        IF(p%x .LT. xmin) THEN
            p%x = p%x + Lx
        END IF
        IF(p%x .GT. xmax) THEN
            p%x = p%x - Lx
        END IF
        IF(p%y .LT. ymin) THEN
            p%y = p%y + Ly
        END IF
        IF(p%y .GT. ymax) THEN
            p%y = p%y - Ly
        END IF

    END IF

    ! outflow boundary conditions (in this case, need to distribute ions and electrons uniformly across boundary
    ! to avoid net charge buildup)
    IF(bndry .EQ. 2) Then       
       ! set out of bounds flag
       p%oob = .TRUE.
     
    END IF
    

END SUBROUTINE particle_boundary



END MODULE particleMover_mod