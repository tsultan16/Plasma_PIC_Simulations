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
    INTEGER, PARAMETER :: npass = 0
	REAL*8 :: p(3), x, y, z, vxe, vye, vze, vxi, vyi, vzi
    
    

    !GO TO 111
    
    DO l = 1, npass
    
    ! distribute particles across grid (1 pair per cell)    
    DO k = 1, nz, 1 
        DO j = 1, ny, 1
            DO i = 1, nx, 1
	
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + p(1)
            y = j + p(2) 
            z = k + p(3) 
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            vxe = v0x + vth_e * (vth_e + SUM(p))
                                        
            CALL RAND_NUM(p)                                       
            p = p - 0.5
            vye = vth_e * SUM(p)

            CALL RAND_NUM(p)                                       
            p = p - 0.5                
            vze = vth_e * SUM(p)
               
                     
            CALL create_particle(1,x,y,z,vxe,vye,vze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
               
            IF(ns .EQ. 2) THEN


                CALL RAND_NUM(p)                                      

                p = p - 0.5
                


                vxi = v0x + vth_i*(vth_i + SUM(p))
                
                CALL RAND_NUM(p)                                      
                p = p - 0.5
                vyi = vth_i * SUM(p)

 

                CALL RAND_NUM(p)                    
                p = p - 0.5                  
                vzi = vth_i * SUM(p)
  
            

                CALL create_particle(2,x,y,z,vxi,vyi,vzi)
    
                Ntot_inj(2) = Ntot_inj(2) + 1

            END IF
            
            END DO
        END DO
	END DO


    END DO
    
    !111 CONTINUE     
   
    IF(mycoord(1) .EQ. 1) THEN
        CALL create_particle(1,DBLE(nx/2),DBLE(ny/2),DBLE(nz/2),0.2d0,0.d0,0.d0)
        CALL create_particle(2,DBLE(nx/2),DBLE(ny/2),DBLE(nz/2),0.d0,0.d0,0.d0)
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

	

END MODULE init_mod