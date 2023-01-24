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
    INTEGER, PARAMETER :: npass = 20
	REAL*8 :: p(3), x, y, z, vxe, vye, vze, vxi, vyi, vzi
    
    

    !GO TO 111
    
    DO l = 1, npass
    
    ! place particles uniformly across grid with approximately Maxwellian velocity distribution    
    DO k = 1, nz, 1 
        DO j = 1, ny, 1
            DO i = 1, nx, 1
	
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + p(1)
            y = j + p(2) 
            z = k !+ p(3) 
    
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
   
    !IF(mycoord(1) .EQ. 0) THEN
    !    CALL create_particle(1,DBLE(nx-2),DBLE(ny/2),DBLE(nz),0.2d0,0.d0,0.d0)
    !    CALL create_particle(2,DBLE(nx-2),DBLE(ny/2),DBLE(nz),0.d0,0.d0,0.d0)
    !END IF
   
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



END MODULE init_mod