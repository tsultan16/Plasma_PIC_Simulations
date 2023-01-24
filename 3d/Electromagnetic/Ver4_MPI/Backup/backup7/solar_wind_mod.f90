MODULE solar_wind_mod

USE constants_mod
USE data_mod
USE particleMover_mod

CONTAINS


SUBROUTINE add_ring_current()

    INTEGER :: n, i, j
    REAL*8 ::  curr, s(-1:1,-1:1,2)
    
    ! radius of current loop
    Rc = 1.0

   
    ! sinusoidal riung current 
    curr = o !*SIN(twopi * i_sim/ 50.d0) 
    
    ! decrement smoothed ring current from electric field
    GO TO 123
    IF(current_filter_on) THEN
        DO n = 1, 27
            Ex(ms(n)+ie,je,ke) = Ex(ms(n)+ie,je,ke) - (sm(n)*curr)
            Ex(ms(n)+ie,je+1,ke) = Ex(ms(n)+ie,je+1,ke) + (sm(n)*curr)
            Ey(ms(n)+ie,je,ke) = Ey(ms(n)+ie,je,ke) + (sm(n)*curr)
            Ey(ms(n)+ie+1,je,ke) = Ey(ms(n)+ie+1,je,ke) - (sm(n)*curr)
        END DO
   ELSE
        n = 14
        Ex(ms(n)+ie,je,ke) = Ex(ms(n)+ie,je,ke) - 8.d0*(sm(n)*curr)
        Ex(ms(n)+ie,je+1,ke) = Ex(ms(n)+ie,je+1,ke) + 8.d0*(sm(n)*curr)
        Ey(ms(n)+ie,je,ke) = Ey(ms(n)+ie,je,ke) + 8.d0*(sm(n)*curr)
        Ey(ms(n)+ie+1,je,ke) = Ey(ms(n)+ie+1,je,ke) - 8.d0*(sm(n)*curr)
    END IF         
    123 CONTINUE
    
    
    ! loop on xy plane
    !Jx(ie,je,ke) = o
    !Jy(ie,je,ke) = -o
    !Jy(ie+1,je,ke) = o
    !Jx(ie,je+1,ke) = -o

    ! loop on xz plane
    Jx(ie,je,ke) = o
    Jz(ie,je,ke) = -o
    Jz(ie+1,je,ke) = o
    Jx(ie,je,ke+1) = -o


    s(:,-1,1) = (/ 1., 1., 0. /)
    s(:,0,1) = (/ 0., 0., 0. /)
    s(:,1,1) = (/ -1., -1., 0. /)
    
    s(:,-1,2) = (/ -1., 0., 1. /)
    s(:,0,2) = (/ -1., 0., 1. /)
    s(:,1,2) = (/ 0., 0., 0. /)
    
    
    DO j = -1, 1
        DO i = -1, 1
            
    !        Jx(ie+i,je+j,ke) = s(i,j,1) * o
    !        Jz(ie+i,je+j,ke) = s(i,j,2) * o
            
        END DO
    END DO

    !gradually ramp up ring current
    !IF(ABS(o1) .GT. 0.d0) THEN
    !    o2 = o2 + o3
    !    o1 = o1 + o2
    !    o = o + o1
    !END IF
    
        
    IF (o1 .GT. 0.d0) THEN
        o2 = o2+o3
        o1 = o1+o2
        o = o+o1
    END IF
        
    PRINT*,'Ring Current = ',o

    
END SUBROUTINE add_ring_current



SUBROUTINE inject_solarwind()

    REAL*8 :: x, y, z, vxe, vye, vze, vxi, vyi, vzi, p(3), p1, p2, p3
    INTEGER :: i, j
    LOGICAL, PARAMETER :: lateral_flux_On = .TRUE.
    INTEGER, PARAMETER :: de = 8, di = 16
    !*****************************************
    ! inject particle pairs through x- face
    !*****************************************            
    IF(mycoord(1) .EQ. 0) THEN
    
    DO j = 1, ny, 2     
        DO i = 1, nz, 2 
           
   
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = 1 + p(1) 
            y = j + 2*p(2) 
            z = i + 2*p(3) 
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            vxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                       
            p = p - 0.5
            vye = vth_e*(p(1) + p(2) + p(3))

            CALL RAND_NUM(p)                                       
            p = p - 0.5                
            vze = vth_e*(p(1) + p(2) + p(3)) 
               

                        
            CALL create_particle(1,x,y,z,vxe,vye,vze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
            
            
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vyi = vth_i*(p(1) + p(2) + p(3))

            CALL RAND_NUM(p)                    
            p = p - 0.5                  
            vzi = vth_i*(p(1) + p(2) + p(3))
            
            
            CALL create_particle(2,x,y,z,vxi,vyi,vzi)
    
            Ntot_inj(2) = Ntot_inj(2) + 1
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            

            
        END DO
    END DO

    END IF
    
    
    IF(lateral_flux_On) THEN
    !*********************************************************************
    ! inject particles through z- face (roughly 4 electrons for every ion)
    !*********************************************************************
    
    ! first electrons
    DO j = de, ny-de+1, de         
        DO i = de, nx-de+1, de
            
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + de*p(1)
            y = j + de*p(2) 
            z = 1 + p(3)
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            vxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                       
            p = p - 0.5
            vye = vth_e*(p(1) + p(2) + p(3))

            CALL RAND_NUM(p)                                       
            vze = vth_e*table(1+64*INT(p(1))) 

            
                        
            CALL create_particle(1,x,y,z,vxe,vye,vze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = di, ny-di+1, di 
        DO i = di, nx-di+1, di
                
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + di*p(1)
            y = j + di*p(2) 
            z = 1 + p(3) 
        
            
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vyi = vth_i*(p(1) + p(2) + p(3))

            CALL RAND_NUM(p)                    
            vzi = vth_i*table(1+64*INT(p(1))) 

            
            CALL create_particle(2,x,y,z,vxi,vyi,vzi)
    
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
    DO j = de, ny-de+1, de 
        DO i = de, nx-de+1, de

            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + de*p(1)
            y = j + de*p(2) 
            z = nz  + p(3)
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            vxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                       
            p = p - 0.5
            vye = vth_e*(p(1) + p(2) + p(3))

            CALL RAND_NUM(p)                                       
            vze = -vth_e*table(1+64*INT(p(1))) 
               
                        
            CALL create_particle(1,x,y,z,vxe,vye,vze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = di, ny-di+1, di
        DO i = di, nx-di+1, di
       
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + di*p(1)
            y = j + di*p(2) 
            z = nz + p(3) 
        
            
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vyi = vth_i*(p(1) + p(2) + p(3))

            CALL RAND_NUM(p)                    
            vzi = -vth_i*table(1+64*INT(p(1))) 
            
            CALL create_particle(2,x,y,z,vxi,vyi,vzi)
    
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
    
    IF(mycoord(2) .EQ. 0) THEN
    
    ! first electrons
    DO j = de, nz-de+1, de
        DO i = de, nx-de+1, de

            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + de*p(1)
            y = 1 + p(2)
            z = j + de*p(3) 
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            vxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                       
            vye = vth_e*table(1+64*INT(p(1))) 

            CALL RAND_NUM(p)                                       
            p = p - 0.5                
            vze = vth_e*((p(1) + p(2) + p(3))) 

                        
            CALL create_particle(1,x,y,z,vxe,vye,vze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = di, nz-di+1, di
        DO i = di, nx-di+1, di
   
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + di*p(1)
            y = 1 + p(2)
            z = j + di*p(3) 
        
            
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                      
            vyi = vth_i*table(1+64*INT(p(1))) 

            CALL RAND_NUM(p)                    
            p = p - 0.5                  
            vzi = vth_i*(p(1) + p(2) + p(3))

            CALL create_particle(2,x,y,z,vxi,vyi,vzi)
    
            Ntot_inj(2) = Ntot_inj(2) + 1
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            
                    
        END DO
    END DO    
    
    END IF
    
    !*********************************************************************
    ! inject particles through y+ face (roughly 4 electrons for every ion)
    !*********************************************************************
    
    IF(mycoord(2) .EQ. nranks_y-1) THEN

    ! first electrons
    DO j = de, nz-de+1, de 
        DO i = de, nx-de+1, de
  
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + de*p(1)
            y = ny + p(2)
            z = j + de*p(3)  
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            vxe = v0x + vth_e*(vth_e + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                       
            vye = -vth_e*table(1+64*INT(p(1))) 

            CALL RAND_NUM(p)                                       
            p = p - 0.5                
            vze = vth_e*((p(1) + p(2) + p(3))) 

                        
            CALL create_particle(1,x,y,z,vxe,vye,vze)
    
            Ntot_inj(1) = Ntot_inj(1) + 1
    
            !IF(Ntot_inj(1) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF
          
        END DO
    END DO
    
    ! now ions
    DO j = di, nz-di+1, di
        DO i = di, nx-di+1, di
    
            CALL RAND_NUM(p)                                        
            p = p - 0.5
            x = i + di*p(1)
            y = ny + p(2)
            z = j + di*p(3) 
        
            
            CALL RAND_NUM(p)                                      
            p = p - 0.5
            vxi = v0x + vth_i*(vth_i + p(1) + p(2) + p(3))
                                        
            CALL RAND_NUM(p)                                      
            vyi = -vth_i*table(1+64*INT(p(1))) 

            CALL RAND_NUM(p)                    
            p = p - 0.5                  
            vzi = vth_i*(p(1) + p(2) + p(3))

            CALL create_particle(2,x,y,z,vxi,vyi,vzi)
    
            Ntot_inj(2) = Ntot_inj(2) + 1
    
            !IF(Ntot_inj(2) .GE. Np_max) THEN
            !    PRINT*,'Particle array full...terminating program.'
            !    STOP
            !END IF            
                    
        END DO
    END DO    
    
    END IF
    
    END IF
     
    
    PRINT*,''
    PRINT*,'Np_in(electrons) = ',Np_in(1)
    PRINT*,'Np_in(ions) = ',Np_in(2)
    PRINT*,''

END SUBROUTINE inject_solarwind



END MODULE solar_wind_mod