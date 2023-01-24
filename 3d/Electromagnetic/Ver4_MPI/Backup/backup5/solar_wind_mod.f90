MODULE solar_wind_mod

USE constants_mod
USE data_mod
USE particleMover_mod

CONTAINS


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



SUBROUTINE inject_solarwind()

    REAL*8 :: x, y, z, uxe, uye, uze, uxi, uyi, uzi, p(3), p1, p2, p3, &
              gam
    INTEGER :: i, j
    
    !*****************************************
    ! inject particle pairs through x- face
    !*****************************************            
    IF(mycoord(1) .EQ. 0) THEN
    
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

    END IF
    
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

END SUBROUTINE inject_solarwind



END MODULE solar_wind_mod