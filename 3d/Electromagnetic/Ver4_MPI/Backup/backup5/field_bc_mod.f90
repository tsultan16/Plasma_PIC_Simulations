MODULE field_bc_mod ! electromagnetic field boundary condition routines


USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


!******************************************
! Radiation absorbing boundary conditions *
!******************************************



!****************************************************************
! Lowest Order Lindman Method (a bit complicated to implement..)
!****************************************************************


! Radiation absorbing boundary for surface components of magnetic field at the x+ boundary.
! (Lowest order Lindman method.) 
SUBROUTINE radiation_b_xp()

    REAL*8 :: rs, s, os
    INTEGER :: i, j, k 
    
    ! Lindman method constants
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
    
    ! First half  update of normal component Bx at x+ boundary

    DO k = bzlow, bzhi !-1, nz+2
        DO j = bylow, byhi !-1, ny+2  
            DO i = exhi, exhi ! nx+3, nx+3    
    
                Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * ((Ey(i,j,k) - Ey(i,j,k+1)) + (Ez(i,j+1,k) - Ez(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update By surface component
    DO k = bzlow, bzhi !-1, nz+2
        DO j = eylow, byhi
            DO i = exhi, exhi ! nx+3, nx+3 

               By(i,j,k) = By(i,j,k) + rs * (By(i-1,j,k) - By(i,j,k) + s * (Bx(i,j,k) - Bx(i,j-1,k))) - &     
                    os * (Ex(i,j,k+1) - Ex(i,j,k)) - (os - c) * (Ex(i-1,j,k+1) - Ex(i-1,j,k)) - &
                    c * (Ez(i,j,k) - Ez(i-1,j,k))
   
            END DO
        END DO
    END DO
   
    ! Update Bz surface component
    DO k = ezlow, bzhi !0, nz+2
        DO j = bylow, byhi ! -1, ny+2
            DO i = exhi, exhi ! nx+3, nx+3 

                Bz(i,j,k) = Bz(i,j,k) + rs * (Bz(i-1,j,k) - Bz(i,j,k) + s * (Bx(i,j,k) - Bx(i,j,k-1))) + &
                    os * (Ex(i,j+1,k) - Ex(i,j,k)) + (os - c) *(Ex(i-1,j+1,k) - Ex(i-1,j,k)) + &
                    c * (Ey(i,j,k) - Ey(i-1,j,k))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of Bx
    DO k =  bzlow, bzhi !-1, nz+2
        DO j = bylow, byhi  
            DO i = exhi, exhi ! nx+3, nx+3     
    
                Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * ((Ey(i,j,k) - Ey(i,j,k+1)) + (Ez(i,j+1,k) - Ez(i,j,k)))

            END DO
        END DO
    END DO
   

END SUBROUTINE radiation_b_xp


! Radiation absorbing boundary for surface components of magnetic field at the y+ boundary.
! (Lowest order Lindman method.) 
SUBROUTINE radiation_b_yp()

    REAL*8 :: rs, s, os
    INTEGER :: i, j, k 
    
    ! Lindman method constants
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
    
    ! First half  update of normal component By at y+ boundary

    DO k =  bzlow, bzhi !-1, nz+2
        DO j = eyhi, eyhi !ny+3, ny+3    
            DO i = bxlow, bxhi ! -1,nx+2
    
                By(i,j,k) = By(i,j,k) - 0.5 * c * ((Ez(i,j,k) - Ez(i+1,j,k)) + (Ex(i,j,k+1) - Ex(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Bz surface component
    DO k = ezlow, bzhi !0, nz+2 
        DO j = eyhi, eyhi !ny+3, ny+3  
            DO i = bxlow, bxhi ! -1, nx+2

               Bz(i,j,k) = Bz(i,j,k) + rs * (Bz(i,j-1,k) - Bz(i,j,k) + s * (By(i,j,k) - By(i,j,k-1))) - &     
                    os * (Ey(i+1,j,k) - Ey(i,j,k)) - (os - c) * (Ey(i+1,j-1,k) - Ey(i,j-1,k)) - &
                    c * (Ex(i,j,k) - Ex(i,j-1,k))
   
            END DO
        END DO
    END DO
   
    ! Update Bx surface component
    DO k = bzlow, bzhi ! -1, nz+2 
        DO j = eyhi, eyhi !ny+3, ny+3 
            DO i = exlow, bxhi ! 0, nx+2

                Bx(i,j,k) = Bx(i,j,k) + rs * (Bx(i,j-1,k) - Bx(i,j,k) + s * (By(i,j,k) - By(i-1,j,k))) + &
                    os * (Ey(i,j,k+1) - Ey(i,j,k)) + (os - c) *(Ey(i,j-1,k+1) - Ey(i,j-1,k)) + &
                    c * (Ez(i,j,k) - Ez(i,j-1,k))
   
            END DO
        END DO
    END DO
   
    ! Final half update of By
    DO k =  bzlow, bzhi !-1, nz+2
        DO j = eyhi, eyhi !ny+3, ny+3    
            DO i = bxlow, bxhi ! -1,nx+2
    
                By(i,j,k) = By(i,j,k) - 0.5 * c * ((Ez(i,j,k) - Ez(i+1,j,k)) + (Ex(i,j,k+1) - Ex(i,j,k)))

            END DO
        END DO
    END DO
    

END SUBROUTINE radiation_b_yp




! Radiation absorbing boundary for surface components of magnetic field at the z+ boundary.
! (Lowest order Lindman method.) 
SUBROUTINE radiation_b_zp()

    REAL*8 :: rs, s, os
    INTEGER :: i, j, k 
    
    ! Lindman method constants
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
    
    ! First half  update of normal component Bz at z+ boundary
    DO k = ezhi, ezhi !nz+3, nz+3
        DO j = bylow, byhi ! -1,ny+2
            DO i = bxlow, bxhi  ! -1,nx+2
    
                     Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Bx surface component
    DO k = ezhi,ezhi !nz+3, nz+3
        DO j = bylow, byhi ! -1,ny+2
            DO i = exlow, bxhi ! 0,nx+2

                Bx(i,j,k) = Bx(i,j,k) + rs * (Bx(i,j,k-1) - Bx(i,j,k) + s * (Bz(i,j,k) - Bz(i-1,j,k))) - &     
                    os * (Ez(i,j+1,k) - Ez(i,j,k)) - (os - c) * (Ez(i,j+1,k-1) - Ez(i,j,k-1)) - &
                    c * (Ey(i,j,k) - Ey(i,j,k-1))
   
            END DO
        END DO
    END DO
   
   ! Update By surface component
   DO k =  ezhi,ezhi !nz+3, nz+3
       DO j = eylow, byhi ! 0,ny+2
           DO i = bxlow, bxhi ! -1,nx+2

                By(i,j,k) = By(i,j,k) + rs * (By(i,j,k-1) - By(i,j,k) + s * (Bz(i,j,k) - Bz(i,j-1,k))) + &
                    os * (Ez(i+1,j,k) - Ez(i,j,k)) + (os - c) *(Ez(i+1,j,k-1) - Ez(i,j,k-1)) + &
                    c * (Ex(i,j,k) - Ex(i,j,k-1))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of Bz
   DO k =  ezhi,ezhi ! nz+3, nz+3
       DO j = bylow, byhi ! -1,ny+2
           DO i = bxlow, bxhi ! -1,nx+2
    
                     Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))

            END DO
        END DO
    END DO
   

END SUBROUTINE radiation_b_zp



! Radiation absorbing boundary for surface components of electric field at the x- boundary.
! (Lowest order Lindman method.) 
SUBROUTINE radiation_e_xm()

    REAL*8 :: rs, s, os
    INTEGER :: i, j, k 
    
    ! Lindman method constants
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
    
    ! First half  update of normal component Ex at x- boundary
    DO k = ezlow, ezhi !0, nz+3
        DO j = eylow, eyhi  
            DO i = bxlow, bxlow !-1, -1

                 Ex(i,j,k) = Ex(i,j,k) - 0.5 * c * ((By(i,j,k) - By(i,j,k-1)) + (Bz(i,j-1,k) - Bz(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Ey surface component
    DO k = ezlow, ezhi ! 0, nz+3
        DO j = eylow, byhi
            DO i =  bxlow, bxlow !-1, -1

                Ey(i,j,k) = Ey(i,j,k) + rs * (Ey(i+1,j,k) - Ey(i,j,k) + s * (Ex(i,j,k) - Ex(i,j+1,k))) - &     
                    os * (Bx(i,j,k-1) - Bx(i,j,k)) - (os - c) * (Bx(i+1,j,k-1) - Bx(i+1,j,k)) - &
                    c * (Bz(i,j,k) - Bz(i+1,j,k))
   
            END DO
        END DO
    END DO
   
    ! Update Ez surface component
    DO k = ezlow, bzhi ! 0, nz+2
        DO j = eylow, eyhi
            DO i =  bxlow, bxlow !-1, -1

                Ez(i,j,k) = Ez(i,j,k) + rs * (Ez(i+1,j,k) - Ez(i,j,k) + s * (Ex(i,j,k) - Ex(i,j,k+1))) + &
                    os * (Bx(i,j-1,k) - Bx(i,j,k)) + (os - c) *(Bx(i+1,j-1,k) - Bx(i+1,j,k)) + &
                    c * (By(i,j,k) - By(i+1,j,k))
   
            END DO
        END DO
    END DO   
   
   ! Final half update of Ex
    DO k = ezlow, ezhi ! 0, nz+3
        DO j = eylow, eyhi  
            DO i = bxlow, bxlow ! -1, -1
    
                     Ex(i,j,k) = Ex(i,j,k) - 0.5 * c * ((By(i,j,k) - By(i,j,k-1)) + (Bz(i,j-1,k) - Bz(i,j,k)))

            END DO
        END DO
    END DO
    
    
   

END SUBROUTINE radiation_e_xm


! Radiation absorbing boundary for surface components of electric field at the y- boundary.
! (Lowest order Lindman method.) 
SUBROUTINE radiation_e_ym()

    REAL*8 :: rs, s, os
    INTEGER :: i, j, k 
    
    ! Lindman method constants
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
    

    ! First half  update of normal component Ey at y- boundary
    DO k = ezlow, ezhi ! 0, nz+3  
        DO j = bylow, bylow ! -1, -1
            DO i = exlow, exhi

                 Ey(i,j,k) = Ey(i,j,k) - 0.5 * c * ((Bz(i,j,k) - Bz(i-1,j,k)) + (Bx(i,j,k-1) - Bx(i,j,k)))

            END DO
        END DO
    END DO

 

    ! Update Ez surface component
    DO k = ezlow, bzhi ! 0, nz+2
        DO j = bylow, bylow ! -1, -1
            DO i = exlow, exhi

                Ez(i,j,k) = Ez(i,j,k) + rs * (Ez(i,j+1,k) - Ez(i,j,k) + s * (Ey(i,j,k) - Ey(i,j,k+1))) - &     
                    os * (By(i-1,j,k) - By(i,j,k)) - (os - c) * (By(i-1,j+1,k) - By(i,j+1,k)) - &
                    c * (Bx(i,j,k) - Bx(i,j+1,k))
   
            END DO
        END DO
    END DO


   
    ! Update Ex surface component
    DO k = ezlow, ezhi !0, nz+3
        DO j = bylow, bylow ! -1, -1
            DO i = exlow, bxhi

                Ex(i,j,k) = Ex(i,j,k) + rs * (Ex(i,j+1,k) - Ex(i,j,k) + s * (Ey(i,j,k) - Ey(i+1,j,k))) + &
                    os * (By(i,j,k-1) - By(i,j,k)) + (os - c) *(By(i,j+1,k-1) - By(i,j+1,k)) + &
                    c * (Bz(i,j,k) - Bz(i,j+1,k))
                    
            END DO
        END DO
    END DO   
   
    
   ! Final half update of Ey
    DO k = ezlow, ezhi !0, nz+3  
        DO j = bylow, bylow ! -1, -1
            DO i = exlow, exhi

                 Ey(i,j,k) = Ey(i,j,k) - 0.5 * c * ((Bz(i,j,k) - Bz(i-1,j,k)) + (Bx(i,j,k-1) - Bx(i,j,k)))

            END DO
        END DO
    END DO
       

END SUBROUTINE radiation_e_ym


! Radiation absorbing boundary for surface components of electric field at the z- boundary.
! (Lowest order Lindman method.) 
SUBROUTINE radiation_e_zm()

    REAL*8 :: rs, s, os
    INTEGER :: i, j, k 
    
    ! Lindman method constants
    rs = 2.d0 * c / (1.d0 + c)
    s = .4142136
    os = 0.5d0 * (1.d0 - s) * rs
    
    ! First half  update of normal component Ez at z- boundary
    DO k = bzlow, bzlow !-1, -1
        DO j = eylow, eyhi
            DO i = exlow, exhi  
    
                     Ez(i,j,k) = Ez(i,j,k) - 0.5 * c * ((Bx(i,j,k) - Bx(i,j-1,k)) + (By(i-1,j,k) - By(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Ex surface component
    DO k = bzlow, bzlow ! -1, -1
        DO j = eylow, eyhi
            DO i = exlow, bxhi

                Ex(i,j,k) = Ex(i,j,k) + rs * (Ex(i,j,k+1) - Ex(i,j,k) + s * (Ez(i,j,k) - Ez(i+1,j,k))) - &     
                    os * (Bz(i,j-1,k) - Bz(i,j,k)) - (os - c) * (Bz(i,j-1,k+1) - Bz(i,j,k+1)) - &
                    c * (By(i,j,k) - By(i,j,k+1))
   
            END DO
        END DO
    END DO
   
   ! Update Ey surface component
   DO k = bzlow, bzlow ! -1, -1
       DO j = eylow, byhi
           DO i = exlow, exhi

                Ey(i,j,k) = Ey(i,j,k) + rs * (Ey(i,j,k+1) - Ey(i,j,k) + s * (Ez(i,j,k) - Ez(i,j+1,k))) + &
                    os * (Bz(i-1,j,k) - Bz(i,j,k)) + (os - c) *(Bz(i-1,j,k+1) - Bz(i,j,k+1)) + &
                    c * (Bx(i,j,k) - Bx(i,j,k+1))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of Ez
    DO k = bzlow, bzlow ! -1, -1
        DO j = eylow, eyhi
            DO i = exlow, exhi  
    
                     Ez(i,j,k) = Ez(i,j,k) - 0.5 * c * ((Bx(i,j,k) - Bx(i,j-1,k)) + (By(i-1,j,k) - By(i,j,k)))

            END DO
        END DO
    END DO
    
   

END SUBROUTINE radiation_e_zm







!########################################################
! Mur (1981) Simple radiation absorbing boundary scheme
!########################################################

! For this scheme, we only need to apply boundary conditions to the tangential
! E field components at boundary surfaces of the computational domain.

!Does not work very well....

! First order Mur
SUBROUTINE rad_mur_E()


    INTEGER :: i, j, k
    REAL*8 :: s, temp 
    
    s = (c - 1.d0) / (c + 1.d0)
   

    ! 6 faces, excluding the edges     
     
    !*************  
    ! x- boundary
    !*************
    
    ! update of tangential components Ey, Ez
    DO k = 2, nz-1 
        DO j = 2, ny-1
        
            ! remaining update of interior Ex, Ey, Ez
            Ex(1,j,k) = Ex(1,j,k) + c* ((By(1,j,k-1) - By(1,j,k)) + (Bz(1,j,k) - Bz(1,j-1,k))) 

            temp = Ey(1,j,k)  ! temp variable holds old interior value
             
            ! remaining update of interior Ey
            Ey(1,j,k) = Ey(1,j,k) + c* ((Bx(1,j,k) - Bx(1,j,k-1)) + (Bz(0,j,k) - Bz(1,j,k)))
        
            Ey(0,j,k) = temp + s * (Ey(1,j,k) - Ey(0,j,k))
        
        
            temp = Ez(1,j,k)  ! temp variable holds old interior value
             
            ! remaining update of interior Ez
            Ez(1,j,k) = Ez(1,j,k) + c* ((Bx(1,j-1,k) - Bx(1,j,k)) + (By(1,j,k) - By(0,j,k))) 
        
            Ez(0,j,k) = temp + s * (Ez(1,j,k) - Ez(0,j,k))        
        

        END DO   
    END DO    
    
    
    !*************  
    ! x+ boundary
    !*************
    
     ! update of tangential components Ey, Ez
    DO k = 2, nz-1
        DO j = 2, ny-1
         
            ! remaining update of interior Ex
            Ex(nx,j,k) = Ex(nx,j,k) + c* ((By(nx,j,k-1) - By(nx,j,k)) + (Bz(nx,j,k) - Bz(nx,j-1,k))) 


            temp = Ey(nx,j,k)  ! temp variable holds old interior value
             
            ! remaining update of interior value
            Ey(nx,j,k) = Ey(nx,j,k) + c* ((Bx(nx,j,k) - Bx(nx,j,k-1)) + (Bz(nx-1,j,k) - Bz(nx,j,k)))
        
            Ey(nx+1,j,k) = temp - s * (Ey(nx+1,j,k) - Ey(nx,j,k))
            
            temp = Ez(nx,j,k)  ! temp variable holds old interior value
             
            ! remaining update of interior value
            Ez(nx,j,k) = Ez(nx,j,k) + c* ((Bx(nx,j-1,k) - Bx(nx,j,k)) + (By(nx,j,k) - By(nx-1,j,k))) 
        
            Ez(nx+1,j,k) = temp - s * (Ez(nx+1,j,k) - Ez(nx,j,k))
        
        END DO    
    END DO    
       
       
    !*************  
    ! y- boundary
    !*************
    ! update of tangential components Ex, Ez
    DO k = 2, nz-1
        DO i = 2, nx-1
        
            ! remaining update of interior Ex
            Ey(i,1,k) = Ey(i,1,k) + c* ((Bz(i-1,1,k) - Bz(i,1,k)) + (Bx(i,1,k) - Bx(i,1,k-1))) 

        
            temp = Ez(i,1,k)  ! temp variable holds old interior value
             
            ! remaining update of interior Ey
            Ez(i,1,k) = Ez(i,1,k) + c* ((By(i,1,k) - By(i-1,1,k)) + (Bx(i,0,k) - Bx(i,1,k)))
        
            Ez(i,0,k) = temp + s * (Ez(i,1,k) - Ez(i,0,k))
                     
            temp = Ex(i,1,k)  ! temp variable holds old interior value
             
            ! remaining update of interior Ez
            Ex(i,1,k) = Ex(i,1,k) + c* ((By(i,1,k-1) - By(i,1,k)) + (Bz(i,1,k) - Bz(i,0,k))) 
        
            Ex(i,0,k) = temp + s * (Ex(i,1,k) - Ex(i,0,k))        
        

        END DO    
    END DO    
    
    
    !*************  
    ! y+ boundary
    !*************
    ! update of tangential components Ex, Ez
    DO k = 2, nz-1
        DO i = 2, nx-1
        
            ! remaining update of interior Ex
            Ey(i,ny,k) = Ey(i,ny,k) + c* ((Bz(i-1,ny,k) - Bz(i,ny,k)) + (Bx(i,ny,k) - Bx(i,ny,k-1))) 
  
            temp = Ez(i,ny,k)  ! temp variable holds old interior value
             
            ! remaining update of interior Ey
            Ez(i,ny,k) = Ez(i,ny,k) + c* ((By(i,ny,k) - By(i-1,ny,k)) + (Bx(i,ny-1,k) - Bx(i,ny,k)))
        
            Ez(i,ny+1,k) = temp - s * (Ez(i,ny+1,k) - Ez(i,ny,k))
        
        
            temp = Ex(i,ny,k)  ! temp variable holds old interior value
             
            ! remaining update of interior value
            Ex(i,ny,k) = Ex(i,ny,k) + c* ((By(i,ny,k-1) - By(i,ny,k)) + (Bz(i,ny,k) - Bz(i,ny-1,k))) 
        
            Ex(i,ny+1,k) = temp - s * (Ex(i,ny+1,k) - Ex(i,ny,k))
        

        END DO    
    END DO    


   
    !*************  
    ! z- boundary
    !*************
    ! update of tangential components Ex, Ey
    DO j = 2, ny-1
        DO i = 2, nx-1
        
            ! remaining update of interior Ez
            Ez(i,j,1) = Ez(i,j,1) + c* ((Bx(i,j-1,1) - Bx(i,j,1)) + (By(i,j,1) - By(i-1,j,1))) 

        
            temp = Ex(i,j,1)  ! temp variable holds old interior value
             
            ! remaining update of interior Ex
            Ex(i,j,1) = Ex(i,j,1) + c* ((Bz(i,j,1) - Bz(i,j-1,1)) + (By(i,j,0) - By(i,j,1)))
        
            Ex(i,j,0) = temp + s * (Ex(i,j,1) - Ex(i,j,0))
                     
            temp = Ey(i,j,1)  ! temp variable holds old interior value
             
            ! remaining update of interior Ey
            Ey(i,j,1) = Ey(i,j,1) + c* ((Bz(i-1,j,1) - Bz(i,j,1)) + (Bx(i,j,1) - Bx(i,j,0))) 
        
            Ey(i,j,0) = temp + s * (Ey(i,j,1) - Ey(i,j,0))        
        

        END DO    
    END DO    


    !*************  
    ! z+ boundary
    !*************
    ! update of tangential components Ex, Ey
    DO j = 2, ny-1
        DO i = 2, nx-1
        
            ! remaining update of interior Ez
            Ez(i,j,nz) = Ez(i,j,nz) + c* ((Bx(i,j-1,nz) - Bx(i,j,nz)) + (By(i,j,nz) - By(i-1,j,nz))) 

        
            temp = Ex(i,j,nz)  ! temp variable holds old interior value
             
            ! remaining update of interior Ex
            Ex(i,j,nz) = Ex(i,j,nz) + c* ((Bz(i,j,nz) - Bz(i,j-1,nz)) + (By(i,j,nz-1) - By(i,j,nz)))
                            
            Ex(i,j,nz+1) = temp - s * (Ex(i,j,nz+1) - Ex(i,j,nz))

                    
            temp = Ey(i,j,nz)  ! temp variable holds old interior value
             
            ! remaining update of interior Ey
            Ey(i,j,nz) = Ey(i,j,nz) + c* ((Bz(i-1,j,nz) - Bz(i,j,nz)) + (Bx(i,j,nz) - Bx(i,j,nz-1))) 
        
            Ey(i,j,nz+1) = temp - s * (Ey(i,j,nz+1) - Ey(i,j,nz))        
        

        END DO    
    END DO    



   ! Now update the 12 edges...
   
   
   


END SUBROUTINE rad_mur_E






END MODULE field_bc_mod