MODULE field_bc_mod ! electromagnetic field boundary condition routines


USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


!******************************************
! Radiation absorbing boundary conditions *
!******************************************




!########################################################
! Mur (1981) Simple radiation absorbing boundary scheme
!########################################################

! For this scheme, we only need to apply boundary conditions to the tangential
! E field components at boundary surfaces of the computational domain.



!###############################
! Lindman (somewhat complicated)
!###############################



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

    DO k = -1, nz+2
        DO j = bylow, byhi  
            DO i = nx+3, nx+3    
    
                Bx(i,j,k) = Bx(i,j,k) - 0.5 * c * ((Ey(i,j,k) - Ey(i,j,k+1)) + (Ez(i,j+1,k) - Ez(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update By surface component
    DO k = -1, nz+2
        DO j = eylow, byhi
            DO i = nx+3, nx+3

               By(i,j,k) = By(i,j,k) + rs * (By(i-1,j,k) - By(i,j,k) + s * (Bx(i,j,k) - Bx(i,j-1,k))) - &     
                    os * (Ex(i,j,k+1) - Ex(i,j,k)) - (os - c) * (Ex(i-1,j,k+1) - Ex(i-1,j,k)) - &
                    c * (Ez(i,j,k) - Ez(i-1,j,k))
   
            END DO
        END DO
    END DO
   
    ! Update Bz surface component
    DO k = 0, nz+2
        DO j = bylow, byhi
            DO i = nx+3, nx+3

                Bz(i,j,k) = Bz(i,j,k) + rs * (Bz(i-1,j,k) - Bz(i,j,k) + s * (Bx(i,j,k) - Bx(i,j,k-1))) + &
                    os * (Ex(i,j+1,k) - Ex(i,j,k)) + (os - c) *(Ex(i-1,j+1,k) - Ex(i-1,j,k)) + &
                    c * (Ey(i,j,k) - Ey(i-1,j,k))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of Bx
    DO k = -1, nz+2
        DO j = bylow, byhi  
            DO i = nx+3, nx+3    
    
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

    DO k = -1, nz+2  
        DO j = ny+3, ny+3    
            DO i = bxlow, bxhi
    
                By(i,j,k) = By(i,j,k) - 0.5 * c * ((Ez(i,j,k) - Ez(i+1,j,k)) + (Ex(i,j,k+1) - Ex(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Bz surface component
    DO k = 0, nz+2
        DO j = ny+3, ny+3
            DO i = bxlow, bxhi

               Bz(i,j,k) = Bz(i,j,k) + rs * (Bz(i,j-1,k) - Bz(i,j,k) + s * (By(i,j,k) - By(i,j,k-1))) - &     
                    os * (Ey(i+1,j,k) - Ey(i,j,k)) - (os - c) * (Ey(i+1,j-1,k) - Ey(i,j-1,k)) - &
                    c * (Ex(i,j,k) - Ex(i,j-1,k))
   
            END DO
        END DO
    END DO
   
    ! Update Bx surface component
    DO k = -1, nz+2
        DO j = ny+3, ny+3
            DO i = exlow, bxhi

                Bx(i,j,k) = Bx(i,j,k) + rs * (Bx(i,j-1,k) - Bx(i,j,k) + s * (By(i,j,k) - By(i-1,j,k))) + &
                    os * (Ey(i,j,k+1) - Ey(i,j,k)) + (os - c) *(Ey(i,j-1,k+1) - Ey(i,j-1,k)) + &
                    c * (Ez(i,j,k) - Ez(i,j-1,k))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of By
    DO k = 0, nz+2
        DO j = ny+3, ny+3
            DO i = bxlow, bxhi  
    
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
    DO k = nz+3, nz+3
        DO j = bylow, byhi
            DO i = bxlow, bxhi  
    
                     Bz(i,j,k) = Bz(i,j,k) - 0.5 * c * ((Ex(i,j,k) - Ex(i,j+1,k)) + (Ey(i+1,j,k) - Ey(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Bx surface component
    DO k = nz+3, nz+3
        DO j = bylow, byhi
            DO i = exlow, bxhi

                Bx(i,j,k) = Bx(i,j,k) + rs * (Bx(i,j,k-1) - Bx(i,j,k) + s * (Bz(i,j,k) - Bz(i-1,j,k))) - &     
                    os * (Ez(i,j+1,k) - Ez(i,j,k)) - (os - c) * (Ez(i,j+1,k-1) - Ez(i,j,k-1)) - &
                    c * (Ey(i,j,k) - Ey(i,j,k-1))
   
            END DO
        END DO
    END DO
   
   ! Update By surface component
   DO k = nz+3, nz+3
       DO j = eylow, byhi
           DO i = bxlow, bxhi

                By(i,j,k) = By(i,j,k) + rs * (By(i,j,k-1) - By(i,j,k) + s * (Bz(i,j,k) - Bz(i,j-1,k))) + &
                    os * (Ez(i+1,j,k) - Ez(i,j,k)) + (os - c) *(Ez(i+1,j,k-1) - Ez(i,j,k-1)) + &
                    c * (Ex(i,j,k) - Ex(i,j,k-1))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of Bz
   DO k = nz+3, nz+3
       DO j = bylow, byhi
           DO i = bxlow, bxhi
    
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
    DO k = 0, nz+3
        DO j = eylow, eyhi  
            DO i = -1, -1

                 Ex(i,j,k) = Ex(i,j,k) - 0.5 * c * ((By(i,j,k) - By(i,j,k-1)) + (Bz(i,j-1,k) - Bz(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Ey surface component
    DO k = 0, nz+3
        DO j = eylow, byhi
            DO i = -1, -1

                Ey(i,j,k) = Ey(i,j,k) + rs * (Ey(i+1,j,k) - Ey(i,j,k) + s * (Ex(i,j,k) - Ex(i,j+1,k))) - &     
                    os * (Bx(i,j,k-1) - Bx(i,j,k)) - (os - c) * (Bx(i+1,j,k-1) - Bx(i+1,j,k)) - &
                    c * (Bz(i,j,k) - Bz(i+1,j,k))
   
            END DO
        END DO
    END DO
   
    ! Update Ez surface component
    DO k = 0, nz+2
        DO j = eylow, eyhi
            DO i = -1, -1

                Ez(i,j,k) = Ez(i,j,k) + rs * (Ez(i+1,j,k) - Ez(i,j,k) + s * (Ex(i,j,k) - Ex(i,j,k+1))) + &
                    os * (Bx(i,j-1,k) - Bx(i,j,k)) + (os - c) *(Bx(i+1,j-1,k) - Bx(i+1,j,k)) + &
                    c * (By(i,j,k) - By(i+1,j,k))
   
            END DO
        END DO
    END DO   
   
   ! Final half update of Ex
    DO k = 0, nz+3
        DO j = eylow, eyhi  
            DO i = -1, -1
    
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
    DO k = 0, nz+3  
        DO j = -1, -1
            DO i = exlow, exhi

                 Ey(i,j,k) = Ey(i,j,k) - 0.5 * c * ((Bz(i,j,k) - Bz(i-1,j,k)) + (Bx(i,j,k-1) - Bx(i,j,k)))

            END DO
        END DO
    END DO

 

    ! Update Ez surface component
    DO k = 0, nz+2
        DO j = -1, -1
            DO i = exlow, exhi

                Ez(i,j,k) = Ez(i,j,k) + rs * (Ez(i,j+1,k) - Ez(i,j,k) + s * (Ey(i,j,k) - Ey(i,j,k+1))) - &     
                    os * (By(i-1,j,k) - By(i,j,k)) - (os - c) * (By(i-1,j+1,k) - By(i,j+1,k)) - &
                    c * (Bx(i,j,k) - Bx(i,j+1,k))
   
            END DO
        END DO
    END DO


   
    ! Update Ex surface component
    DO k = 0, nz+3
        DO j = -1, -1
            DO i = exlow, bxhi

                Ex(i,j,k) = Ex(i,j,k) + rs * (Ex(i,j+1,k) - Ex(i,j,k) + s * (Ey(i,j,k) - Ey(i+1,j,k))) + &
                    os * (By(i,j,k-1) - By(i,j,k)) + (os - c) *(By(i,j+1,k-1) - By(i,j+1,k)) + &
                    c * (Bz(i,j,k) - Bz(i,j+1,k))
                    
            END DO
        END DO
    END DO   
   
    
   ! Final half update of Ey
    DO k = 0, nz+3  
        DO j = -1, -1
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
    DO k = -1, -1
        DO j = eylow, eyhi
            DO i = exlow, exhi  
    
                     Ez(i,j,k) = Ez(i,j,k) - 0.5 * c * ((Bx(i,j,k) - Bx(i,j-1,k)) + (By(i-1,j,k) - By(i,j,k)))

            END DO
        END DO
    END DO
    
    ! Update Ex surface component
    DO k = -1, -1
        DO j = eylow, eyhi
            DO i = exlow, bxhi

                Ex(i,j,k) = Ex(i,j,k) + rs * (Ex(i,j,k+1) - Ex(i,j,k) + s * (Ez(i,j,k) - Ez(i+1,j,k))) - &     
                    os * (Bz(i,j-1,k) - Bz(i,j,k)) - (os - c) * (Bz(i,j-1,k+1) - Bz(i,j,k+1)) - &
                    c * (By(i,j,k) - By(i,j,k+1))
   
            END DO
        END DO
    END DO
   
   ! Update Ey surface component
   DO k = -1, -1
       DO j = eylow, byhi
           DO i = exlow, exhi

                Ey(i,j,k) = Ey(i,j,k) + rs * (Ey(i,j,k+1) - Ey(i,j,k) + s * (Ez(i,j,k) - Ez(i,j+1,k))) + &
                    os * (Bz(i-1,j,k) - Bz(i,j,k)) + (os - c) *(Bz(i-1,j,k+1) - Bz(i,j,k+1)) + &
                    c * (Bx(i,j,k) - Bx(i,j,k+1))
   
            END DO
        END DO
    END DO
   
   
   ! Final half update of Ez
    DO k = -1, -1
        DO j = eylow, eyhi
            DO i = eylow, eyhi  
    
                     Ez(i,j,k) = Ez(i,j,k) - 0.5 * c * ((Bx(i,j,k) - Bx(i,j-1,k)) + (By(i-1,j,k) - By(i,j,k)))

            END DO
        END DO
    END DO
    
   

END SUBROUTINE radiation_e_zm




!##############################################################################
!! EDGE ROUTINES BROKEN, NEED TO FIX. MOST LIKELY SOMETHING WRONG IN edge_prel.
!******************************************************************************


! Radiation absorbing x+ edges 
SUBROUTINE radiation_b_xedge_prel()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    
    
    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO k = nz+3, nz+3
        DO j = 0, ny+2
            DO i = nx+3, nx+3

                by(i,j,k) = by(i-1,j,k-1) + (1.d0-4.d0*t)*by(i,j,k) + (1.d0-2.d0*t)*(by(i,j,k-1) + by(i-1,j,k)) + &
                            s * t * (bz(i,j,k) - bz(i,j-1,k) + bz(i-1,j,k) - bz(i-1,j-1,k) + &
                            bx(i,j,k) - bx(i,j-1,k) + bx(i,j,k-1) - bx(i-1,j,k-1))
                
            END DO
        END DO    
    END DO
      
       

    DO k = -1, nz+2
        DO j = -1, -1
            DO i = nx+3, nx+3

                by(i,j,k) = (1.d0-c*r)*(by(i,j,k)+by(i,j+1,k))+by(i-1,j,k)+by(i-1,j+1,k)-r*(by(i,j,k)+   &
                c*((1.d0-s)*(ey(i,j,k+1)-ey(i,j,k))+(1.+s)*0.25d0*(ex(i,j,k+1)-ex(i,j,k)+ex(i,j+1,k+1) - &
                ex(i,j+1,k)+ex(i-1,j,k+1)-ex(i-1,j,k)+ex(i-1,j+1,k+1)-ex(i-1,j+1,k))))

            END DO
        END DO    
    END DO
    
    
    DO k = nz+3, nz+3
        DO j = -1, -1
            DO i = -1, nx+2
    
                by(i,j,k) = (1.d0 - c*r) * (by(i,j,k) + by(i,j+1,k)) + by(i,j,k-1) + by(i,j+1,k-1) - & 
                r * (bz(i,j,k) - c*((1.d0-s) * (ey(i+1,j,k) - ey(i,j,k)) + &
                (1.d0+s)*0.25d0 * (ez(i+1,j,k) - ez(i,j,k) + ez(i+1,j+1,k) - &
                ez(i,j+1,k) + ez(i+1,j,k-1) - ez(i,j,k-1) + ez(i+1,j+1,k-1) - ez(i,j+1,k-1))))

            END DO
        END DO    
    END DO
    
    DO k = nz+3, nz+3
        DO j = -1, ny+2 
            DO i = nx+3, nx+3

                temp = bx(i,j,k) - 0.5d0 * c * (1.d0 - s) * (ez(i,j+1,k) - ez(i,j,k) + ez(i,j+1,k-1) - ez(i,j,k-1))
                bx(i,j,k) = bx(i,j,k-1) - bx(i,j,k) + p * temp + q * bz(i,j,k)
                bz(i,j,k) = bz(i-1,j,k) - bz(i,j,k) + p * bz(i,j,k) + q * temp
   
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_b_xedge_prel


SUBROUTINE radiation_b_xedge_post()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    
    
    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO k = nz+3, nz+3
        DO j = 0, ny+2
            DO i = nx+3, nx+3

                by(i,j,k) = by(i,j,k) - (1.d0-4.d0*t)*by(i-1,j,k-1) - (1.d0-2.d0*t)*(by(i,j,k-1)+by(i-1,j,k)) + &
                            s * t * (bz(i,j,k) - bz(i,j-1,k) + bz(i-1,j,k) - bz(i-1,j-1,k) + &
                            bx(i,j,k) - bx(i,j-1,k) + bx(i,j,k-1) - bx(i,j-1,k-1))
                      
           
            END DO
        END DO    
    END DO
      
       

    DO k = -1, nz+2
        DO j = -1, -1
            DO i = nx+3, nx+3

                by(i,j,k) = by(i,j,k) - by(i,j+1,k) - (1.d0 - c*r) * (by(i-1,j,k) + by(i-1,j+1,k)) + r * bx(i,j,k)
              
            END DO
        END DO    
    END DO
    
    
    DO k = nz+3, nz+3
        DO j = -1, -1
            DO i = -1, nx+2
    
                by(i,j,k) = by(i,j,k) - by(i,j+1,k) - (1.d0 - c*r) * (by(i,j,k-1) + by(i,j+1,k-1)) + r * bz(i,j,k)            

            END DO
        END DO    
    END DO
    
    DO k = nz+3, nz+3
        DO j = -1, ny+2 
            DO i = nx+3, nx+3

                temp = bz(i-1,j,k)-.5*c*(1.-s)*(ex(i,j+1,k)-ex(i,j,k)+ex(i-1,j+1,k)-ex(i-1,j,k))
                bx(i,j,k) = bx(i,j,k-1) + bx(i,j,k) - q * temp - p * bx(i,j,k-1)
                bz(i,j,k) = bz(i-1,j,k) + bz(i,j,k) - q * bx(i,j,k-1)- p * temp
               
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_b_xedge_post



! Radiation absorbing y+ edges
SUBROUTINE radiation_b_yedge_prel()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k

    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO k = 0, nz+2
        DO j = ny+3, ny+3
            DO i = nx+3, nx+3

                bz(i,j,k) = bz(i-1,j-1,k) + (1.d0-4.d0*t)*bz(i,j,k) + (1.d0-2.d0*t)*(bz(i-1,j,k) + bz(i,j-1,k)) + &
                            s * t * (bx(i,j,k) - bx(i,j,k-1) + bx(i,j-1,k) - bx(i,j-1,k-1) + &
                            by(i,j,k) - by(i,j,k-1) + by(i-1,j,k) - by(i-1,j-1,k))
                
            END DO
        END DO    
    END DO
             

    DO k = -1, -1
        DO j = ny+3, ny+3
            DO i = -1, nx+2

                bz(i,j,k) = (1.d0-c*r)*(bz(i,j,k)+bz(i,j,k+1))+bz(i,j-1,k)+bz(i,j-1,k+1)-r*(bz(i,j,k)+   &
                c*((1.d0-s)*(ez(i+1,j,k)-ez(i,j,k))+(1.+s)*0.25d0*(ey(i+1,j,k)-ey(i,j,k)+ey(i+1,j,k+1) - &
                ey(i,j,k+1)+ey(i+1,j-1,k)-ey(i,j-1,k)+ey(i+1,j-1,k+1)-ey(i,j-1,k+1))))

            END DO
        END DO    
    END DO
    
    
    DO k = -1, -1
        DO j = -1, ny+2
            DO i = nx+3, nx+3
    
                bz(i,j,k) = (1.d0 - c*r) * (bz(i,j,k) + bz(i,j,k+1)) + bz(i-1,j,k) + bz(i-1,j,k+1) - & 
                r * (bx(i,j,k) - c*((1.d0-s) * (ez(i,j+1,k) - ez(i,j,k)) + &
                (1.d0+s)*0.25d0 * (ex(i,j+1,k) - ex(i,j,k) + ex(i,j+1,k+1) - &
                ex(i,j,k+1) + ex(i-1,j+1,k) - ex(i-1,j,k) + ex(i-1,j+1,k+1) - ex(i-1,j,k+1))))

            END DO
        END DO    
    END DO
    
    
    DO k = -1, nz+2 
        DO j = ny+3, ny+3
            DO i = nx+3, nx+3

                temp = by(i,j,k) - 0.5d0 * c * (1.d0 - s) * (ex(i,j,k+1) - ex(i,j,k) + ex(i-1,j,k+1) - ex(i-1,j,k))
                by(i,j,k) = by(i-1,j,k) - by(i,j,k) + p * temp + q * bx(i,j,k)
                bx(i,j,k) = bx(i,j-1,k) - bx(i,j,k) + p * bx(i,j,k) + q * temp
   
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_b_yedge_prel



SUBROUTINE radiation_b_yedge_post()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k

    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO k = 0, nz+2
        DO j = ny+3, ny+3
            DO i = nx+3, nx+3

                bz(i,j,k) = bz(i,j,k) - (1.d0-4.d0*t)*bz(i-1,j-1,k) - (1.d0-2.d0*t)*(bz(i-1,j,k)+bz(i,j-1,k)) + &
                            s * t * (bx(i,j,k) - bx(i,j,k-1) + bx(i,j-1,k) - bx(i,j-1,k-1) + &
                            by(i,j,k) - by(i,j,k-1) + by(i-1,j,k) - by(i-1,j,k-1))               
                
            END DO
        END DO    
    END DO
             

    DO k = -1, -1
        DO j = ny+3, ny+3
            DO i = -1, nx+2

                bz(i,j,k) = bz(i,j,k) - bz(i,j,k+1) - (1.d0 - c*r) * (bz(i,j-1,k) + bz(i,j-1,k+1)) + r * by(i,j,k)
                
            END DO
        END DO    
    END DO
    
    
    DO k = -1, -1
        DO j = -1, ny+2
            DO i = nx+3, nx+3
    
                bz(i,j,k) = bz(i,j,k) - bz(i,j,k+1) - (1.d0 - c*r) * (bz(i-1,j,k) + bz(i-1,j,k+1)) + r * bx(i,j,k)            
                
            END DO
        END DO    
    END DO
    
    
    DO k = -1, nz+2 
        DO j = ny+3, ny+3
            DO i = nx+3, nx+3
             
                temp = bx(i,j-1,k)-.5*c*(1.-s)*(ey(i,j,k+1)-ey(i,j,k)+ey(i,j-1,k+1)-ey(i,j-1,k))
                       by(i,j,k) = by(i-1,j,k) + by(i,j,k) - q * temp - p * by(i-1,j,k)
                       bx(i,j,k) = bx(i,j-1,k) + bx(i,j,k) - q * by(i-1,j,k)- p * temp
               
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_b_yedge_post



! Radiation absorbing z+ edges
SUBROUTINE radiation_b_zedge_prel()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    

    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
      
    DO k= nz+3, nz+3
        DO j = ny+3, ny+3
            DO i = 0, nx+2
                bx(i,j,k) = bx(i,j-1,k-1) + (1.d0-4.d0*t)*bx(i,j,k) + (1.d0-2.d0*t)*(bx(i,j-1,k) + bx(i,j,k-1)) + &
                            s * t * (by(i,j,k) - by(i-1,j,k) + by(i,j,k-1) - by(i-1,j,k-1) + &
                            bz(i,j,k) - bz(i-1,j,k) + bz(i,j-1,k) - bz(i-1,j-1,k))
                
            END DO
        END DO    
    END DO
      
       

    DO k= nz+3, nz+3
        DO j = -1, ny+2
            DO i = -1, -1
                bx(i,j,k) = (1.d0-c*r)*(bx(i,j,k)+bx(i+1,j,k))+bx(i,j,k-1)+bx(i+1,j,k-1)-r*(bz(i,j,k)+   &
                c*((1.d0-s)*(ex(i,j+1,k)-ex(i,j,k))+(1.+s)*0.25d0*(ez(i,j+1,k)-ez(i,j,k)+ez(i+1,j+1,k) - &
                ez(i+1,j,k)+ez(i,j+1,k-1)-ez(i,j,k-1)+ez(i+1,j+1,k-1)-ez(i+1,j,k-1))))

            END DO
        END DO    
    END DO
    
    
    DO k= -1, nz+2
        DO j = ny+3, ny+3
            DO i = -1, -1
    
                bx(i,j,k) = (1.d0 - c*r) * (bx(i,j,k) + bx(i+1,j,k)) + bx(i,j-1,k) + bx(i+1,j-1,k) - & 
                r * (by(i,j,k) - c*((1.d0-s) * (ex(i,j,k+1) - ex(i,j,k)) + &
                (1.d0+s)*0.25d0 * (ey(i,j,k+1) - ey(i,j,k) + ey(i+1,j,k+1) - &
                ey(i+1,j,k) + ey(i,j-1,k+1) - ey(i,j-1,k) + ey(i+1,j-1,k+1) - ey(i+1,j-1,k))))

            END DO
        END DO    
    END DO
    
    DO k= nz+3, nz+3
        DO j = ny+3, ny+3
            DO i = -1, nx+2 
                temp = bz(i,j,k) - 0.5d0 * c * (1.d0 - s) * (ey(i+1,j,k) - ey(i,j,k) + ey(i+1,j-1,k) - ey(i,j-1,k))
                bz(i,j,k) = bz(i,j-1,k) - bz(i,j,k) + p * temp + q * by(i,j,k)
                by(i,j,k) = by(i,j,k-1) - by(i,j,k) + p * by(i,j,k) + q * temp
   
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_b_zedge_prel


SUBROUTINE radiation_b_zedge_post()

    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i,j,k

      
    s = 0.4142136
    p = (1.d0 + c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    t = c / (2.d0 * c + 1.d0 + s)
    r = 4.d0 / (2.d0 * c + 2.d0 + s)

      
    DO k= nz+3, nz+3
        DO j = ny+3, ny+3
            DO i = 0, nx+2

                bx(i,j,k) = bx(i,j,k) - (1.d0-4.d0*t)*bx(i,j-1,k-1) - (1.d0-2.d0*t)*(bx(i,j-1,k)+bx(i,j,k-1)) + &
                        s * t * (by(i,j,k) - by(i-1,j,k) + by(i,j,k-1) - by(i-1,j,k-1) + &
                        bz(i,j,k) - bz(i-1,j,k) + bz(i,j-1,k) - bz(i-1,j-1,k))
                      
            END DO
        END DO    
    END DO
      
       

    DO k= nz+3, nz+3
        DO j = -1, ny+2
            DO i = -1, -1
            
                bx(i,j,k) = bx(i,j,k) - bx(i+1,j,k) - (1.d0 - c*r) * (bx(i,j,k-1) + bx(i+1,j,k-1)) + r * bz(i,j,k)

            END DO
        END DO    
    END DO
    
    
    DO k= -1, nz+2
        DO j = ny+3, ny+3
            DO i = -1, -1
    
                 bx(i,j,k) = bx(i,j,k) - bx(i+1,j,k) - (1.d0 - c*r) * (bx(i,j-1,k) + bx(i+1,j-1,k)) + r * by(i,j,k)
                   
            END DO
        END DO    
    END DO
    
    DO k= nz+3, nz+3
        DO j = ny+3, ny+3
            DO i = -1, nx+2 

                temp = by(i,j,k-1)-.5*c*(1.-s)*(ez(i+1,j,k)-ez(i,j,k)+ez(i+1,j,k-1)-ez(i,j,k-1))
                bz(i,j,k) = bz(i,j-1,k) + bz(i,j,k) - q * temp - p * bz(i,j-1,k)
                by(i,j,k) = by(i,j,k-1) + by(i,j,k) - q * bz(i,j-1,k)- p * temp

           END DO
        END DO    
    END DO
     
END SUBROUTINE radiation_b_zedge_post



! Radiation absorbing x- edges 
SUBROUTINE radiation_e_xedge_prel()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    
    
    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO k = -1, -1
        DO j = 0, ny+2
            DO i = -1, -1
    
                ey(i,j,k) = ey(i+1,j,k+1) + (1.d0-4.d0*t)*ey(i,j,k) + (1.d0-2.d0*t)*(ey(i,j,k+1) + ey(i+1,j,k)) + &
                            s * t * (ez(i,j,k) - ez(i,j+1,k) + ez(i+1,j,k) - ez(i+1,j+1,k) + &
                            ex(i,j,k) - ex(i,j+1,k) + ex(i,j,k+1) - ex(i+1,j,k+1))
                
            END DO
        END DO    
    END DO
      
       

    DO k = 0, nz+3
        DO j = ny+3, ny+3
            DO i = -1, -1
            
                ey(i,j,k) = (1.d0-c*r)*(ey(i,j,k)+ey(i,j-1,k))+ey(i+1,j,k)+ey(i+1,j-1,k)-r*(ey(i,j,k)+   &
                c*((1.d0-s)*(by(i,j,k-1)-by(i,j,k))+(1.+s)*0.25d0*(bx(i,j,k-1)-bx(i,j,k)+bx(i,j-1,k-1) - &
                bx(i,j-1,k)+bx(i+1,j,k-1)-bx(i+1,j,k)+bx(i+1,j-1,k-1)-bx(i+1,j-1,k))))

            END DO
        END DO    
    END DO
    
    
    DO k = -1, -1
        DO j = ny+3, ny+3
            DO i= 0, nx+3       
    
                ey(i,j,k) = (1.d0 - c*r) * (ey(i,j,k) + ey(i,j,k+1)) + ey(i,j,k+1) + ey(i,j-1,k+1) - & 
                r * (ez(i,j,k) - c*((1.d0-s) * (by(i-1,j,k) - by(i,j,k)) + &
                (1.d0+s)*0.25d0 * (bz(i-1,j,k) - bz(i,j,k) + bz(i-1,j-1,k) - &
                bz(i,j-1,k) + bz(i-1,j,k+1) - bz(i,j,k+1) + bz(i-1,j-1,k+1) - bz(i,j-1,k+1))))

            END DO
        END DO    
    END DO
    
    DO k = -1, -1
        DO j = 0, ny+3 
            DO i= -1, -1

                temp = ex(i,j,k) - 0.5d0 * c * (1.d0 - s) * (bz(i,j-1,k) - bz(i,j,k) + bz(i,j-1,k+1) - bz(i,j,k+1))
                ex(i,j,k) = ex(i,j,k+1) - ex(i,j,k) + p * temp + q * ez(i,j,k)
                ez(i,j,k) = ez(i+1,j,k) - ez(i,j,k) + p * ez(i,j,k) + q * temp
   
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_e_xedge_prel



SUBROUTINE radiation_e_xedge_post()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    
    
    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO k = -1, -1
        DO j = 0, ny+2
            DO i = -1, -1
    
                ey(i,j,k) = ey(i,j,k) - (1.d0-4.d0*t)*ey(i+1,j,k+1) - (1.d0-2.d0*t)*(ey(i,j,k+1)+ey(i+1,j,k)) + &
                            s * t * (ez(i,j,k) - ez(i,j+1,k) + ez(i+1,j,k) - ez(i+1,j+1,k) + &
                            ex(i,j,k) - ex(i,j+1,k) + ex(i,j,k+1) - ex(i,j+1,k+1))
                
            END DO
        END DO    
    END DO
      
       

    DO k = 0, nz+3
        DO j = ny+3, ny+3
            DO i = -1, -1
            
                ey(i,j,k) = ey(i,j,k) - ey(i,j-1,k) - (1.d0 - c*r) * (ey(i+1,j,k) + ey(i+1,j-1,k)) + r * ex(i,j,k)

            END DO
        END DO    
    END DO
    
    
    DO k = -1, -1
        DO j = ny+3, ny+3
            DO i= 0, nx+3       
    
                ey(i,j,k) = ey(i,j,k) - ey(i,j,k+1) - (1.d0 - c*r) * (ey(i,j,k+1) + ey(i,j-1,k+1)) + r * ez(i,j,k)            


            END DO
        END DO    
    END DO
    
    DO k = -1, -1
        DO j = 0, ny+3 
            DO i= -1, -1

                temp = ez(i+1,j,k)-.5*c*(1.-s)*(bx(i,j-1,k)-bx(i,j,k)+bx(i+1,j-1,k)-bx(i+1,j,k))
                ex(i,j,k) = ex(i,j,k+1) + ex(i,j,k) - q * temp - p * ex(i,j,k+1)
                ez(i,j,k) = ez(i+1,j,k) + ez(i,j,k) - q * ex(i,j,k+1)- p * temp
               
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_e_xedge_post


! Radiation absorbing y- edges 
SUBROUTINE radiation_e_yedge_prel()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    
    
    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO i = -1, -1
        DO k = 0, nz+2
            DO j = -1, -1
    
                 ez(i,j,k) = ez(i+1,j+1,k) + (1.d0-4.d0*t)*ez(i,j,k) + (1.d0-2.d0*t)*(ez(i+1,j,k) + ez(i,j+1,k)) + &
                        s * t * (ex(i,j,k) - ex(i,j,k+1) + ex(i,j+1,k) - ex(i,j+1,k+1) + &
                        ey(i,j,k) - ey(i,j,k+1) + ey(i+1,j,k) - ey(i+1,j+1,k))
        
            END DO
        END DO    
    END DO
      
       

    DO i = 0, nx+3
        DO k = nz+3, nz+3
            DO j = -1, -1
            
                ez(i,j,k) = (1.d0-c*r)*(ez(i,j,k)+ez(i,j,k-1))+ez(i,j+1,k)+ez(i,j+1,k-1)-r*(ez(i,j,k)+   &
                c*((1.d0-s)*(bz(i-1,j,k)-bz(i,j,k))+(1.+s)*0.25d0*(by(i-1,j,k)-by(i,j,k)+by(i-1,j,k-1) - &
                by(i,j,k-1)+by(i-1,j+1,k)-by(i,j+1,k)+by(i-1,j+1,k-1)-by(i,j+1,k-1))))
                
            END DO
        END DO    
    END DO
    
    
    DO i = -1, -1
        DO k = nz+3, nz+3
            DO j = 0, ny+3       
    
                ez(i,j,k) = (1.d0 - c*r) * (ez(i,j,k) + ez(i,j,k-1)) + ez(i+1,j,k) + ez(i+1,j,k-1) - & 
                r * (ex(i,j,k) - c*((1.d0-s) * (bz(i,j-1,k) - bz(i,j,k)) + &
                (1.d0+s)*0.25d0 * (bx(i,j-1,k) - bx(i,j,k) + bx(i,j-1,k-1) - &
                bx(i,j,k-1) + bx(i+1,j-1,k) - bx(i+1,j,k) + bx(i+1,j-1,k-1) - bx(i+1,j,k-1))))  

            END DO
        END DO    
    END DO
    
    DO i = -1, -1
        DO k = 0, nz+3 
            DO j = -1, -1

                temp = ey(i,j,k) - 0.5d0 * c * (1.d0 - s) * (bx(i,j,k-1) - bx(i,j,k) + bx(i+1,j,k-1) - bx(i+1,j,k))
                ey(i,j,k) = ey(i+1,j,k) - ey(i,j,k) + p * temp + q * ex(i,j,k)
                ex(i,j,k) = ex(i,j+1,k) - ex(i,j,k) + p * ex(i,j,k) + q * temp
       
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_e_yedge_prel


SUBROUTINE radiation_e_yedge_post()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    
    
    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
    
     
    DO i = -1, -1
        DO k = 0, nz+2
            DO j = -1, -1
    
                ez(i,j,k) = ez(i,j,k) - (1.d0-4.d0*t)*ez(i+1,j+1,k) - (1.d0-2.d0*t)*(ez(i+1,j,k)+ez(i,j+1,k)) + &
                            s * t * (ex(i,j,k) - ex(i,j,k+1) + ex(i,j+1,k) - ex(i,j+1,k+1) + &
                            ey(i,j,k) - ey(i,j,k+1) + ey(i+1,j,k) - ey(i+1,j,k+1))               
        
            END DO
        END DO    
    END DO
      
       

    DO i = 0, nx+3
        DO k = nz+3, nz+3
            DO j = -1, -1
            
                ez(i,j,k) = ez(i,j,k) - ez(i,j,k-1) - (1.d0 - c*r) * (ez(i,j+1,k) + ez(i,j+1,k-1)) + r * ey(i,j,k)

                
            END DO
        END DO    
    END DO
    
    
    DO i = -1, -1
        DO k = nz+3, nz+3
            DO j = 0, ny+3       
    
                ez(i,j,k) = ez(i,j,k) - ez(i,j,k-1) - (1.d0 - c*r) * (ez(i+1,j,k) + ez(i+1,j,k-1)) + r * ex(i,j,k)            


            END DO
        END DO    
    END DO
    
    DO i = -1, -1
        DO k = 0, nz+3 
            DO j = -1, -1

                temp = ex(i,j+1,k)-.5*c*(1.-s)*(by(i,j,k-1)-by(i,j,k)+by(i,j+1,k-1)-by(i,j+1,k))
                       ey(i,j,k) = ey(i+1,j,k) + ey(i,j,k) - q * temp - p * ey(i+1,j,k)
                       ex(i,j,k) = ex(i,j+1,k) + ex(i,j,k) - q * ey(i+1,j,k)- p * temp
               
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_e_yedge_post


! Radiation absorbing z- edges
SUBROUTINE radiation_e_zedge_prel()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    

    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
      
    DO k = -1, -1
        DO j = -1, -1
            DO i = 0, nx+2
            
                ex(i,j,k) = ex(i,j+1,k+1) + (1.d0-4.d0*t)*ex(i,j,k) + (1.d0-2.d0*t)*(ex(i,j+1,k) + ex(i,j,k+1)) + &
                            s * t * (ey(i,j,k) - ey(i+1,j,k) + ey(i,j,k+1) - ey(i+1,j,k+1) + &
                            ez(i,j,k) - ez(i+1,j,k) + ez(i,j+1,k) - ez(i+1,j+1,k))
                
            END DO
        END DO    
    END DO
      
       

    DO k = -1, -1
        DO j = 0, ny+3
            DO i = nx+3, nx+3
            
                ex(i,j,k) = (1.d0-c*r)*(ex(i,j,k)+bx(i-1,j,k))+ex(i,j,k+1)+ex(i-1,j,k+1)-r*(ez(i,j,k)+   &
                c*((1.d0-s)*(bx(i,j-1,k)-bx(i,j,k))+(1.+s)*0.25d0*(bz(i,j-1,k)-bz(i,j,k)+bz(i-1,j-1,k) - &
                bz(i-1,j,k)+bz(i,j-1,k+1)-bz(i,j,k+1)+bz(i-1,j-1,k+1)-bz(i-1,j,k+1))))

            END DO
        END DO    
    END DO
    
    
    DO k= 0, nz+3
        DO j = -1, -1
            DO i = nx+3, nx+3
    
                ex(i,j,k) = (1.d0 - c*r) * (ex(i,j,k) + ex(i-1,j,k)) + ex(i,j+1,k) + ex(i-1,j+1,k) - & 
                r * (ey(i,j,k) - c*((1.d0-s) * (bx(i,j,k-1) - bx(i,j,k)) + &
                (1.d0+s)*0.25d0 * (by(i,j,k-1) - by(i,j,k) + by(i-1,j,k-1) - &
                by(i-1,j,k) + by(i,j+1,k-1) - by(i,j+1,k) + by(i-1,j+1,k-1) - by(i-1,j+1,k))))

            END DO
        END DO    
    END DO
    
    DO k= -1, -1
        DO j = -1, -1
            DO i = 0, nx+3 
            
                temp = ez(i,j,k) - 0.5d0 * c * (1.d0 - s) * (by(i-1,j,k) - by(i,j,k) + by(i-1,j+1,k) - by(i,j+1,k))
                ez(i,j,k) = ez(i,j+1,k) - ez(i,j,k) + p * temp + q * ey(i,j,k)
                ey(i,j,k) = ey(i,j,k+1) - ey(i,j,k) + p * ey(i,j,k) + q * temp
   
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_e_zedge_prel


SUBROUTINE radiation_e_zedge_post()


    REAL*8 :: p, q, r, s, t, temp
    INTEGER :: i, j, k
    

    s = 0.4142136
    t = c / (2.d0 * c + 1.0d0 + s)
    r = 4.d0 /(2.d0 * c + 2.d0 + s)    
    p = (1.d0+c) * 2.d0 / (1.d0 + 2.d0 * c * (1.d0 + c * s))
    q = c * s * 2.d0 / (1.d0+2.d0*c*(1.d0+c*s))
      
    DO k = -1, -1
        DO j = -1, -1
            DO i = 0, nx+2
     
                ex(i,j,k) = ex(i,j,k) - (1.d0-4.d0*t)*ex(i,j+1,k+1) - (1.d0-2.d0*t)*(ex(i,j+1,k)+ex(i,j,k+1)) + &
                        s * t * (ey(i,j,k) - ey(i+1,j,k) + ey(i,j,k+1) - ey(i+1,j,k+1) + &
                        ez(i,j,k) - ez(i+1,j,k) + ez(i,j+1,k) - ez(i+1,j+1,k))       
                
            END DO
        END DO    
    END DO
      
       

    DO k = -1, -1
        DO j = 0, ny+3
            DO i = nx+3, nx+3
            
                ex(i,j,k) = ex(i,j,k) - ex(i-1,j,k) - (1.d0 - c*r) * (ex(i,j,k+1) + ex(i-1,j,k+1)) + r * ez(i,j,k)
        
            END DO
        END DO    
    END DO
    
    
    DO k= 0, nz+3
        DO j = -1, -1
            DO i = nx+3, nx+3
    
                 ex(i,j,k) = ex(i,j,k) - ex(i-1,j,k) - (1.d0 - c*r) * (ex(i,j+1,k) + ex(i+1,j+1,k)) + r * ey(i,j,k)

            END DO
        END DO    
    END DO
    
    DO k= -1, -1
        DO j = -1, -1
            DO i = 0, nx+3 
            
       
                temp = ey(i,j,k+1)-.5*c*(1.-s)*(bz(i-1,j,k)-bz(i,j,k)+bz(i-1,j,k-1)-bz(i,j,k+1))
                ez(i,j,k) = ez(i,j+1,k) + ez(i,j,k) - q * temp - p * ez(i,j+1,k)
                ey(i,j,k) = ey(i,j,k+1) + ey(i,j,k) - q * ez(i,j+1,k)- p * temp     
   
           END DO
        END DO    
    END DO
     

END SUBROUTINE radiation_e_zedge_post


END MODULE field_bc_mod