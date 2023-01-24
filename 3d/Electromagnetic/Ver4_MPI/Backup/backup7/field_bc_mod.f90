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



! Routines for spreading out lateral currents at the boundary faces to mitigate charge build-up efects

SUBROUTINE boundary_conduct_x(s)

    REAL*8, INTENT(IN) :: s
    INTEGER :: i, j, k
    REAL*8 :: temp


    IF(mycoord(1)  .EQ. 0) THEN
    
        DO j = 1, ny
            DO k = 1, nz
            
                temp=.25*s*ey(1,j,k)
                ey(0,j,k)=ey(0,j,k)-temp
                ey(2,j,k)=ey(2,j,k)-temp
                ey(1,j,k)=ey(1,j,k)-2.*temp
      
                temp=.25*s*ez(1,j,k)
                ez(0,j,k)=ez(0,j,k)-temp
                ez(2,j,k)=ez(2,j,k)-temp
                ez(1,j,k)=ez(1,j,k)-2.*temp
   
            END DO    
        END DO
    
    END IF

    IF(mycoord(1)  .EQ. nranks_x-1) THEN

        DO j = 1, ny
            DO k = 1, nz
               
                temp=.25*s*ey(nx+1,j,k)
                ey(nx,j,k)=ey(nx,j,k)-temp
                ey(nx+2,j,k)=ey(nx+2,j,k)-temp
                ey(nx+1,j,k)=ey(nx+1,j,k)-2.*temp
      
                temp=.25*s*ez(nx+1,j,k)
                ez(nx,j,k)=ez(nx,j,k)-temp
                ez(nx+2,j,k)=ez(nx+2,j,k)-temp
                ez(nx+1,j,k)=ez(nx+1,j,k)-2.*temp
                
            END DO    
        END DO

    END IF
    
    
END SUBROUTINE boundary_conduct_x


SUBROUTINE boundary_conduct_y(s)

    REAL*8, INTENT(IN) :: s
    INTEGER :: i, j, k
    REAL*8 :: temp

    IF(ndims .EQ. 2) THEN
    
        IF(mycoord(2)  .EQ. 0) THEN
    
        DO k = 1, nz
            DO i = 1, nx
            
                temp=.25*s*ez(i,1,k)
                ez(i,0,k)=ez(i,0,k)-temp
                ez(i,2,k)=ez(i,2,k)-temp
                ez(i,1,k)=ez(i,1,k)-2.*temp
      
                temp=.25*s*ex(i,1,k)
                ex(i,0,k)=ex(i,0,k)-temp
                ex(i,2,k)=ex(i,2,k)-temp
                ex(i,1,k)=ex(i,1,k)-2.*temp
                
            END DO    
        END DO

        END IF
    
    ELSE   
        
        DO k = 1, nz
            DO i = 1, nx
            
                temp=.25*s*ez(i,1,k)
                ez(i,0,k)=ez(i,0,k)-temp
                ez(i,2,k)=ez(i,2,k)-temp
                ez(i,1,k)=ez(i,1,k)-2.*temp
      
                temp=.25*s*ex(i,1,k)
                ex(i,0,k)=ex(i,0,k)-temp
                ex(i,2,k)=ex(i,2,k)-temp
                ex(i,1,k)=ex(i,1,k)-2.*temp
                
            END DO    
        END DO
    
    END IF
    
    
    
    IF(ndims .EQ. 2) THEN
    
        IF(mycoord(2)  .EQ. nranks_y-1) THEN
    
        DO k = 1, nz
            DO i = 1, nx

                temp=.25*s*ez(i,ny+1,k)
                ez(i,ny,k)=ez(i,ny,k)-temp
                ez(i,ny+2,k)=ez(i,ny+2,k)-temp
                ez(i,ny+1,k)=ez(i,ny+1,k)-2.*temp
                
                temp=.25*s*ex(i,ny+1,k)
                ex(i,ny,k)=ex(i,ny,k)-temp
                ex(i,ny+2,k)=ex(i,ny+2,k)-temp
                ex(i,ny+1,k)=ex(i,ny+1,k)-2.*temp
                
            END DO    
        END DO

        END IF
        
    ELSE

        DO k = 1, nz
            DO i = 1, nx

                temp=.25*s*ez(i,ny+1,k)
                ez(i,ny,k)=ez(i,ny,k)-temp
                ez(i,ny+2,k)=ez(i,ny+2,k)-temp
                ez(i,ny+1,k)=ez(i,ny+1,k)-2.*temp
                
                temp=.25*s*ex(i,ny+1,k)
                ex(i,ny,k)=ex(i,ny,k)-temp
                ex(i,ny+2,k)=ex(i,ny+2,k)-temp
                ex(i,ny+1,k)=ex(i,ny+1,k)-2.*temp
                              
            END DO    
        END DO

    END IF

END SUBROUTINE boundary_conduct_y



SUBROUTINE boundary_conduct_z(s)

    REAL*8, INTENT(IN) :: s
    INTEGER :: i, j, k
    REAL*8 :: temp

    DO i = 1, nx
            DO j = 1, ny
            
                temp=.25*s*ex(i,j,1)
                ex(i,j,0)=ex(i,j,0)-temp
                ex(i,j,2)=ex(i,j,2)-temp
                ex(i,j,1)=ex(i,j,1)-2.*temp
      
                temp=.25*s*ey(i,j,1)
                ey(i,j,0)=ey(i,j,0)-temp
                ey(i,j,2)=ey(i,j,2)-temp
                ey(i,j,1)=ey(i,j,1)-2.*temp
                
        END DO    
    END DO
    
    
    
    DO i = 1, nx
            DO j = 1, ny

                temp=.25*s*ex(i,j,nz+1)
                ex(i,j,nz)=ex(i,j,nz)-temp
                ex(i,j,nz+2)=ex(i,j,nz+2)-temp
                ex(i,j,nz+1)=ex(i,j,nz+1)-2.*temp
                
                temp=.25*s*ey(i,j,nz+1)
                ey(i,j,nz)=ey(i,j,nz)-temp
                ey(i,j,nz+2)=ey(i,j,nz+2)-temp
                ey(i,j,nz+1)=ey(i,j,nz+1)-2.*temp
                
        END DO    
    END DO
    
    
END SUBROUTINE boundary_conduct_z



END MODULE field_bc_mod