MODULE field_bc_mod ! electromagnetic field boundary condition routines


USE constants_mod
USE data_mod
IMPLICIT NONE


CONTAINS


!*******************************
! Periodic boundary conditions *
!*******************************

SUBROUTINE periodic_boundary_B()

    
    buffer = Bx 
    
    IF(nranks_x .EQ. 1) THEN
        !Bx(1,:,:) = buffer(nx,:,:) 
        Bx(0,:,:) = buffer(nx,:,:) 
        Bx(-1,:,:) = buffer(nx-1,:,:)

        !Bx(nx,:,:) = buffer(1,:,:) 
        Bx(nx+1,:,:) = buffer(1,:,:) 
        Bx(nx+2,:,:) = buffer(2,:,:)
        
    END IF
                                
            
    IF(nranks_y .EQ. 1) THEN
        !Bx(:,1,:) = buffer(:,ny,:) 
        Bx(:,0,:) = buffer(:,ny,:) 
        Bx(:,-1,:) = buffer(:,ny-1,:)

        !Bx(:,ny,:) = buffer(:,1,:) 
        Bx(:,ny+1,:) = buffer(:,1,:) 
        Bx(:,ny+2,:) = buffer(:,2,:)
            
    END IF
    
    IF(nz .GT. 1) THEN
        !Bx(:,:,1) = buffer(:,:,nz) 
        Bx(:,:,0) = buffer(:,:,nz) 
        Bx(:,:,-1) = buffer(:,:,nz-1)

        !Bx(:,:,nz) = buffer(:,:,1) 
        Bx(:,:,nz+1) = buffer(:,:,1) 
        Bx(:,:,nz+2) = buffer(:,:,2)
    END IF 
       
     
    buffer = By
  
     IF(nranks_x .EQ. 1) THEN
        !By(1,:,:) = buffer(nx,:,:) 
        By(0,:,:) = buffer(nx,:,:) 
        By(-1,:,:) = buffer(nx-1,:,:)

        !By(nx,:,:) = buffer(1,:,:) 
        By(nx+1,:,:) = buffer(1,:,:) 
        By(nx+2,:,:) = buffer(2,:,:)
        
    END IF
                                
            
    IF(nranks_y .EQ. 1) THEN
        !By(:,1,:) = buffer(:,ny,:) 
        By(:,0,:) = buffer(:,ny,:) 
        By(:,-1,:) = buffer(:,ny-1,:)

        !By(:,ny,:) = buffer(:,1,:) 
        By(:,ny+1,:) = buffer(:,1,:) 
        By(:,ny+2,:) = buffer(:,2,:)
            
    END IF
    
    IF(nz .GT. 1) THEN
        !By(:,:,1) = buffer(:,:,nz) 
        By(:,:,0) = buffer(:,:,nz) 
        By(:,:,-1) = buffer(:,:,nz-1)

        !By(:,:,nz) = buffer(:,:,1) 
        By(:,:,nz+1) = buffer(:,:,1) 
        By(:,:,nz+2) = buffer(:,:,2)
    END IF
    
    buffer = Bz
  
     IF(nranks_x .EQ. 1) THEN
        !Bz(1,:,:) = buffer(nx,:,:) 
        Bz(0,:,:) = buffer(nx,:,:) 
        Bz(-1,:,:) = buffer(nx-1,:,:)

        !Bz(nx,:,:) = buffer(1,:,:) 
        Bz(nx+1,:,:) = buffer(1,:,:) 
        Bz(nx+2,:,:) = buffer(2,:,:)
        
    END IF
                                
            
    IF(nranks_y .EQ. 1) THEN
        !Bz(:,1,:) = buffer(:,ny,:) 
        Bz(:,0,:) = buffer(:,ny,:) 
        Bz(:,-1,:) = buffer(:,ny-1,:)

        !Bz(:,ny,:) = buffer(:,1,:) 
        Bz(:,ny+1,:) = buffer(:,1,:) 
        Bz(:,ny+2,:) = buffer(:,2,:)
            
    END IF
    
    IF(nz .GT. 1) THEN
        !Bz(:,:,1) = buffer(:,:,nz) 
        Bz(:,:,0) = buffer(:,:,nz) 
        Bz(:,:,-1) = buffer(:,:,nz-1)

        !Bz(:,:,nz) = buffer(:,:,1) 
        Bz(:,:,nz+1) = buffer(:,:,1) 
        Bz(:,:,nz+2) = buffer(:,:,2)
    END IF
    

END SUBROUTINE periodic_boundary_B



SUBROUTINE periodic_boundary_E()
    
     buffer = Ex 
    
    IF(nranks_x .EQ. 1) THEN
        !Ex(1,:,:) = buffer(nx,:,:) 
        Ex(0,:,:) = buffer(nx,:,:) 
        Ex(-1,:,:) = buffer(nx-1,:,:)

        !Ex(nx,:,:) = buffer(1,:,:) 
        Ex(nx+1,:,:) = buffer(1,:,:) 
        Ex(nx+2,:,:) = buffer(2,:,:)
        
    END IF
                                
            
    IF(nranks_y .EQ. 1) THEN
        !Ex(:,1,:) = buffer(:,ny,:) 
        Ex(:,0,:) = buffer(:,ny,:) 
        Ex(:,-1,:) = buffer(:,ny-1,:)

        !Ex(:,ny,:) = buffer(:,1,:) 
        Ex(:,ny+1,:) = buffer(:,1,:) 
        Ex(:,ny+2,:) = buffer(:,2,:)
            
    END IF
    
    IF(nz .GT. 1) THEN
        !Ex(:,:,1) = buffer(:,:,nz) 
        Ex(:,:,0) = buffer(:,:,nz) 
        Ex(:,:,-1) = buffer(:,:,nz-1)

        !Ex(:,:,nz) = buffer(:,:,1) 
        Ex(:,:,nz+1) = buffer(:,:,1) 
        Ex(:,:,nz+2) = buffer(:,:,2)
    END IF 
       
     
    buffer = Ey
  
     IF(nranks_x .EQ. 1) THEN
        !Ey(1,:,:) = buffer(nx,:,:) 
        Ey(0,:,:) = buffer(nx,:,:) 
        Ey(-1,:,:) = buffer(nx-1,:,:)

        !Ey(nx,:,:) = buffer(1,:,:) 
        Ey(nx+1,:,:) = buffer(1,:,:) 
        Ey(nx+2,:,:) = buffer(2,:,:)
        
    END IF
                                
            
    IF(nranks_y .EQ. 1) THEN
        !Ey(:,1,:) = buffer(:,ny,:) 
        Ey(:,0,:) = buffer(:,ny,:) 
        Ey(:,-1,:) = buffer(:,ny-1,:)

        !Ey(:,ny,:) = buffer(:,1,:) 
        Ey(:,ny+1,:) = buffer(:,1,:) 
        Ey(:,ny+2,:) = buffer(:,2,:)
            
    END IF
    
    IF(nz .GT. 1) THEN
        !Ey(:,:,1) = buffer(:,:,nz) 
        Ey(:,:,0) = buffer(:,:,nz) 
        Ey(:,:,-1) = buffer(:,:,nz-1)

        !Ey(:,:,nz) = buffer(:,:,1) 
        Ey(:,:,nz+1) = buffer(:,:,1) 
        Ey(:,:,nz+2) = buffer(:,:,2)
    END IF
    
    buffer = Ez
  
     IF(nranks_x .EQ. 1) THEN
        !Ez(1,:,:) = buffer(nx,:,:) 
        Ez(0,:,:) = buffer(nx,:,:) 
        Ez(-1,:,:) = buffer(nx-1,:,:)

        !Ez(nx,:,:) = buffer(1,:,:) 
        Ez(nx+1,:,:) = buffer(1,:,:) 
        Ez(nx+2,:,:) = buffer(2,:,:)
        
    END IF
                                
            
    IF(nranks_y .EQ. 1) THEN
        !Ez(:,1,:) = buffer(:,ny,:) 
        Ez(:,0,:) = buffer(:,ny,:) 
        Ez(:,-1,:) = buffer(:,ny-1,:)

        !Ez(:,ny,:) = buffer(:,1,:) 
        Ez(:,ny+1,:) = buffer(:,1,:) 
        Ez(:,ny+2,:) = buffer(:,2,:)
            
    END IF
    
    IF(nz .GT. 1) THEN
        !Ez(:,:,1) = buffer(:,:,nz) 
        Ez(:,:,0) = buffer(:,:,nz) 
        Ez(:,:,-1) = buffer(:,:,nz-1)

        !Ez(:,:,nz) = buffer(:,:,1) 
        Ez(:,:,nz+1) = buffer(:,:,1) 
        Ez(:,:,nz+2) = buffer(:,:,2)
    END IF
    
        
END SUBROUTINE periodic_boundary_E



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
            
                temp=.25*s*Jy(1,j,k)
                Jy(0,j,k)=Jy(0,j,k)-temp
                Jy(2,j,k)=Jy(2,j,k)-temp
                Jy(1,j,k)=Jy(1,j,k)-2.*temp
      
                temp=.25*s*Jz(1,j,k)
                Jz(0,j,k)=Jz(0,j,k)-temp
                Jz(2,j,k)=Jz(2,j,k)-temp
                Jz(1,j,k)=Jz(1,j,k)-2.*temp
   
            END DO    
        END DO
    
    END IF

    IF(mycoord(1)  .EQ. nranks_x-1) THEN

        DO j = 1, ny
            DO k = 1, nz
               
                temp=.25*s*Jy(nx+1,j,k)
                Jy(nx,j,k)=Jy(nx,j,k)-temp
                Jy(nx+2,j,k)=Jy(nx+2,j,k)-temp
                Jy(nx+1,j,k)=Jy(nx+1,j,k)-2.*temp
      
                temp=.25*s*Jz(nx+1,j,k)
                Jz(nx,j,k)=Jz(nx,j,k)-temp
                Jz(nx+2,j,k)=Jz(nx+2,j,k)-temp
                Jz(nx+1,j,k)=Jz(nx+1,j,k)-2.*temp
                
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
            
                temp=.25*s*Jz(i,1,k)
                Jz(i,0,k)=Jz(i,0,k)-temp
                Jz(i,2,k)=Jz(i,2,k)-temp
                Jz(i,1,k)=Jz(i,1,k)-2.*temp
      
                temp=.25*s*Jx(i,1,k)
                Jx(i,0,k)=Jx(i,0,k)-temp
                Jx(i,2,k)=Jx(i,2,k)-temp
                Jx(i,1,k)=Jx(i,1,k)-2.*temp
                
            END DO    
        END DO

        END IF
    
    ELSE   
        
        DO k = 1, nz
            DO i = 1, nx
            
                temp=.25*s*Jz(i,1,k)
                Jz(i,0,k)=Jz(i,0,k)-temp
                Jz(i,2,k)=Jz(i,2,k)-temp
                Jz(i,1,k)=Jz(i,1,k)-2.*temp
      
                temp=.25*s*Jx(i,1,k)
                Jx(i,0,k)=Jx(i,0,k)-temp
                Jx(i,2,k)=Jx(i,2,k)-temp
                Jx(i,1,k)=Jx(i,1,k)-2.*temp
                
            END DO    
        END DO
    
    END IF
    
    
    
    IF(ndims .EQ. 2) THEN
    
        IF(mycoord(2)  .EQ. nranks_y-1) THEN
    
        DO k = 1, nz
            DO i = 1, nx

                temp=.25*s*Jz(i,ny+1,k)
                Jz(i,ny,k)=Jz(i,ny,k)-temp
                Jz(i,ny+2,k)=Jz(i,ny+2,k)-temp
                Jz(i,ny+1,k)=Jz(i,ny+1,k)-2.*temp
                
                temp=.25*s*Jx(i,ny+1,k)
                Jx(i,ny,k)=Jx(i,ny,k)-temp
                Jx(i,ny+2,k)=Jx(i,ny+2,k)-temp
                Jx(i,ny+1,k)=Jx(i,ny+1,k)-2.*temp
                
            END DO    
        END DO

        END IF
        
    ELSE

        DO k = 1, nz
            DO i = 1, nx

                temp=.25*s*Jz(i,ny+1,k)
                Jz(i,ny,k)=Jz(i,ny,k)-temp
                Jz(i,ny+2,k)=Jz(i,ny+2,k)-temp
                Jz(i,ny+1,k)=Jz(i,ny+1,k)-2.*temp
                
                temp=.25*s*Jx(i,ny+1,k)
                Jx(i,ny,k)=Jx(i,ny,k)-temp
                Jx(i,ny+2,k)=Jx(i,ny+2,k)-temp
                Jx(i,ny+1,k)=Jx(i,ny+1,k)-2.*temp
                              
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
            
                temp=.25*s*Jx(i,j,1)
                Jx(i,j,0)=Jx(i,j,0)-temp
                Jx(i,j,2)=Jx(i,j,2)-temp
                Jx(i,j,1)=Jx(i,j,1)-2.*temp
      
                temp=.25*s*Jy(i,j,1)
                Jy(i,j,0)=Jy(i,j,0)-temp
                Jy(i,j,2)=Jy(i,j,2)-temp
                Jy(i,j,1)=Jy(i,j,1)-2.*temp
                
        END DO    
    END DO
    
    
    
    DO i = 1, nx
            DO j = 1, ny

                temp=.25*s*Jx(i,j,nz+1)
                Jx(i,j,nz)=Jx(i,j,nz)-temp
                Jx(i,j,nz+2)=Jx(i,j,nz+2)-temp
                Jx(i,j,nz+1)=Jx(i,j,nz+1)-2.*temp
                
                temp=.25*s*Jy(i,j,nz+1)
                Jy(i,j,nz)=Jy(i,j,nz)-temp
                Jy(i,j,nz+2)=Jy(i,j,nz+2)-temp
                Jy(i,j,nz+1)=Jy(i,j,nz+1)-2.*temp
                
        END DO    
    END DO
    
    
END SUBROUTINE boundary_conduct_z



END MODULE field_bc_mod